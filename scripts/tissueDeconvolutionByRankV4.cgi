#!/usr/bin/perl -w

use strict;
use warnings;

use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);
use File::Basename;

use lib '/UCSC/Pathways-Auxiliary/';
require 'common_functions_algal.pl';

# limit input size to 100MB

$CGI::POST_MAX = 1024 * 100000;
my $cgi = new CGI;

ThrowInternalError('File size limit (100MB) exceeded.')
    if ( $cgi->cgi_error() );

my %raw_user_input;
for my $key ( $cgi->param() ) {
    $raw_user_input{$key} = $cgi->param($key);
}

# create unique job ID and create local work directory

my ($year, $month, $date, $hour, $min, $sec) = (localtime(time))[5, 4, 3, 2, 1, 0];
my $job_id = sprintf('%04d%02d%02d%02d%02d%02d%04d', $year + 1900, $month + 1, $date, $hour, $min, $sec, int(rand(9999)));

my $scratch_dir    = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space';
my $local_work_dir = "${scratch_dir}/goTeles_tissueDeconvolutionV2_${job_id}";

umask 000;
mkdir($local_work_dir);

# get user-submitted matrix

my $local_matrix_file_name = $job_id . '.txt';
my $local_matrix_path      = "${local_work_dir}/${local_matrix_file_name}";
my $file_name = $raw_user_input{'matrix_file'};

open my $SUBMISSION_DEST, '>', $local_matrix_path or die $!;

if ($file_name) {
    my $SUBMISSION_SRC = $cgi->upload('matrix_file');

    binmode $SUBMISSION_DEST;
    while (my $submission_src_data = <$SUBMISSION_SRC>) {
        $submission_src_data =~ s/
/\n/g;
        print $SUBMISSION_DEST $submission_src_data;
    }
}
elsif (!$file_name && $raw_user_input{'uploadSettings'} ne "server") {
    ThrowInternalError('No input provided.');
}
elsif (!$raw_user_input{'serverFileName'}) {
    ThrowInternalError('No filename provided.');
}
else {
    #replace with actual path
    open my $FILE, '<', '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/iSigDB_uploads/associations.txt' or die "Could not find associations.txt file $!";
    while (my $line = <$FILE>) {
        chomp $line;
        my @assoc_line = split(/\t/,$line);
        if ($assoc_line[0] eq $raw_user_input{'serverFileName'}) {
            #replace with actual path
            $local_matrix_path = "/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/iSigDB_uploads/$assoc_line[1]";
            unless (-e $local_matrix_path) {
                ThrowInternalError('Cannot find associated file');
            }
        }
    }
}

close $SUBMISSION_DEST;

my $heatmap_chosen = !(grep{$_ eq "matrix"} (keys %raw_user_input))
if ($heatmap_chosen) {

    my @selected_sigs = split(/,/,$raw_user_input{'checkedSigs'});

    my $local_abbrev_list_path = "${local_work_dir}/abbrevs.txt";
    open my $ABBREV_LIST, '>', $local_abbrev_list_path or die $!;


    my %full_to_abbrev;
    open my $abbrevs_file, '<', '/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/abbrevs.txt' or die "file not found$!";
    while (my $line = <$abbrevs_file>) {
        chomp $line;
        my ($abbrev, $full) = split /\t/,$line;
        $full_to_abbrev{$full} = $abbrev;
    }

    my $num_selected = 0;

    foreach my $long_name (keys %full_to_abbrev) {
        if (grep {$_ eq $long_name} @selected_sigs)) {
            print $ABBREV_LIST $full_to_abbrev{$long_name} . "\t" . $long_name . "\n";
            $num_selected++;
        }
    }

    my %acceptable_num_genes = (10 => undef,
                                25 => undef,
                                50 => undef,
                                100 => undef,
                                250 => undef,
                                500 => undef,
                                1000 => undef);

    my $num_genes = $raw_user_input{'num_genes'};

    if (!exists $acceptable_num_genes{$num_genes}) {
      ThrowInternalError("Invalid number of genes");
    }
    # launch python script


    my %acceptable_heatmap_metrics = ('rank_avg' => undef,
                                      'rank_delta' => undef,
                                      'log' => undef,
                                      'absolute' => undef);

    my $heatmap_metric_argument = $raw_user_input{'heatmap_metric'};
    if (!exists $acceptable_heatmap_metrics{$heatmap_metric_argument}) {
      ThrowInternalError("Invalid heatmap metric.");
    }

    my $log_argument = '-l';

    if ($heatmap_metric_argument eq 'absolute') {
      $heatmap_metric_argument = 'log';
      $log_argument = '-a';
    }

    my $scale_columns_argument = 'matrix';

    if ($raw_user_input{'scale_columns'} ne 'checked') {
      $scale_columns_argument = 'none';
    }

    my $invert_argument = $raw_user_input{'invert'} eq 'checked' ? 'checked' : 'none';

    my $fixed_argument = $raw_user_input{'fixed'} eq 'checked' ? 'checked' : 'none';

    my %acceptable_metrics = ('pear_cor' => undef,
                              'euclidean' => undef,
                              'none' => undef);

    my $row_metric = $raw_user_input{'row_metric'};
    my $col_metric = $raw_user_input{'col_metric'};
    my $null_metric = $raw_user_input{'null'};

    if ((!exists $acceptable_metrics{$row_metric}) || (!exists $acceptable_metrics{$col_metric})) {
      ThrowInternalError("Invalid clustering metric");
    }
    my $gene_option = 'none';
} else {
    my $heatmap_metric_argument = $raw_user_input{'matrix_metric'};
    my $invert_argument = $raw_user_input{'matrix_invert'};
    my $row_metric = $raw_user_input{'matrix_row_metric'};
    my $col_metric = $raw_user_input{'matrix_col_metric'};
    if ($raw_user_input{'spear_gene'} == 'spearGeneAll'){
        my $num_genes = 'none';
    } elsif ($raw_user_input{'spear_gene'} == 'spearGeneTop'){
        my $num_genes = $raw_user_input{'matrix_num_genes'};
    } elsif ($raw_user_input{'spear_gene'} == 'spearGeneMag'){
        my $num_genes = $raw_user_input{'matrix_mag'};
    }
    my $scale_columns_argument = 'none';
    my $log_argument = '-a';
    my $fixed_argument = 'none';
    my $gene_argument = $raw_user_input{'spear_gene'};
    my $null_argument = 'none';
}

#logging
#TODO: fix log path
my $logpath = '.';
open(my $log, '>>', $logpath) or die $!;
say $log "$ENV{'REMOTE_ADDR'}\t$job_id\t$heatmap_metric_argument\t$scale_columns_argument\t$log_argument\t$row_metric\t$col_metric\t$invert_argument\t$fixed_argument\t$gene_argument";
close $log;
my $exec_path = '/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank';
#ThrowInternalError("python $exec_path/Sig_Avg_Matrix_Derm.RankV4.py -n $num_genes -t $local_matrix_path -s $exec_path/SigGenes.txt -g ${local_work_dir}/abbrevs.txt -j ${job_id} -v $heatmap_metric_argument -z $scale_columns_argument $log_argument -r $row_metric -c $col_metric -i $invert_argument");
my $python_out = `python $exec_path/Sig_Avg_Matrix_Derm.RankV4.py -n $num_genes -t $local_matrix_path -s $exec_path/SIGS -x $exec_path/SigGenes.txt -g ${local_work_dir}/abbrevs.txt -j ${job_id} -v $heatmap_metric_argument -z $scale_columns_argument $log_argument -r $row_metric -c $col_metric -i $invert_argument -f $fixed_argument -o $gene_argument -u $null_argument --server`;

#$python_out =~ s/null device \n          1 //;
#$python_out =~ s/pdf \n  2 //;

# display output

print qq(Content-type: text/html

$python_out);


