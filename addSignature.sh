#!/usr/bin/env bash

#paths
SIGS='/home/mike/workspace/PellegriniResearch/sigdir/SIGS/'
ABBREVS='/home/mike/workspace/PellegriniResearch/sigdir/abbrevs_test.txt'
SIG_SCRIPT='/home/mike/workspace/PellegriniResearch/sigGeneration/Generate_Sig_Genes.py'
SIGGENES='/home/mike/workspace/PellegriniResearch/sigdir/SigGenesTest.txt'
MOUSEDIR='/home/mike/workspace/PellegriniResearch/sigGeneration/MouseHumanTranslation/'
FORMAT='/home/mike/workspace/PellegriniResearch/frontend/makeCheckboxes/signatureFormat.txt'
WEB_SCRIPT='/home/mike/workspace/PellegriniResearch/frontend/makeCheckboxes/generateTDFrontpage.py'
BASE_SITE='/home/mike/workspace/PellegriniResearch/frontend/makeCheckboxes/iSigPreCheckboxes.html'
OUTPUT_SITE='/home/mike/workspace/PellegriniResearch/frontend/iSigHomepage.html'

echo "Before running this script, there are two thins that need to be done"
echo "First, add the signature files to your SIGS directory located at $SIGS"
echo -e "Second, add the new signatures to the hierarchy located at $FORMAT\n\n"

#abbrevs
echo "Enter a comma-separate list of new signature names to add to $ABBREVS"
echo "If you've already added the names to $ABBREVS, just hit return"
read -p "Signature names:" -r sigNames
echo "Adding new signatures to $ABBREVS"
echo "$sigNames" | sed -n 1'p' | tr ',' '\n' | while read sigName; do
    echo -e "$sigName\t$sigName" >> $ABBREVS
done
echo -e "Done adding new signatures to $ABBREVS \n\n"

#SigGenes
echo "Generating $SIGGENES (condensed signature file)"
echo "This could take a while..."
python $SIG_SCRIPT -s $SIGS -n 1000 -m $MOUSEDIR -o $SIGGENES
echo -e "Done generating $SIGGENES \n\n"

#website
echo "Generating the new homepage at $OUTPUT_SITE using the base at $BASE_SITE and the format file $FORMAT"
python $WEB_SCRIPT -f $FORMAT -b $BASE_SITE -o $OUTPUT_SITE
echo -e "Done generating website \n\n"

echo "Done generating everything! Now files just need to be moved to server (unless this is run on server)"
echo "$ABBREVS should replace abbrevs.txt, probably in the scripts directory"
echo "$OUTPUT_SITE should replace whatever the current website is called, probably in the html directory"
echo "$SIGGENES should replace SigGenes.txt, probably in the scripts directory"
