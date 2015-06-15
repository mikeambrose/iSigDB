#!/usr/bin/env python
"""script written by Mike Ambrose, mikeambrose@berkeley.edu"""
import os.path
import os
from optparse import OptionParser

def generateAbbrevs(abbrev):
    """Generates a dictionary of abbreviation : full name
    abbrev is a two-column, tab-delimited file of abbreviation\tfull name
    returns the dictionary of abbrev : full name
    """
    abbrevs = dict()
    with open(abbrev) as all_lines:
        for line in all_lines.read().split('\n'):
            if line:
                abbrev,full = line.split('\t')
                abbrevs[abbrev] = full
    return abbrevs

def generateAToLen(sigs):
    """Generate a dictionary of abbreviation : number of genes in signature
    sigs is all signature files (can include abbrevs.txt, but doesn't have to)
    returns the dictionary of abbrev : len(sig file)
    """
    aToLen = {}
    for sig in sigs:
        abbrev = os.path.basename(sig)
        if '--' in abbrev:
            abbrev = abbrev[:abbrev.index('--')]
        fileLength = sum(1 for line in open(sig))
        assert abbrev not in aToLen
        aToLen[abbrev] = fileLength
    del aToLen['abbrevs.txt']
    return aToLen

def generateGenesChanged(aToLen):
    """Generates the html corresponding to the genesChanged section of the js
    aToLen is a dictionary of abbreviation : number of lines in the file
        If creating this file is too time-consuming and the set of limited genes is already known,
            this can be passed in with only relevant genes
    returns the html corresponding to the whole genesChanged function
    """
    dropdownValues = [10,25,50,100,250,500,1000]
    disableAt = [[],[],[],[],[],[],[]]
    for abbrev in aToLen:
        for i in range(len(dropdownValues)):
            if aToLen[abbrev] <= dropdownValues[i]:
                disableAt[i].append(abbrev)
                break
    baseHTML = "function genesChanged(dropdown) {"
    for i in range(len(disableAt)):
        if not disableAt[i]:
            continue
        baseHTML += "if (dropdown.selectedIndex > " + str(i) + ") {\n    var disableList = " + str(disableAt[i]).replace("\"","'") + """;
    for (var i = 0; i < disableList.length; i++) {
      document.getElementById(disableList[i]).checked = false;
      document.getElementById(disableList[i]).disabled = true;
    }
  }
  else {
    var disableList = """ + str(disableAt[i]).replace('"',"'") + """;
    for (var i = 0; i < disableList.length; i++) {
      document.getElementById(disableList[i]).disabled = false;
    }
  }

"""
    return baseHTML + "\n}"
    
def generateJSTree(formatLines,abbrevs,delim="    "):
    """generates the JSTree portion of the front tree based on the file formatLines"""
    lastInd = 0 #indendation of previous line (number of tabs)
    returnHTML = ''
    for line in formatLines:
        currentInd = line.count(delim)
        if currentInd == lastInd: #no change
            returnHTML += "</li><li id=\"" + line.replace(delim,'') + "\">"+line.replace(delim,'') + '\n'
        elif currentInd > lastInd: #up one
            assert (currentInd-lastInd) == 1
            returnHTML += "<ul><li id=\"" + line.replace(delim,'') + "\">"+line.replace(delim,'') + '\n'
        else: #down some
            spacingDiff = lastInd-currentInd
            returnHTML += "</li>"+"</ul></li>"*spacingDiff + "<li id=\"" + line.replace(delim,'') + "\">"+line.replace(delim,'') + '\n'
        lastInd = currentInd
    #get rid of first </li><li> and last <li> and adding auto-select
    returnHTML = returnHTML[5:-9]
    return returnHTML

def generateHTML(abbrevs, formatLines, aToLen, output):
    """Generates the HTML frontend (currently called tissueDeconvolutionV2/3
    abbrevs is a dictioanry of abbreviation : full name 
        see generateAbbrevs
    groups is a dictionary of group name:a list of all abbreviations in that group
        see generateGroups
    output is the file which the html is written to
    returns nothing"""

    #------header and constant js functions--------

    returnHTML = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>iSigDB</title>
<link rel="stylesheet" href="../vakata-jstree-5bece58/dist/themes/default/style.min.css" />
<link rel="stylesheet" href="../tooltip.css">

<script type="text/javascript">

function toggleVisibility(section) {
    var section = document.getElementById(section);
    if (section.style.display == 'block' || section.style.display == ''){
        section.style.display = 'none';
    } else {
        section.style.display = 'block';
    }
}

function submitCheckboxes() {
    var checkedIds = $('#sigs').jstree("get_selected");
    document.getElementById('checkedSigs').value = checkedIds.join(",");
}

</script>

<br>
<form id="tissueDeconvolution" name="tissueDeconvolution" method="post" action="/cgi-bin/goTeles/tissueDeconvolutionByRank.cgi" ENCTYPE="multipart/form-data">
<strong class = "tooltip">Submit a file:</strong> <br>
<input type="radio" name="uploadSettings" value="client" checked="checked" onchange="toggleVisibility('clientFile');toggleVisibility('serverFile');">Upload a file from your computer <a href="#" class="tooltip"><img src="../questionmark.png" width=15 height=15><span><strong>File guidelines</strong><br/>The uploaded file must be tab-separated. Columns should correspond to samples and rows should correspond to genes. There is a size limit of 100mb.</span></a>  <br>
<input type="radio" name="uploadSettings" value="server" onchange="toggleVisibility('clientFile');toggleVisibility('serverFile');">Use a file from the server <br>
<div id="clientFile" name="clientFile">
File: <input type="file" name="matrix_file" size="30">
</div>
<div id="serverFile" name="serverFile" style="display:none">
Filename: <select id="serverFileName" name="serverFileName">
{FILES GO HERE}
</select>
or <a href="./serverFileUpload.html">upload a file</a>
</div>
<hr>

<strong>Heatmap options<a href="OptionDescriptions.html">[Help]</a>:</strong><br /><br />
Compute heatmap based on:&nbsp;
<select id="heatmap_metric" name="heatmap_metric">
<option value="rank_avg">Rank Average</option>
<option value="rank_delta">Rank Delta</option>
<option value="log">Log10-Transformed Values</option>
<option value="absolute">Values (Do Not Log-Transform)</option>
</select> <a href="#" class="tooltip"><img src="../questionmark.png" width=15 height=15><span><strong>Rank Average</strong> computes similarity between genes based on ordering them by importance<br><strong>Rank Delta</strong> computes similarity using the difference between a gene's rank in one sample and its mean rank across all samples<br><strong>Log10-Transformed Values</strong> uses log-transformed raw expression values instead of gene rank<br><strong>Values (Do Not Log-Transform)</strong> uses non-log-transformed raw values instead of gene rank</span></a>

<br /><br />
<input type="checkbox" name="scale_columns" id="scale_columns" value="checked"
checked> Scale heatmap <a href="#" class="tooltip"><img src="../questionmark.png" width=15 height=15><span>Scaling a heatmap replaces output values with their z-scores taken across the entire matrix</span></a>

<br /><br />

<input type="checkbox" name="invert" id="invert" value="checked" checked> Signatures on vertical axis

<br /><br />

<input type="checkbox" name="fixed" id="fixed" value="checked"> Keep colors constant
<a href="#" class="tooltip"><img src="../questionmark.png" width=15 height=15><span>If unchecked, the colors in the heatmap will rescale to fit the input data better. However, if you wish to compare different heatmaps, it is better to check this option so that their colors correspond to the same values</span></a>
<br> <br>

Metric for sample clustering: &nbsp;
<select id="row_metric" name="row_metric">
<option value="euclidean">Euclidean Distance</option>
<option value="pear_cor">Pearson Correlation</option>
<option value="none">None (Do Not Cluster)</option>
</select>
<br />

Metric for signature clustering: &nbsp;
<select id="col_metric" name="col_metric">
<option value="euclidean">Euclidean Distance</option>
<option value="pear_cor">Pearson Correlation</option>
<option value="none">None (Do Not Cluster)</option>
</select><br />


<br />

Number of genes per signature:&nbsp;<select
name="num_genes" onchange="genesChanged(this);"><option value="10">10</option><option
value="25">25</option><option value="50" selected="selected">50</option><option
value="100">100</option><option value="250">250</option><option
value="500">500</option><option value="1000">1000</option></select>
<a href="#" class="tooltip"><img src="../questionmark.png" width=15 height=15><span>How many genes from each signature will be compared against the input file.</span></a>
<hr>

<strong>Select signatures:</strong><br />
Click on a category to select specific signatures within that category <br />
<div id="sigs">
<ul>
"""
    #--------checkboxes for each group---------
    returnHTML += generateJSTree(formatLines,abbrevs)
    #--------small ending section--------------
    returnHTML += """</ul></div>
<input type="submit" id="tissueDeconvolutionSubmit" name="tissueDeconvolutionSubmit" onclick="submitCheckboxes()" value="Submit" />

</form>
<script src="../js/jquery-1.11.3.min.js"></script>
	<script src="../vakata-jstree-5bece58/dist/jstree.min.js"></script>
    <script>
        $(function () {
            $('#sigs').jstree({
                "checkbox" : {
                    "keep_selected_style" : false
                },
                "plugins" : [ "checkbox" ],
                "core": {
                    "themes":{
                        "icons":false
                    }
                }
            });
        window.onload = function (){$('#sigs').jstree('close_all');$('#sigs').jstree('check_all');};
        });

    </script>
</body>
</html>"""
    f_write = open(output,'w')
    f_write.write(returnHTML)
    f_write.close()
    print("successfully generated")

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-f","--format",dest="format",help="file with formatting")
    parser.add_option("-a","--abbrevs",dest="abbrevs",help="abbrevs file")
    parser.add_option("-o","--output",dest="output",help="output filename")
    parser.add_option("-s","--sigs",dest="sigs",help="directory of all signatures")
    options,_ = parser.parse_args()
    #abbrevs = generateAbbrevs(options.abbrevs)
    #aToLen = generateAToLen([options.sigs + name for name in os.listdir(options.sigs)])
    abbrevs,aToLen = None,None
    generateHTML(abbrevs,open(options.format).read().split('\n'),aToLen,options.output)
