#!/usr/bin/env python
from optparse import OptionParser
from datetime import datetime
import os
import random
import subprocess
HTML_BASE = """<!DOCTYPE HTML>
<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
        <title>Gene Heatmap</title>

        <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.8.2/jquery.min.js"></script>
        <style type="text/css">
.highcharts-tooltip>span {
    background: rgba(255,255,255,0.85);
    border: 1px solid silver;
    border-radius: 3px;
    box-shadow: 1px 1px 2px #888;
    padding: 8px;
    z-index: 2;
}
        </style> 
"""
import numpy

def removeNA(data):
    """Removes all lines of 'N/A' from the data and returns the
    genes which were removed"""
    i = 1
    NALines = []
    while i < len(data):
        if any('N' in x for x in data[i][1:]):
            NALines.append(data[i][0])
            del data[i]
        else:
            i += 1
    return NALines

def writeMatrix(data,output):
    """Writes an R-readable tab-delimited matrix from data"""
    f = open(output,'w')
    for line in data:
        for element in line[:-1]:
            f.write(element+"\t")
        assert "\n" in line[-1]
        f.write(line[-1])
    """
    for yLabel in yLabels[:-1]:
        f.write(yLabel+"\t")
    f.write(xLabels[-1]+"\n")
    for i in range(len(xLabels)):
        f.write(xLabels[i]+"\t")
        for dval in data[i][:-1]:
            f.write(str(dval)+"\t")
        f.write(str(data[i][-1])+"\n")"""
    f.close()

def readMatrix(infile):
    f = open(infile).read().split('\n')[:-1]
    f = [x.split(',') for x in f]
    f[0] = [l[1:-1] for l in f[0]]
    for i in range(1,len(f)):
        f[i][0] = f[i][0][1:-1]
    return f
    

def generateCanvas(dataFile,outFile,signature,invert,rowClustMethod,colClustMethod,seed,isClient):
    f = open(dataFile)
    for line in f:
        if "Signature: " + signature in line:
            break
    data = []
    for line in f:
        if line != '\n' and "Signature:" not in line:
            data.append(line)
        else:
            break
    f.close()
    if not data:
        print "Content-type: text/html\n\n<!doctype html><html><body>line not found</body></html>"
        exit(1)
    data = [x.split('\t') for x in data]

    skipLines = removeNA(data)
    #write to file
    if isClient:
        scratch = "/home/mike/workspace/PellegriniResearch/scripts/scratch"
    else:
        scratch = "/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_"+seed
    filename = scratch+"/"+seed+".geneHeatmapMatrix.txt"
    writeMatrix(data,filename)
    #pass to R script
    FNULL=open(os.devnull, 'w')
    rscript = "/UCSC/Pathways-Auxiliary/UCLApathways-R-3.1.1/R-3.1.1/bin/Rscript" if not isClient else 'Rscript'
    heatsigLoc = "/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/heatsigV4.R" if not isClient else '/home/mike/workspace/PellegriniResearch/scripts/heatsigV4.R'
    detailedHeatmapLoc = '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/goTeles_tissueDeconvolution_'+seed+'RDetailHeatmap.pdf' if not isClient else '/home/mike/workspace/PellegriniResearch/output/RDetailHeatmap.pdf'
    subprocess.call([rscript,heatsigLoc,filename,detailedHeatmapLoc,"none",rowClustMethod,colClustMethod,scratch+"/matrixForHM.txt"], stdout=FNULL, stderr=FNULL)
    data = readMatrix(scratch+"/matrixForHM.txt")
    xLabels = data[0][1:]
    yLabels = [data[x][0] for x in range(1,len(data))]
    htmlText = HTML_BASE
    #add the right location for scripts (use outFile != None to check if we are on server)
    if isClient:
        base_loc = '.'
    else:
        base_loc = ''
    htmlText += """
        <script src="""+'"'+base_loc+"""/heatmap_canvasV4.js">
    </script>
    </head>
    <body>
<script src=""" + '"'+base_loc+"""/HighCharts/js/highcharts.js"></script>
<script src="""+ '"'+base_loc+"""/HighCharts/js/modules/data.js"></script>
<script src="""+'"'+base_loc+"""/HighCharts/js/modules/heatmap.js"></script>
<script src="""+'"'+base_loc+"""/HighCharts/js/modules/exporting.js"></script>
<div id="container"></div> <br />
<a href="http://pathways-pellegrini.mcdb.ucla.edu/submit/img/goTeles_tissueDeconvolution_"""+seed+'RDetailHeatmap.pdf"'+""">R Heatmap Output</a><br>
<b>Signature genes not found in input</b>:<br>
"""
    maxVal = float('-inf')
    minVal = float('+inf')
    for skipLine in skipLines:
        htmlText += "&emsp;"+ skipLine + "<br>\n"
    htmlText = htmlText[:-1] + """
</body>
</html>
<pre id="csv" style="display: none">a
"""

    for j in range(1,len(data[-1])):
        for i in range(1,len(data)):
            if "N/A" not in data[i][j]:
                maxVal = max(maxVal,float(data[i][j]))
                minVal = min(minVal,float(data[i][j]))
                if not invert:
                    htmlText += str(j-1) + "," + str(i-1) + "," + data[i][j] + "\n"
                else:
                    htmlText += str(i-1) + "," + str(j-1) + "," + data[i][j] + "\n"
    htmlText += """</pre>
    <pre id = "legendLabel" style="display:none">""" + "Value" + """</pre>
    <pre id="max" style="display:none">"""+str(maxVal)+"""</pre>
    <pre id="min" style="display:none">""" + str(minVal)+"""</pre>
    <pre id="title" style="display:none">genes similar to """ + signature + """ signature</pre>
    <pre id= "xlabels" style="display:none">
    """
    for label in xLabels if not invert else yLabels:
        htmlText += label+","
    htmlText = htmlText[:-1] + """</pre>
<pre id="ylabels" style="display:none">
"""
    for label in yLabels if not invert else xLabels:
        htmlText += label.replace('"','')+','
    htmlText = htmlText[:-1] + """</pre>
</body>
</html>"""
    if not outFile:
        print "Content-type: text/html\n\n"+htmlText
    else:
        with open(outFile,'w') as out:
            out.write(htmlText)
"""
import cgi
import cgitb
cgitb.enable()
form = cgi.FieldStorage()
signature = form["signatureName"].value
scaleArg = "scale" in form
invertArg = "invert" in form
rowArg = form["row_metric"].value
colArg = form["col_metric"].value
seed = form["seed"].value
fileSource = "/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_" + seed + '/' + seed + ".matrix.txt.full.txt"
generateCanvas(fileSource,None,signature,invertArg,rowArg,colArg,seed,False)"""
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input")
    parser.add_option("-o", "--output", dest="output")
    parser.add_option("-s", "--sig", dest="sig")
    parser.add_option("-v", "--invert", dest="invert")
    parser.add_option("-r", "--row", dest="row")
    parser.add_option("-c", "--col", dest="col")
    (options, args) = parser.parse_args()
    generateCanvas(options.input,options.output,options.sig,options.invert!="none",options.row,options.col,'-1',True)

