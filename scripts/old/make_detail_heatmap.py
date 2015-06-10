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
        <script src="./heatmap_canvas.js">
    </script>
    </head>
    <body>
<script src="./HighCharts/js/highcharts.js"></script>
<script src="./HighCharts/js/modules/data.js"></script>
<script src="./HighCharts/js/modules/heatmap.js"></script>
<script src="./HighCharts/js/modules/exporting.js"></script>
<div id="container"></div> <br />
Genes with no value: 
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


def scale(data):
    """Replaces each entry of data with it's z-score"""
    dcopy = data[::]
    i = 1
    d = numpy.matrix([[float(x) for x in line[1:]] for line in dcopy[1:]])
    m = d.mean()
    s = d.std()
    for i in range(1,len(data)):
        for j in range(1,len(data[0])):
            if "N" not in data[i][j]:
                data[i][j] = str((float(data[i][j])-m)/s)

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
    

def generateCanvas(dataFile,signature,outFile,invert,normalize,bcolumn,rowClustMethod,colClustMethod):
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
        print "line not found"
        exit(1)
    data = [x.split('\t') for x in data]
    skipLines = removeNA(data)
    print data,len(data)
    #generate job ID
    randSeed = str(random.randint(0,9999))
    while len(randSeed) != 4:
        randSeed = "0" + randSeed
    seed = datetime.now().strftime('%Y%m%d%H%M%S') + randSeed
    #write to file
    scratch = "./scratch/goTeles_geneHeatmap"+seed
    filename = scratch+"/"+seed+".geneHeatmapMatrix.txt"
    os.mkdir(scratch)
    print data, "before"
    writeMatrix(data,filename)
    #pass to R script
    subprocess.call(["Rscript","heatsigV2.R",filename,scratch+"/ignore",bcolumn,rowClustMethod,colClustMethod,scratch+"/matrixForHM.txt",scratch+"/ignore",scratch+"/ignore"])
    data = readMatrix(scratch+"/matrixForHM.txt")
    xLabels = data[0][1:]
    yLabels = [data[x][0] for x in range(1,len(data))]
    print data, "after"
    if normalize:
        scale(data)
    htmlText = HTML_BASE
    maxVal = float('-inf')
    minVal = float('+inf')
    for skipLine in skipLines:
        htmlText += skipLine + ","
    htmlText = htmlText[:-1] + """
</body>
</html>
<pre id="csv" style="display: none">"""

    for j in range(1,len(data[-1])):
        for i in range(1,len(data)):
            if "N/A" not in data[i][j]:
                maxVal = max(maxVal,float(data[i][j]))
                minVal = min(minVal,float(data[i][j]))
                if not invert:
                    htmlText += str(j-1) + "," + str(i-1) + "," + data[i][j] + "\n"
                else:
                    htmlText += str(i-1) + "," + str(j-1) + "," + data[i][j] + "\n"
    if normalize:
        maxVal = max(abs(maxVal),abs(minVal))
        minVal = -maxVal
    htmlText += """</pre>
    <pre id = "legendLabel" style="display:none">""" + ("Z-score" if normalize else "") + """</pre>
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
        htmlText += label[1:-1]+","
    htmlText = htmlText[:-1] + """</pre>
</body>
</html>"""
    f_out = open(outFile,'w')
    f_out.write(htmlText)

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i","--input",dest="dataFile",help="data source")
    parser.add_option("-o","--output",dest="outFile",help="output html location")
    parser.add_option("-s","--signature",dest="signature",help="signature generating")
    parser.add_option("-b","--bcolumn",dest="bcolumn")
    parser.add_option("-x",dest="horiz")
    parser.add_option("-v",dest="vert")
    (options,args) = parser.parse_args()
    generateCanvas(options.dataFile,options.signature,options.outFile,True,False,options.bcolumn,options.horiz,options.vert)
