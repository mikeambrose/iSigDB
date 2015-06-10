"""A script to generate the html file which will display the heatmap"""
from optparse import OptionParser
#a constant string used for the body of the html
HTML_BASE = """<!DOCTYPE HTML>
<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
        <title>Heatmap</title>

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
<b>Downloads: </b> <br />
<a href="TODO">Processed data</a> <br />
<a href="TODO">R heatmap (pdf)</a> <br />
<a href="TODO">Signature similarity dendogram</a> <br />
<a href="TODO">Expression similarity dendogram</a> <br /> <br />
Look in more detail at one signature:
<form id="detailedHeatmap" name="detailedHeatmap" method="post" action="./sig_heatmap.cgi" ENCTYPE="multipart/form-data">
<select>
"""


def generateCanvas(dataFile,outFile,legendLabel='',invertHeatmap=True,centerAroundZero=False,givenMinVal=None,givenMaxVal=None):
    f = open(dataFile).read().split('\n')
    f = [x.split(',') for x in f]
    while not f[-1] or not any(f[-1]):    del f[-1]
    # the first row is the x labels
    xLabels = f[0][1:]
    # the first column is the y labels
    yLabels = [f[x][0] for x in range(1,len(f))]
    #now we generate the actual html
    #todo: there's probably a way of doing this that isn't terrible string manipulation
    htmlText = HTML_BASE
    #add values for option
    for label in sorted(xLabels):
        htmlText += "<option value=" + label + ">" + label[1:-1] + "</option>\n"
    htmlText += """</select> <br/>
<input type="checkbox" value="invert"> Invert output <br/>

Metric for row (sample) clustering: &nbsp;
<select id="row_metric" name="row_metric">
<option value="euclidean">Euclidean Distance</option>
<option value="pear_cor">Pearson Correlation</option>
<option value="none">None (Do Not Cluster)</option>
</select><br />

Metric for column (signature) clustering: &nbsp;
<select id="col_metric" name="col_metric">
<option value="euclidean">Euclidean Distance</option>
<option value="pear_cor">Pearson Correlation</option>
<option value="none">None (Do Not Cluster)</option>
</select><br />
<input type="submit" value="Go">
    </body>
</form>
</html>
<pre id="csv" style="display: none">
ignored line
"""
    #next add values
    minVal = float('+inf')
    maxVal = float('-inf')
    for j in range(1,len(f[-1])):
        for i in range(1,len(f)):
            minVal = min(float(f[i][j]),minVal)
            maxVal = max(float(f[i][j]),maxVal)
            #print len(f),len(xLabels),len(f[0]),len(yLabels),f[i][j]
            if not invertHeatmap:
                htmlText += str(j-1) + "," + str(i-1) + "," + f[i][j] + "\n" 
            else:
                htmlText += str(i-1) + "," + str(j-1) + "," + f[i][j] + "\n"
        htmlText = htmlText[:-1] + '\n'
    if centerAroundZero:
        maxVal = max(abs(maxVal),abs(minVal))
        minVal = -maxVal
    minVal = givenMinVal if (givenMinVal != None) else minVal
    maxVal = givenMaxVal if (givenMaxVal != None) else maxVal
    htmlText += """</pre>
<pre id="legendLabel" style="display:none">"""+legendLabel
    htmlText += """</pre>
<pre id="title" style="display:none">Heatmap</pre>
<pre id="max" style="display:none">""" + str(maxVal)
    htmlText += """</pre>
<pre id="min" style="display:none">""" + str(minVal)
    htmlText += """</pre>
<pre id="xlabels" style="display: none">
"""
    for label in xLabels if not invertHeatmap else yLabels:
        htmlText += label[1:-1]+","
    htmlText = htmlText[:-1] + """</pre>
<pre id="ylabels" style="display: none">
"""
    for label in yLabels if not invertHeatmap else xLabels:
        htmlText += label[1:-1]+","
    htmlText = htmlText[:-1] + """</pre>
    </body>
</html>"""
    f_out = open(outFile,'w')
    f_out.write(htmlText)

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i","--input",dest="inputFile",help="data source")
    parser.add_option("-o","--output",dest="outputFile",help="output html location")
    (options,args) = parser.parse_args()
    generateCanvas(options.inputFile,options.outputFile)
    
