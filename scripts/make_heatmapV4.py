"""A script to generate the html file which will display the heatmap"""
from optparse import OptionParser
import os.path
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
"""



def generateCanvas(dataFile,outFile,legendLabel='',invertHeatmap=True,centerAroundZero=False,givenMinVal=None,givenMaxVal=None,baseFile=''):
    f = open(dataFile).read().split('\n')
    f = [x.split(',') for x in f]
    while not f[-1] or not any(f[-1]):    del f[-1]
    # the first row is the x labels
    xLabels = f[0][1:]
    # the first column is the y labels
    yLabels = [f[x][0] for x in range(1,len(f))]
    #now we generate the actual html
    htmlText = HTML_BASE
    #add the right location for scripts (use outFile != None to check if we are on server)
    if outFile != None: #we are in client
        base_loc = '.'
    else:
        base_loc = ''
    htmlText += """
        <script src="""+'"'+base_loc+"""/heatmap_canvas.js">
    </script>
    </head>
    <body>
<script src=""" + '"'+base_loc+"""/HighCharts/js/highcharts.js"></script>
<script src="""+ '"'+base_loc+"""/HighCharts/js/modules/data.js"></script>
<script src="""+'"'+base_loc+"""/HighCharts/js/modules/heatmap.js"></script>
<script src="""+'"'+base_loc+"""/HighCharts/js/modules/exporting.js"></script>

<div id="container"></div> <br />
If you can't see the heatmap, make sure adblock is disabled and try again.<br>
<b>Downloads: </b> <br />
<a href="""
    #add links and hidden file location
    baseFile = 'http://pathways-pellegrini.mcdb.ucla.edu/submit/img/' + os.path.basename(baseFile)
    seed = baseFile[-18:]
    htmlText += '"' + baseFile+'''heatmap.pdf">R heatmap output</a> <br />
Look in more detail at one signature:
<form id="detailedHeatmap" name="detailedHeatmap" method="post" action="/cgi-bin/goTeles/sig_heatmap.cgi" ENCTYPE="multipart/form-data">
<input type="hidden" id="seed" name="seed" value="''' + seed + '''">
<select name="signatureName" id="signatureName">'''
    #add values for option
    for label in sorted(xLabels, key=lambda s: s.lower()):
        htmlText += "<option value=" + label + ">" + label.replace('"','') + "</option>\n"
    htmlText += """</select> <br/>
<input type="checkbox" name="invert" id="invert" value="invert"> Invert output <br/>

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
    if outFile:
        with open(outFile,'w') as f_out:
            f_out.write(htmlText)
    else:
        print htmlText
    
