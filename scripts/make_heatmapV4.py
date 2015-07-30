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

def generateCanvas(dataFile,outFile,legendLabel='',invertHeatmap=True,centerAroundZero=False,givenMinVal=None,givenMaxVal=None,baseFile='',includeDetailed=True,nullFilename=None,inpHistFilename=None,rDownloadFilename=None,tooltips=None,color=None,optionStr=''):
    """Generates the html for the heatmap page (for both signature- and matrix- based but not detailed)
    dataFile is the r-style output
    outFile is a debug argument - it should always be called with None on the server
        if given a filename, instead of printing the html, it will write it to that file
    legendLabel is the label for the color legend
    invertHeatmap dictates whether the signatures are on the y- or x-axis
    centerAroundZero is true when delta, normalize, and a few other options are selected
        it makes min and max on either side of the axis
    givenMin/MaxVal set the color axis values
    baseFile is the R heatmap output and also how we get the seed
    includeDetailed dictates whether or not the detailed signature view is included
        it is false whenever called from the matrix tool, since there is no detail there
    nullFilename is the name of the file containing the null distribution
        if it is '' or None, that file is not attached
    inpHistFilename is the name of the file containing the input distribution
        same rules as nullFilename
    rDownloadFilename has the data download from R
        same rules as nullFilename
    tooltips is a dictionary of xLabel:yLabel:p-value
        if none, not attached
    color is a color string which dictates what colors the color axis holds
        options are rwb (red->white->blue), wr (white->red), rw (red->white)
    optionStr is a human-readable description of the options selected to make the heatmap
    """
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
        <script src="""+'"'+base_loc+"""/heatmap_canvasV4.js">
    </script>
<script src=""" + '"'+base_loc+"""/HighCharts/js/highcharts.js"></script>
<script src="""+ '"'+base_loc+"""/HighCharts/js/modules/data.js"></script>
<script src="""+'"'+base_loc+"""/HighCharts/js/modules/heatmap.js"></script>
<script src="""+'"'+base_loc+"""/HighCharts/js/modules/exporting.js"></script>
    </head>
    <body>
    <h1>Results</h1>
<div id="container"></div> <br />
If you can't see the heatmap, make sure adblock is disabled and try again.<br>"""
    if optionStr != '':
        htmlText += "Selected options: {0}<br>".format(optionStr)
    htmlText += """
<b>Downloads: </b> <br />
<a target="_blank" href="""
    #add links and hidden file location
    baseFile = 'http://pathways-pellegrini.mcdb.ucla.edu/submit/img/' + os.path.basename(baseFile)
    seed = baseFile[-30:-12]
    htmlText += '"' + baseFile+'''">R heatmap output</a> <br />'''
    if nullFilename:
        htmlText += "\n<a target=\"_blank\" href=" + nullFilename + ">Null model output</a><br>\n"
    if inpHistFilename:
        htmlText += "\n<a target=\"_blank\" href=" + inpHistFilename + ">Input distribution</a><br>\n"
    if rDownloadFilename:
        htmlText += "\n<a target=\"_blank\" href=" + rDownloadFilename + ">Raw values</a><br>\n"
    if includeDetailed:
        htmlText += '''
Look in more detail at one signature:
<form id="detailedHeatmap" name="detailedHeatmap" method="post" action="/cgi-bin/goTeles/sig_heatmapV4.cgi" ENCTYPE="multipart/form-data">
<input type="hidden" id="seed" name="seed" value="''' + seed + '''">
<select name="signatureName" id="signatureName">'''
        #add values for option
        for label in sorted(xLabels, key=lambda s: s.lower()):
            htmlText += "<option value=" + label + ">" + label.replace('"','') + "</option>\n"
        htmlText += """</select> <br/>
<input type="checkbox" name="invert" id="invert" value="invert"> Invert output <br/>

Metric for sample clustering: &nbsp;
<select id="col_metric" name="col_metric">
<option value="euclidean">Euclidean Distance</option>
<option value="pear_cor">Pearson Correlation</option>
<option value="none">None (Do Not Cluster)</option>
</select><br />

Metric for gene clustering: &nbsp;
<select id="row_metric" name="row_metric">
<option value="euclidean">Euclidean Distance</option>
<option value="pear_cor">Pearson Correlation</option>
<option value="none">None (Do Not Cluster)</option>
</select><br />
<input type="submit" value="Go">
</form>
</body>
"""
    htmlText += """
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
    #and add tooltips
    htmlText += """</pre>
<pre id="tooltips" style="display:none">"""
    if tooltips:
        for j in range(len(yLabels)):
            for i in range(len(xLabels)):
                if not invertHeatmap:
                    htmlText += str(j)+","+str(i)+","+str(tooltips[yLabels[j]][xLabels[i]])+"\n"
                else:
                    htmlText += str(i)+","+str(j)+","+str(tooltips[yLabels[j]][xLabels[i]])+"\n"
    if (centerAroundZero and color==None) or color=='bwr':
        maxVal = max(abs(maxVal),abs(minVal))
        minVal = -maxVal
        c0,c25,c50,c75,c100 = '#00005C,#3060cf,#ffffff,#c4463a,#800000'.split(',')
    elif color==None or color=='wr':
        c0,c25,c50,c75,c100 = '#ffffff,#e7b5b0,#c4463a,#62231d,#800000'.split(',')
    elif color=='rw':
        c0,c25,c50,c75,c100 = '#800000,#62231d,#c4463a,#e7b5b0,#ffffff'.split(',')
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
        htmlText += label+","
    htmlText = htmlText[:-1] + """</pre>
<pre id="ylabels" style="display: none">
"""
    for label in yLabels if not invertHeatmap else xLabels:
        htmlText += label+","
    htmlText = htmlText[:-1] + """</pre>
<pre id="color0" style="display:none">{0}</pre>
<pre id="color25" style="display:none">{1}</pre>
<pre id="color50" style="display:none">{2}</pre>
<pre id="color75" style="display:none">{3}</pre>
<pre id="color100" style="display:none">{4}</pre>
</html>""".format(c0,c25,c50,c75,c100)
    if outFile:
        with open(outFile,'w') as f_out:
            f_out.write(htmlText)
    else:
        print htmlText
    
