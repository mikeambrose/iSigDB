"""A script to generate the html file which will display the heatmap"""
from optparse import OptionParser
import os.path

def replaceSections(html,replacements):
    """html is the raw text of the base html file
    replacements is a dictionary of html key to replacement text
    a replacement key a would look like {a} in the html text"""
    for replacement in replacements:
        assert "{"+replacement+"}" in html
        html = html.replace("{"+replacement+"}",replacements[replacement])
    return html

def generateCanvas(C,outFile,legendLabel='',invertHeatmap=True,centerAroundZero=False,givenMinVal=None,givenMaxVal=None,includeDetailed=True,showNull=False,showInput=False,showRDownload=False,tooltips=None,color=None,optionStr=''):
    """Generates the html for the heatmap page (for both signature- and matrix- based but not detailed)
    C is the set of constants
    outFile is a debug argument - it should always be called with None on the server
        if given a filename, instead of printing the html, it will write it to that file
    legendLabel is the label for the color legend
    invertHeatmap dictates whether the signatures are on the y- or x-axis
    centerAroundZero is true when delta, normalize, and a few other options are selected
        it makes min and max on either side of the axis
    givenMin/MaxVal set the color axis values
    includeDetailed dictates whether or not the detailed signature view is included
        it is false whenever called from the matrix tool, since there is no detail there
    show{Null,Input,RDownload} dictate whether or not their respective elements are linked to
    tooltips is a dictionary of xLabel:yLabel:p-value
        if none, not attached
    color is a color string which dictates what colors the color axis holds
        options are rwb (red->white->blue), wr (white->red), rw (red->white)
    optionStr is a human-readable description of the options selected to make the heatmap
    """
    replacements = {}
    f = open(C.R_TXT).read().split('\n')
    f = [x.split(',') for x in f]
    while not f[-1] or not any(f[-1]):    del f[-1]
    # the first row is the x labels
    xLabels = f[0][1:]
    # the first column is the y labels
    yLabels = [f[x][0] for x in range(1,len(f))]
    #now we generate the actual html
    #add the right location for scripts (use outFile != None to check if we are on server)
    if outFile != None: #we are in client
        base_loc = '.'
    else:
        base_loc = ''
    replacements["baseloc"] = base_loc
    if optionStr != '':
        replacements["options"] = "Selected options: {0}<br>".format(optionStr)
    else:
        replacements["options"] = ''
    #add links and hidden file location
    baseFile = 'http://pathways-pellegrini.mcdb.ucla.edu/submit/img/' + os.path.basename(C.R_HEATMAP)
    seed = baseFile[-30:-12]
    replacements['rPdf'] = baseFile
    if showNull:
        nullFilename="{0}img/{1}".format(C.ACCESIBLE_LINK,os.path.basename(C.NULL_PDF))
        replacements['null'] = "<a target=\"_blank\" href=" + nullFilename + ">Null model output</a><br>\n"
    else:
        replacements['null'] = ''
    if showInput:
        inpHistFilename="{0}img/{1}".format(C.ACCESIBLE_LINK,os.path.basename(C.INPUT_DIST_PDF))
        replacements["inpHist"] = "<a target=\"_blank\" href=" + inpHistFilename + ">Input distribution</a><br>\n"
    else:
        replacements["inpHist"] = ''
    if showRDownload:
        rDownloadableFilename = "{0}/data/{1}".format(C.ACCESIBLE_LINK,os.path.basename(C.R_DOWNLOAD))
        replacements["rDownload"] = "\n<a target=\"_blank\" href=" + rDownloadableFilename + ">Raw values</a><br>\n"
    else:
        replacements["rDownload"] = ''
    if includeDetailed:
        detailedText = '''
Look in more detail at one signature:
<form id="detailedHeatmap" name="detailedHeatmap" method="post" action="/cgi-bin/goTeles/sig_heatmapV4.cgi" ENCTYPE="multipart/form-data">
<input type="hidden" id="seed" name="seed" value="''' + seed + '''">
<select name="signatureName" id="signatureName">'''
        #add values for option
        for label in sorted(xLabels, key=lambda s: s.lower()):
            detailedText += "<option value=" + label + ">" + label.replace('"','') + "</option>\n"
        detailedText += """</select> <br/>
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
</form>"""
        replacements["includeDetailed"] = detailedText
    else:
        replacements["includeDetailed"] = ''
    #next add values
    minVal = float('+inf')
    maxVal = float('-inf')
    valText = ""
    for j in range(1,len(f[-1])):
        for i in range(1,len(f)):
            minVal = min(float(f[i][j]),minVal)
            maxVal = max(float(f[i][j]),maxVal)
            if not invertHeatmap:
                valText += str(j-1) + "," + str(i-1) + "," + f[i][j] + "\n" 
            else:
                valText += str(i-1) + "," + str(j-1) + "," + f[i][j] + "\n"
    replacements['vals'] = valText
    #and add tooltips
    if tooltips:
        tooltipText=""
        for j in range(len(yLabels)):
            for i in range(len(xLabels)):
                if not invertHeatmap:
                    tooltipText += str(j)+","+str(i)+","+str(tooltips[yLabels[j]][xLabels[i]])+"\n"
                else:
                    tooltipText += str(i)+","+str(j)+","+str(tooltips[yLabels[j]][xLabels[i]])+"\n"
        replacements['tooltips'] = tooltipText
    else:
        replacements['tooltips'] = ''
    if centerAroundZero:
        maxVal = max(abs(maxVal),abs(minVal))
        minVal = -maxVal
    if (centerAroundZero and color==None) or color=='bwr':
        c0,c25,c50,c75,c100 = '#00005C,#3060cf,#ffffff,#c4463a,#800000'.split(',')
    elif color==None or color=='wr':
        c0,c25,c50,c75,c100 = '#ffffff,#e7b5b0,#c4463a,#62231d,#800000'.split(',')
    elif color=='rw':
        c0,c25,c50,c75,c100 = '#800000,#62231d,#c4463a,#e7b5b0,#ffffff'.split(',')
    minVal = givenMinVal if (givenMinVal != None) else minVal
    maxVal = givenMaxVal if (givenMaxVal != None) else maxVal
    replacements['min'] = str(minVal)
    replacements['max'] = str(maxVal)
    replacements['legendLabel'] = legendLabel
    replacements['title'] = optionStr

    xLabelText = ""
    for label in xLabels if not invertHeatmap else yLabels:
        xLabelText += label+","
    replacements['xlabels'] = xLabelText

    yLabelText = ""
    for label in yLabels if not invertHeatmap else xLabels:
        yLabelText += label+","
    replacements['ylabels'] = yLabelText
    replacements['c0'] = c0
    replacements['c25'] = c25
    replacements['c50'] = c50
    replacements['c75'] = c75
    replacements['c100'] = c100
    htmlText = open(C.HEATMAP_TEMPLATE).read()
    htmlText = replaceSections(htmlText,replacements)
    if outFile:
        with open(outFile,'w') as f_out:
            f_out.write(htmlText)
    else:
        print htmlText
    
