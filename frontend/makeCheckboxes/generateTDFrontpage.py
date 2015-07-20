#!/usr/bin/env python
"""script written by Mike Ambrose, mikeambrose@berkeley.edu"""
from optparse import OptionParser

def generateJSTree(formatLines,delim="    "):
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
    returnHTML = returnHTML[5:-11]
    return returnHTML

def generateHTML(formatLines, baseSite, output):
    """Generates the HTML frontend (currently called tissueDeconvolutionV2/3
    formatLines is the format file
    output is the file which the html is written to
    returns nothing"""
    #returnHTML = open(baseSite).read().format(checkboxes=generateJSTree(formatLines),files="<option name=\"disabled\">")
    returnHTML = open(baseSite).read().replace("{checkboxes}",generateJSTree(formatLines))
    with open(output,'w') as f_write:
        f_write.write(returnHTML)
    print("successfully generated")

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-f","--format",dest="format",help="file with formatting")
    parser.add_option("-b","--base",dest="base",help="base html site pre-checkboxes")
    parser.add_option("-o","--output",dest="output",help="output filename")
    options,_ = parser.parse_args()
    generateHTML(open(options.format).read().split('\n'),options.base,options.output)
