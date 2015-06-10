#!/usr/bin/env python
"""Written by Mike Ambrose, mikeambrose@berkeley.edu
Takes the formatting file and writes html corresponding to a jsTree with the same structure
"""

def getJSTree(lines,output,delim="    "):
    lastInd = 0 #indendation of previous line (number of tabs)
    returnHTML = ''
    for line in lines:
        currentInd = line.count(delim)
        if currentInd == lastInd: #no change
            returnHTML += "</li><li>"+line.replace(delim,'') + '\n'
        elif currentInd > lastInd: #up one
            assert (currentInd-lastInd) == 1
            returnHTML += "<ul><li>"+line.replace(delim,'') + '\n'
        else:
            spacingDiff = lastInd-currentInd
            returnHTML += "</li>"+"</ul></li>"*spacingDiff + "<li>"+line.replace(delim,'') + '\n'
        lastInd = currentInd
    #get rid of first </li><li> and last <li> and adding auto-select
    returnHTML = "<li data-jstree='{\"selected\":true}'>"+returnHTML[9:-5]
    with open(output,'w') as out:
        out.write(returnHTML)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-i","--input",dest="input",help="formatting file")
    parser.add_option("-o","--output",dest="output",help="output file")
    (options,args) = parser.parse_args()
    getJSTree(open(options.input).read().split('\n'),options.output)
