#!/usr/bin/env python
import cgi
basePath = '/home/mike/workspace/PellegriniResearch/frontend/baseSiteV4.html'
assocPath = '/home/mike/workspace/PellegriniResearch/UserUploadedFiles/associations.txt'
replaceSection = '{FILES GO HERE}'
fileHTML = open(basePath).read()
loc = fileHTML.index(replaceSection)
fileHTML = fileHTML[:loc] + fileHTML[loc+len(replaceSection):]
selects = ""
with open(assocPath) as assocs:
    assocs.readline()
    for line in assocs:
        fname = line.split('\t')[0]
        selects = selects + "<option value=\"{0}\">{0}</option>\n".format(fname)
fileHTML = fileHTML[:loc] + selects + fileHTML[loc:]
#for debug only
#with open('output.html','w') as out:
#    out.write(fileHTML)
print "Content-type:text/html\n\n"+fileHTML
