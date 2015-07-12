#!/usr/bin/env python
import os
import cgi
cgi.maxlen = 100 * 1024**2
#TODO: add to python path
import sys
sys.path.append('/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank')
import iSigDBUtilities as util
import SigAvg as sigComp

form = cgi.FieldStorage()

if 'serverFile' in form and form['serverFile'].value != '':
    #the user has selected a file from the server
    with open('/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/iSigDB_uploads/associations.txt')\
        as assoc:
        assoc.readline() #skip header line
        userFile = None
        for line in assoc:
            if line.split('\t')[0] == form['serverFile'].value:
                userFileLoc = line.replace("\n","").split('\t')[1]
                userFile = open('/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/iSigDB_uploads/'+\
                            userFileLoc)
                break
        if not userFile:
            util.displayErrorMessage("No such file found",True)
else:
    #get the file from the regular upload
    userFile = form["matrix_file"].file
    if not userFile:
        util.displayErrorMessage("No file uploaded",True)

jobID = util.getJobID()
workDir = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_{0}'.format(jobID)
os.makedirs(workDir)
outputFile = workDir + '/{0}.txt'.format(jobID)
util.copyFile(userFile,outputFile)

#write signatures to local abbrevs file
allAbbrevs = open('/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/abbrevs_fixed.txt').read()
fullToAbbrev = {}
for line in allAbbrevs.split('\n'):
    if not line:
      continue
    abbrev,full=line.split('\t')
    fullToAbbrev[full]=abbrev

if form["checkedSigs"].value.count(',') == 0:
    util.displayErrorMessage("No signatures selected",True)

with open("{0}/abbrevs.txt".format(workDir),'w') as localAbbrevs:
    selectedFulls = form["checkedSigs"].value.split(',')
    for full in selectedFulls:
        if full in fullToAbbrev: #filters out headers like 'mouse', 'human'
            localAbbrevs.write("{0}\t{1}\n".format(fullToAbbrev[full],full))
version = ''
for option in ["log","delta","rank"]:
    #add option string to version if it's in the form
    version = version + option * (option in form)
zTransform = "matrix" if "scale_columns" in form else "none"
compNull = "null" in form
numIter = int(form["nullNumIter"].value)
invert = "invert" in form
rowMetric = form["row_metric"].value
colMetric = form["col_metric"].value
n = int(form["num_genes"].value)
mn,mx = None,None
if "scale" in form:
    mn = form["mn"].value
    mx = form["mx"].value

#logging this call

userIP = cgi.escape(os.environ["REMOTE_ADDR"])
logFileDir = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/iSigDB_uploads/useLog.txt' #TODO: add directory
with open(logFileDir,'a') as logFile:
    logFile.write('\t'.join([str(x) for x in [userIP,version,zTransform,compNull,numIter,invert,rowMetric,colMetric,n,mn,mx]]))

print """Content-type: text/html

"""

sigComp.generateHeatmap(outputFile,'/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/SigGenes.txt',"{0}/abbrevs.txt".format(workDir),n,version,zTransform,jobID,rowMetric,colMetric,invert,compNull,False,numIter,mn,mx)
