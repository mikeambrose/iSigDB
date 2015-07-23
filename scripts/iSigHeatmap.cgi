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
if not 1 <= numIter <= 100000:
    util.displayErrorMessage("The number of iterations of the null distribution must be between 1 and 100000")
invert = "invert" in form
rowMetric = form["row_metric"].value
colMetric = form["col_metric"].value
rowColAcceptable = ["euclidean","pear_cor","none"]
for met in rowMetric,colMetric:
    if met not in rowColAcceptable:
        util.displayErrorMessage("Not a valid clustering metric: {0}".format(met),True)
n = form["num_genes"].value
acceptableN = [str(x) for x in [10,25,50,100,250,500,1000]]
if n not in acceptableN:
    util.displayErrorMessage("Not a valid number of genes: {0}".format(n),True)
n = int(n)
species = form["species"].value
acceptableSpecies = ["human","mouse"]
if species not in acceptableSpecies:
    util.displayErrorMessage("Not a valid species: {0}".format(species),True)
mn,mx = None,None
if "scale" in form:
    try:
        mn = float(form["mn"].value)
        mx = float(form["mx"].value)
    except Exception:
        util.displayErrorMessage("Min and max color range must be numbers")

#logging this call

userIP = cgi.escape(os.environ["REMOTE_ADDR"])
logFileDir = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/iSigDB_uploads/useLog.txt' #TODO: add directory
with open(logFileDir,'a') as logFile:
    logFile.write('\t'.join([str(x) for x in [userIP,version,zTransform,compNull,numIter,invert,rowMetric,colMetric,n,mn,mx,species]]))

print """Content-type: text/html

"""

sigComp.generateHeatmap(outputFile,'/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/SigGenes.txt',"{0}/abbrevs.txt".format(workDir),n,version,zTransform,jobID,rowMetric,colMetric,invert,compNull,False,numIter,mn,mx,species)
