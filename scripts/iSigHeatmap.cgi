#!/usr/bin/env python
import os
import cgi
cgi.maxlen = 100 * 1024**2
#TODO: add to python path
import iSigDBUtilities as util
#TODO: fix filename - don't think python allows .
import Sig_Avg_Matrix_Derm.RankV4 as sigComp
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

form = cgi.FieldStorage()

jobID = util.getJobID()
workDir = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_corrMatrix_{0}'.format(jobID)
os.makedirs(workDir)
outputFile = workDir + '/{0}.txt'.format(jobID)
with open(userFile) as userInput:
    util.copyFile(userInput,outputFile)

#write signatures to local abbrevs file
allAbbrevs = open('/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/abbrevs_fixed.txt').read()
fullToAbbrev = {}
for line in allAbbrevs.split('\n'):
    abbrev,full=line.split('\t')
    fullToAbbrev[full]=abbrev

if form["checkedSigs"].count(',') == 0:
    util.displayErrorMessage("No signatures selected",True)

with open("{0}/abbrevs.txt".format(workDir),'w') as localAbbrevs:
    selectedFulls = form["checkedSigs"].value.split(',')
    for full in selectedFulls:
        if full in fullToAbbrev: #filters out headers like 'mouse', 'human'
            localAbbrevs.write("{0}\t{1}\n".format(fullToAbbrev[full],full))

version = form["heatmap_metric"].value
zTransform = "matrix" if "scale_columns" in form else "none"
fixed = "fixed" in form
compNull = "null" in form
invert = "invert" in form
rowMetric = form["row_metric"].value
colMetric = form["col_metric"].value
n = int(form["num_genes"].value)
sigComp.generateHeatmap(outputFile,'/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/SigGenes.txt',"{0}/abbrevs.txt".format(workDir),n,version,zTransform,jobID,rowMetric,colMetric,invert,fixed,compNull,False)
