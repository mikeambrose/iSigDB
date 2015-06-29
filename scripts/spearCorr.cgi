#!/usr/bin/env python
import os
import corrMatrix
import cgi
cgi.maxlen = 100 * 1024**2 #100mb
#TODO: add this to the python path
import iSigDBUtilities as util
form = cgi.FieldStorage()
#the default file is named matrix_file
#the server file is named serverFile
if 'serverFile' in form and form['serverFile'].value != '':
    #the user has selected a file from the server
    with open('/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/iSigDB_uploads/associations.txt')\
        as assoc:
        assoc.readline() #skip header line
        for line in assoc:
            if line.split('\t')[0] == form['serverFile'].value:
                userFileLoc = line.replace("\n","").split('\t')[1]
                userFile = open('/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/iSigDB_uploads/'+\
                            userFileLoc)
else:
    #get the file from the regular upload
    userFile = form["matrix_file"].file
    if not userFile:
        util.displayErrorMessage("No file uploaded",True)

#make directory
job_id = util.getJobID()
work_dir = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_corrMatrix_{0}'.format(job_id)
os.makedirs(work_dir)
output_file = work_dir+'/{0}.txt'.format(job_id)
with open(userFile) as userInput:
    util.copyFile(userInput,output_file)

#the heatmap metric is called heatmap_metric
version = form["heatmap_metric"].value
#the row and column metrics are called row_metric and col_metric
rowMetric = form["row_metric"].value
colMetric = form["col_metric"].value
#the invert metric is called invert
invertMetric = "invert" in form
#the gene selection metric is called spear_gene and has values spearGeneAll, spearGeneTop, spearGeneMag
geneMetric = form["spear_gene"].value
#spearGeneAll -> nothing
if geneMetric == 'spearGeneAll':
    geneVal = None
#spearGeneTop -> matrix_num_genes
if geneMetric == 'spearGeneTop':
    geneVal = int(form["matrix_num_genes"].value)
#spearGeneMag -> matrix_mag
if geneMetric == 'spearGeneMag':
    geneVal = int(form["matrix_mag"].value)

#the signature matrix is in matrix
matrix_selected = form["matrix"].value

#logging

userIP = cgi.escape(os.environ["REMOTE_ADDR"])
logFileDir = '' #TODO: add directory
with open(logFileDir,'a') as logFile:
    logFile.write('\t'.join([str(x) for x in [userIP,version,interMetric,geneMetric,geneVal,rowMetric,colMetric]]))

#TODO: get path to matrix_selected working somehow
print "Content-type: text/html\n\n"
corrMatrix.runCorrelation(output_file,version,invertMetric,rowMetric,colMetric,geneMetric,geneVal,matrix_selected,False,job_id,None)
