#!/usr/bin/env python
import os
import datetime
import corrMatrix
import cgi
cgi.maxlen = 100 * 1024**2 #100mb
#TODO: add this to the python path
import iSigDBUtilities as util
form = cgi.FieldStorage()
#gets relevant forms from the spearman html
#the default file is named matrix_file
#the server file is named serverFile
if 'serverFile' in form and form['serverFile'] != '':
    #the user has selected a file from the server
    #TODO
    pass
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
version = form["heatmap_metric"]
#the row and column metrics are called row_metric and col_metric
rowMetric = form["row_metric"]
colMetric = form["col_metric"]
#the gene selection metric is called spear_gene and has values spearGeneAll, spearGeneTop, spearGeneMag
geneMetric = form["spear_gene"]
#spearGeneAll -> nothing
if geneMetric == 'spearGeneAll':
    geneVal = None
#spearGeneTop -> matrix_num_genes
if geneMetric == 'spearGeneTop':
    geneVal = int(form["matrix_num_genes"])
#spearGeneMag -> matrix_mag
if geneMetric == 'spearGeneMag':
    geneVal = int(form["matrix_mag"])

#the signature matrix is in matrix
matrix_selected = form["matrix"]
#TODO: we should probably have some sort of abbrevs for matrices (or keep the naming consistent)
#TODO: fix abbrevs, give invert option on HTML
corrMatrix.runCorrelation(output_file,version,False,rowMetric,colMetric,geneMetric,geneVal,matrix_selected,False,job_id,None)
