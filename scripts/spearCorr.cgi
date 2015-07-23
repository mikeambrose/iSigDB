#!/UCSC/Pathways-Auxiliary/UCLApathways-PyPy-20150702/pypy-2.6-linux_x86_64-portable/bin/pypy
import os
import sys
sys.path.append('/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank')
import corrMatrix
import cgi
cgi.maxlen = 200 * 1024**2 #100mb
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
work_dir = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_{0}'.format(job_id)
os.makedirs(work_dir)
output_file = work_dir+'/{0}.txt'.format(job_id)
util.copyFile(userFile,output_file)

#the heatmap metric is called heatmap_metric
version = form["heatmap_metric"].value
#the row and column metrics are called row_metric and col_metric
rowMetric = form["row_metric"].value
colMetric = form["col_metric"].value
#the invert metric is called invert
invertMetric = "invert" in form
#set min and max if the scale checkbox is checked
if "scale" in form:
    mn = form["minVal"].value
    mx = form["maxVal"].value
else:
    mn,mx = None,None
#the gene selection metric is called spear_gene and has values spearGeneAll, spearGeneTop, spearGeneMag
geneMetric = form["spear_gene"].value
#depending on the geneMetric, we look at different forms for the gene value
geneVal = {'all':None,'top':int(form["matrix_num_genes"].value,'mag':int(form["matrix_mag"].value,\
            'cov':int(form["matrix_cov"].value)}[geneMetric]

#the signature matrix is in matrix
matrix_selected = form["matrix"].value
matrix_abbrevs = util.loadAbbrevs('/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/matrixAssociations.txt')
matrix_file = matrix_abbrevs[matrix_selected]

#logging

userIP = cgi.escape(os.environ["REMOTE_ADDR"])
logFileDir = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/iSigDB_uploads/useLog.txt' #TODO: add directory
with open(logFileDir,'a') as logFile:
    logFile.write('\t'.join([str(x) for x in [userIP,version,invertMetric,geneMetric,geneVal,rowMetric,colMetric]]))

print "Content-type: text/html\n\n"
corrMatrix.runCorrelation(output_file,version,invertMetric,mn,mx,rowMetric,colMetric,geneMetric,geneVal,matrix_file,False,job_id)
