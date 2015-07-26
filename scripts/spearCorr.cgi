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
versionRestrictions = ["pearson","spearman","deconv"]
if version not in versionRestrictions:
    util.displayErrorMessage("Invalid version selected - {0}".format(version),True)
#the row and column metrics are called row_metric and col_metric
rowMetric = form["row_metric"].value
colMetric = form["col_metric"].value
rowColRestrictions = ["euclidean","pear_cor","none"]
    if met not in rowColRestrictions:
        util.displayErrorMessage("Invalid metric selected - {0}".format(met),True)
#the invert metric is called invert
invertMetric = "invert" in form
#set min and max if the scale checkbox is checked
if "scale" in form:
    mn = form["minVal"].value
    mx = form["maxVal"].value
    for val in mn,mx:
        try:
            x = float(val)
            if not -1 <= x <= 1:
                util.displayErrorMessage("Invalid range for color axes - must range from -1 to 1",True)
        except Exception:
            util.displayErrorMessage("Color axis ranges must be a number",True)

else:
    mn,mx = None,None
#the gene selection metric is called spear_gene and has values spearGeneAll, spearGeneTop, spearGeneMag
geneMetric = form["spear_gene"].value
#depending on the geneMetric, we look at different forms for the gene value
geneValues = {'all':None,'top':int(form["matrix_num_genes"].value),'mag':int(form["matrix_mag"].value),\
            'cov':int(form["matrix_cov"].value)}
geneRestrictions = {'all':[None],'top':[10,25,50,100,250,1000],'mag':[2,5,10,50],'cov':[500,1000,2500,5000]}
geneVal = geneVales[geneMetric]
if geneVal not in geneRestrictions[geneMetric]:
    util.displayErrorMessage("Invalid options selected - number of genes cannot be {0} when in mode {1}".format(geneVal,geneRestrictions),True)

#the signature matrix is in matrix
matrix_selected = form["matrix"].value
matrix_abbrevs = util.loadAbbrevs('/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/matrixAssociations.txt')
matrix_file = matrix_abbrevs[matrix_selected]

#logging

userIP = cgi.escape(os.environ["REMOTE_ADDR"])
logFileDir = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/iSigDB_uploads/useLog.txt'
with open(logFileDir,'a') as logFile:
    logFile.write('\t'.join([str(x) for x in ['matrix',userIP,version,invertMetric,geneMetric,geneVal,rowMetric,colMetric]]))

print "Content-type: text/html\n\n"
corrMatrix.runCorrelation(output_file,version,invertMetric,mn,mx,rowMetric,colMetric,geneMetric,geneVal,matrix_file,False,job_id,False)
