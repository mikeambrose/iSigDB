
from optparse import OptionParser
import numpy as np
import glob
import os
import gzip
import subprocess
import scipy.stats
S_GENE_COUNT = '50'
S_VERSION = 'rank_delta'
S_ZSCORE = 'column'

def replaceLineEndings(f):
    s = open(f).read()
    unix = '\n' in s
    dos = '\r' in s
    if unix and dos:
        s = s.replace('\r\n','\n')
        if '\r' in s: #we didn't successfully rid ourselves
            errorMessage("File has non-standard line endings")
    elif dos:
        s = s.replace('\r','\n')
    if dos:
        with open(f,'w') as fout:
            fout.write(s)

def errorMessage(s):
    print "<!doctype html> <html> <body>" + s + "</body> </html>"
    exit()

def checkForErrors(f):
    """Checks that f is vaguely in the format we expect
    Checks that it appears to be tab-delimited with a consistent number of columns
    and that all entries are decimals.
    If not, prints out html corresponding to an error message and exits execution"""
    f = open(f).read().split('\n')
    if len(f) <= 1:
        errorMessage("Formatting error: Only one line detected")
    numTabs = f[0].count('\t')
    if numTabs == 0:
        errorMessage("Not a tab-separated file")
    for i in range(1,len(f)-1):
        line = f[i]
        if line.count('\t') != numTabs:
            errorMessage("Inconsistent number of columns around line " + str(i+1))
        for x in line.split('\t')[1:]:
            try:
                float(x)
            except:
                errorMessage("Non-decimal alue around line " + str(i+1))

def getSigGenes(strSigFile,nGeneCount):
    dSigToGene = {}
    with open(strSigFile) as sigToGenes:
        for line in sigToGenes:
            splitLine = line.strip("\n").split("\t")
            sig = splitLine[0]
            dSigToGene[splitLine[0]] = splitLine[1:nGeneCount+1]
    #import pdb; pdb.set_trace()
    return dSigToGene
        
def getSamToGeneToRank(dSamToCountToGene):

    dSamToGeneToRank = {}
    dGeneToMeanRank = {}

    #process rank
    for strCurSam in dSamToCountToGene:
        nRank = 1
        dSamToGeneToRank[strCurSam] = {}
        for fCount in sorted(dSamToCountToGene[strCurSam].keys()):
            for strCurGene in sorted(dSamToCountToGene[strCurSam][fCount]):
                dSamToGeneToRank[strCurSam][strCurGene] = nRank
                nRank+=1

    #process mean gene rank
    for strCurGene in sorted(dSamToGeneToRank[dSamToGeneToRank.keys()[0]].keys()):
        lRowRanks = []
        for strCurSam in sorted(dSamToGeneToRank.keys()):
            lRowRanks.append(dSamToGeneToRank[strCurSam][strCurGene])
            dGeneToMeanRank[strCurGene] = np.average(lRowRanks)
    #QC print rankings

    #fout = open('./test_rank.log','w')
    #strHeader = 'GENE'
    #for strCurSam in sorted(dSamToGeneToRank.keys()):
    #    strHeader+='\t'+strCurSam
    #fout.write(strHeader+'\tAVG\n')
    #for strCurGene in sorted(dSamToGeneToRank[dSamToGeneToRank.keys()[0]].keys()):
    #    strRow = strCurGene
    #    for strCurSam in sorted(dSamToGeneToRank.keys()):
    #        strRow+='\t%d'%dSamToGeneToRank[strCurSam][strCurGene]
    #    fout.write(strRow+'\t%f\n'%dGeneToMeanRank[strCurGene])
    #fout.close()

    return dSamToGeneToRank,dGeneToMeanRank

def procGeneCountMatrix(strGeneCount,dSigToGenes,lSigs,strOutFile,strVersion, bIsLog,sigNames):
    dColIDToColLbl = {}
    dSamToGeneToCount = {}
    dSamToSigToLVals = {}
    dSamToCountToGene = {}
    lGlobalAvg = []

    #load Sample To Gene To Count
    fopen = None
    if strGeneCount.endswith('.gz'):
        fopen = gzip.open(strGeneCount,'rb')
    else:
        fopen = open(strGeneCount,'r')

    for strLine in fopen:
        lcols = strLine.rstrip().split('\t')
        if len(dColIDToColLbl.keys())==0:
            for i in range(1,len(lcols)):
                dColIDToColLbl[i] = lcols[i]
                dSamToGeneToCount[lcols[i]] = {}
                dSamToSigToLVals[lcols[i]] = {}
                dSamToCountToGene[lcols[i]] = {}
                for strSig in dSigToGenes.keys():
                    dSamToSigToLVals[lcols[i]][strSig] = []
        else:
            strCurGene = lcols[0]
            strCurGene = strCurGene.upper()
            for i in range(1,len(lcols)):
                dSamToGeneToCount[dColIDToColLbl[i]][strCurGene] = float(lcols[i])
                if dSamToCountToGene[dColIDToColLbl[i]].has_key(float(lcols[i])) == False:
                    dSamToCountToGene[dColIDToColLbl[i]][float(lcols[i])] = []
                dSamToCountToGene[dColIDToColLbl[i]][float(lcols[i])].append(strCurGene)
    fopen.close()

    #qc print gene
    #===================================================================================================================
    #write avg log gene count
    #===================================================================================================================
    if strVersion == 'log':
        for strSig in dSigToGenes.keys():
            for strCurSigGene in dSigToGenes[strSig]:
                for strSam in dSamToGeneToCount.keys():
                    if(dSamToGeneToCount[strSam].has_key(strCurSigGene) == False):
                        dSamToGeneToCount[strSam][strCurSigGene] = 1
                    #dSamToSigToLVals[strSam][strSig].append(np.log10(dSamToGeneToCount[strSam][strCurSigGene]+1))
                    if bIsLog:
                        dSamToSigToLVals[strSam][strSig].append(np.log10(dSamToGeneToCount[strSam][strCurSigGene]+1))
                    else:
                        dSamToSigToLVals[strSam][strSig].append(dSamToGeneToCount[strSam][strCurSigGene])
                    #print '%s\t%s\t%f'%(strSam,strCurSigGene,dSamToGeneToCount[strSam][strCurSigGene])

        # TODO: extra output here
        fout2 = open(strOutFile+'.full.txt', 'w')
        for strSig in sorted(dSigToGenes.keys()):
            # import pdb; pdb.set_trace()
            fout2.write('# Signature: %s\n' % sigNames[strSig] if strSig in sigNames else strSig)
            for strSam in sorted(dSamToGeneToCount.keys()):
                fout2.write('\t%s' % strSam)
            fout2.write('\n')
            for strCurSigGene in dSigToGenes[strSig]:
                fout2.write('%s' % strCurSigGene)
                for strSam in sorted(dSamToGeneToCount.keys()):
                    fout2.write('\t%f' % dSamToGeneToCount[strSam][strCurSigGene])
                fout2.write('\n')
            fout2.write('\n')

        #print header
        fout = open(strOutFile,'w')
        strHeader = 'SAMPLE'
        for i in range(len(lSigs)):
            strHeader+='\t'+(sigNames[lSigs[i]] if lSigs[i] in sigNames else lSigs[i])
        fout.write(strHeader+'\n')
        #print sample Row Averages
        for strSam in sorted(dSamToGeneToCount.keys()):
            strRow = strSam
            for j in range(len(lSigs)):
                strCurSig = lSigs[j]
                strRow += '\t%f'%np.average(dSamToSigToLVals[strSam][strCurSig])
            fout.write(strRow+'\n')
        fout.close()

    #===================================================================================================================
    #write delta rank metric
    #===================================================================================================================
    else:
        dSamToGeneToRank,dGeneToMeanRank = getSamToGeneToRank(dSamToCountToGene)

        #dSamToSigToLVals DELTA MEAN MATRIX
        for strSig in dSigToGenes.keys():
            for strCurSigGene in dSigToGenes[strSig]:
                for strSam in dSamToGeneToRank.keys():
                    if(dSamToGeneToRank[strSam].has_key(strCurSigGene) == False):
                        continue
                    if strVersion.startswith('rank_avg'):
                        dSamToSigToLVals[strSam][strSig].append(dSamToGeneToRank[strSam][strCurSigGene])
                        lGlobalAvg.append(dSamToGeneToRank[strSam][strCurSigGene])
                    elif strVersion.startswith('rank_delta'):
                        dSamToSigToLVals[strSam][strSig].append(dSamToGeneToRank[strSam][strCurSigGene] - dGeneToMeanRank[strCurSigGene])
                        lGlobalAvg.append(dSamToGeneToRank[strSam][strCurSigGene] - dGeneToMeanRank[strCurSigGene])
        #lLoadGlobalAvg
        '''
        for strSam in sorted(dSamToGeneToCount.keys()):
            strRow = strSam
            for j in range(len(lSigs)):
                strCurSig = lSigs[j]
                lGlobalAvg.append(np.average(dSamToSigToLVals[strSam][strCurSig]))
        '''

        # TODO: extra output here
        fout2 = open(strOutFile+'.full.txt', 'w')
        for strSig in sorted(dSigToGenes.keys()):
            fout2.write('# Signature: %s\n' % sigNames[strSig] if strSig in sigNames else strSig)
            for strSam in sorted(dSamToGeneToRank.keys()):
                fout2.write('\t%s' % strSam)
            fout2.write('\n')
            for strCurSigGene in dSigToGenes[strSig]:
                fout2.write('%s' % strCurSigGene)
                for strSam in sorted(dSamToGeneToRank.keys()):
                    if(dSamToGeneToRank[strSam].has_key(strCurSigGene) == False):
                        fout2.write('\tN/A')
                    else:
                        fout2.write('\t%f' % dSamToGeneToRank[strSam][strCurSigGene])
                fout2.write('\n')
            fout2.write('\n')

        #find signature gene
        fout = open(strOutFile,'w')
        strHeader = 'SAMPLE'
        for i in range(len(lSigs)):
            strHeader+='\t'+(sigNames[lSigs[i]] if lSigs[i] in sigNames else lSigs[i])
        fout.write(strHeader+'\n')
        #print sample Row Averages
        for strSam in sorted(dSamToGeneToCount.keys()):
            strRow = strSam
            for j in range(len(lSigs)):
                strCurSig = lSigs[j]
                if strVersion.endswith('_z'):
                    fMetric = (np.average(dSamToSigToLVals[strSam][strCurSig]) - np.average(lGlobalAvg))/np.std(lGlobalAvg)
                else:
                    fMetric = np.average(dSamToSigToLVals[strSam][strCurSig])
                strRow += '\t%f'%fMetric
            fout.write(strRow+'\n')
        fout.close()

#===============================
#functions for Spearman
#==============================

def spearmanRank(sigs,sams,geneNames,strOutFile,sigNames):
    """Uses the spearman rank coefficient to evaluate the similarity between signatures and samples
    sigs is a dictionary where each entry has the n genes chosen for analysis
    sams is in the same format as sigDict, but has samples as keys instead of signatures
    geneNames is a list of gene names in the same order as they appear in sigs and sams
    sigNames is a dictionary of sig abbreviation : full name
    Similar to other functions, writes to two files:
        strOutFile contains the basic information used by the R script
        strOutFile + ".matrix.txt" has the detailed information used in generating the detailed heatmaps
    No value is returned"""
    #stored as a dictionary of (sig,sam):coef
    sigSamCoefficient = {}
    #dictionary of signature:gene:list of d values corresponding to (sig,gene,sam)
    geneDetail = {}
    for sig in sigs:
        geneDetail[sig] = {}
        for gene in geneNaes:
            geneDetail[sig][gene] = []
    for sig in sigs:
        for sam in sams:
            coef,rankDiffs = getSpearmanVals(sigs[sig],sams[sam])
            for i in range(len(geneNames)):
                geneDetail[sig][geneNames[i]].append(coef[i])
            sigSamCoefficient[(sig,sam)] = coef
    #write detailed output
    with open(strOutFile+".matrix.txt","w") as detailedOut:
        for sig in sigs:
            detailedOut.write('# Signature: ' + sigNames[sig] if sig in sigNames else sig+'\n')
            #write sample names
            for sam in sams:
                detailedOut.write('\t'+sam)
            #write gene detail
            for gene in geneNames:
                detailedOut.write(gene)
                samVals = geneDetail[sig][gene]
                assert len(samVals) == len(sams)
                for val in samVals:
                    detailedOut.write('\t'+val)
                detailedOut.write('\n')
            detailedOut.write('\n')
    #write regular output
    with open(strOutFile,"w") as basicOut:
        #write header
        basicOut.write("SAMPLE")
        for sig in sigs:
            #use the full name if it exists, otherwise use the abbreviation
            basicOut.write("\t"+sigNames[sig] if sig in sigNames else sig)
            basicOut.write("\n")
        for sam in sams:
            basicOut.write(sam)
            for sig in sigs:
                basicOut.write("\t"+sigSamCoefficient[(sig,sam)])
            basicOut.write("\n")

def getRanks(lst):
    """Returns a dictionary of entry : rank"""
    ranks = {}
    sorted_lst = sorted(lst)
    i = 0
    while i < len(lst):
        #detect if there is a run
        if i+1<len(lst) and lst[i+1] == lst[i]:
            run_length = 2
            j = i+1
            while j+1<len(lst) and lst[j+1] == lst[j]:
                run_length += 1
                j += 1
            rankVal = i+(run_length-1)/2.0
            for k in range(i,j+1):
                ranks[lst[k]] = rankVal
        else:
            ranks[lst[i]] = i
    return ranks

def getSpearmanVals(v1,v2):
    """Returns the spearman correlation coefficient between v1 and v2
    as well as the list of rank differentials between v1 and v2"""
    n = len(v1)
    v1Ranks,v2Ranks = getRanks(v1),getRanks(v2)
    d = [v1Ranks[v1[i]] - v2Ranks[v2[i]] for i in range(n)]
    rho = 1 - 6*sum(x**2 for x in d)/float(n*(n**2-1))
    return rho,d
    
def getSpearmanGenes():
    return ["DYNAP","CHRNA1","POPDC3"]

def getSpearmanDict(inputDict,genes):
    """Contructs a matrix of key : values for each gene"""
    returnDict = {}
    for line in inputDict:
        returnDict[line] = []
        for gene in genes:
            returnDict[line].append(inputDict[line][gene])
    return returnDict

def getSampleDict(f):
    """returns a dictionary of sample:gene:value"""
    samDict = {}
    f = open(f).read().split('\n')
    sams = f[0].split('\t')[1:]
    #for line in 

#===============================
#end functions for Spearman
#==============================

def loadGroupInfo(strPathToGroupFile):
    dGroupToLTisDesc = {}
    strGroup = strPathToGroupFile.split('.')[0]
    dGroupToLTisDesc[strGroup] = []
    for strLine in open(strPathToGroupFile,'r'):
        lcols = strLine.rstrip().split('\t')
        strTissue = lcols[0]
        strDesc = lcols[1]
        dGroupToLTisDesc[strGroup].append((strTissue,strDesc))
    return dGroupToLTisDesc

def loadSigDictionary(strPathToGroupFile):
    returnDict = {}
    f = open(strPathToGroupFile,'r')
    for line in f:
        abbrev, real = line.split('\t')
        returnDict[abbrev] = real[:-2]
    return returnDict

def createHeatMap(strMatrixFile,strOutFile,strVersion,strColumnZ,strRowMetric,strColMetric,invert,fixed):
    #print strMatrixFile
    #print strOutFile
    strTxtOutFile = 'routput.txt'
    strDend1OutFile = 'dend1.pdf'
    strDend2OutFile = 'dend2.pdf'
    subprocess.call(["Rscript",'heatsig.R',strMatrixFile,strOutFile,strColumnZ,strRowMetric,strColMetric,strTxtOutFile, strDend1OutFile, strDend2OutFile])
    import make_heatmap
    centerAroundZero = (strColumnZ == 'matrix') or (strVersion == "rank_delta")
    maxVal,minVal = None,None
    if fixed:
        if strColumnZ == 'matrix':
            minVal = -5
            maxVal = 5
        else:
            if strVersion == "rank_delta":
                maxVal = 10000
                minVal = -10000
            elif strVersion == "rank_avg":
                maxVal = 30000
                minVal = 0
            elif strVersion == "log":
                maxVal = 10
                minVal = 0
            else:
                maxVal = 30
                minVal = 0
    make_heatmap.generateCanvas(strTxtOutFile,'heatmapoutput.html','Matrix Z-Score' if strColumnZ == 'matrix' else 'Value',
                                    invert,centerAroundZero,minVal,maxVal)
    #print'Done'

#----------------------------------------------------------------------------
# main function call
#----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--gene_count", dest="gene_count", help="table of gene counts", metavar="TAB", default="")
    parser.add_option("-s", "--sig", dest="sig", help="signature gene list", metavar="LIST", default="")
    parser.add_option("-g", "--group", dest="group", help="path to group files", metavar="PATH", default="")
    parser.add_option("-n", "--ngene_count", dest="ngene_count", help="number of genes to use", metavar="NUM", default=S_GENE_COUNT)
    parser.add_option("-v", "--version", dest="version", help="metric version (rank_avg,rank_delta,log)", metavar="VER", default=S_VERSION)
    parser.add_option("-z", "--zscore", dest="zscore", help="default no column z-score", metavar="ZSCORE", default=S_ZSCORE)
    parser.add_option("-j", "--job_id", dest="job_id", help="job id, for output path", default="-1")
    parser.add_option("-a", "--absolute", dest="bIsLog", action="store_false", default=True, help="no log transformation")
    parser.add_option("-l", "--log_transform", dest="bIsLog", action="store_true", help="log transform")
    parser.add_option("-r", "--row_metric", dest="row_metric", help="metric for clustering rows (samples)", default="pear_cor")
    parser.add_option("-c", "--col_metric", dest="col_metric", help="metric for clustering columns (signatures)", default="pear_cor")
    parser.add_option("-i", "--invert",default=True,dest="invert",help="heatmap columns inverted")
    parser.add_option("-f", "--fixed", dest="fixed", default="none", help="use fixed color axis")
    (options, args) = parser.parse_args()
    #check that it's in the form we expect and replace non-unix line endings
    replaceLineEndings(options.gene_count)
    checkForErrors(options.gene_count)
    strOutMatrixTxt = 'matrix.txt'   
    if options.version != "spearman":
        #load categories
        dGroupToLTisDesc = loadGroupInfo(options.group)
        #list signatures
        lTis = []
        for strGroup in sorted(dGroupToLTisDesc.keys()):
            for strTis,strDes in dGroupToLTisDesc[strGroup]:
                lTis.append(strTis)
        #load list of genes
        dGroupSigToGene = getSigGenes(options.sig,int(options.ngene_count))

        #process gene count matrix
        procGeneCountMatrix(options.gene_count,dGroupSigToGene,lTis,strOutMatrixTxt,options.version,options.bIsLog,loadSigDictionary(options.group))
    else:
        genes = getSpearmanGenes()
        #load dictionaries for sigs
        sigs = getSpearmanSigs(options.gene_count,genes)
        #load dictionary for samples
        sams = getSpearmanSams(options.sig,genes)
        spearmanRank(sigs,sams,genes,strOutMatrixTxt,loadSigDictionary(options.group))
    #generate heatmap pdf
    createHeatMap(strOutMatrixTxt,'heatmap.pdf',options.version,options.zscore,options.row_metric,options.col_metric,False if options.invert=='none' else True,False if options.fixed=='none' else True)
