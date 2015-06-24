#==============================================================================
# Date Modified: 2015.04.24
#==============================================================================
import os
from optparse import OptionParser
import numpy as np
import glob
import os
import gzip
import subprocess
import math
from collections import OrderedDict
import iSigDBUtilities as util
S_GENE_COUNT = '50'
S_VERSION = 'rank_delta'
S_ZSCORE = 'column'

def getSigGenes(strSigFile,nGeneCount):
    dSigToGene = {}
    with open(strSigFile) as sigToGenes:
        for line in sigToGenes:
            splitLine = line.strip("\n").split("\t")
            sig = splitLine[0]
            dSigToGene[splitLine[0]] = splitLine[1:nGeneCount+1]
    return dSigToGene

def getSamToGeneToRank(dSamToCountToGene):
    dSamToGeneToRank = OrderedDict()
    dGeneToMeanRank = {}

    #process rank
    for strCurSam in dSamToCountToGene:
        nRank = 1
        currentRankVal = None
        currentRunLength = 0
        dSamToGeneToRank[strCurSam] = {}
        for fCount in sorted(dSamToCountToGene[strCurSam].keys()):
            meanRank = nRank + len(dSamToCountToGene[strCurSam][fCount])/2.0
            for strCurGene in sorted(dSamToCountToGene[strCurSam][fCount]):
                dSamToGeneToRank[strCurSam][strCurGene] = meanRank
                nRank+=1

    #process mean gene rank
    for strCurGene in sorted(dSamToGeneToRank[dSamToGeneToRank.keys()[0]].keys()):
        lRowRanks = []
        for strCurSam in sorted(dSamToGeneToRank.keys()):
            lRowRanks.append(dSamToGeneToRank[strCurSam][strCurGene])
            dGeneToMeanRank[strCurGene] = np.average(lRowRanks)

    return dSamToGeneToRank,dGeneToMeanRank

def procGeneCountMatrix(strGeneCount,dSigToGenes,lSigs,strOutFile,strVersion, bIsLog,sigNames):
    dColIDToColLbl = {}
    dSamToGeneToCount = OrderedDict()
    dSamToSigToLVals = OrderedDict()
    dSamToCountToGene = OrderedDict()
    lGlobalAvg = []

    #load Sample To Gene To Count
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
            #processing when more than one gene is in the name (with //)
            strCurGenes = lcols[0].upper().split('//')
            strCurGenes = [strCurGene.strip() for strCurGene in strCurGenes]
            if all(float(val) == 0 for val in lcols[1:]):
                continue
            for strCurGene in strCurGenes:
                for i in range(1,len(lcols)):
                    if strCurGene not in dSamToGeneToCount[dColIDToColLbl[i]]:
                        dSamToGeneToCount[dColIDToColLbl[i]][strCurGene] = float(lcols[i])
                        if dSamToCountToGene[dColIDToColLbl[i]].has_key(float(lcols[i])) == False:
                            dSamToCountToGene[dColIDToColLbl[i]][float(lcols[i])] = []
                        dSamToCountToGene[dColIDToColLbl[i]][float(lcols[i])].append(strCurGene)
                    else:
                        dSamToGeneToCount[dColIDToColLbl[i]][strCurGene] = max(float(lcols[i]),dSamToGeneToCount[dColIDToColLbl[i]][strCurGene])
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
                        continue
                    if bIsLog:
                        dSamToSigToLVals[strSam][strSig].append(np.log10(dSamToGeneToCount[strSam][strCurSigGene]+1))
                    else:
                        dSamToSigToLVals[strSam][strSig].append(dSamToGeneToCount[strSam][strCurSigGene])
        #write outputs
        util.writeDetailedOutput(dSigToGenes,dSamToGeneToCount,strOutFile+'.full.txt',sigNames)
        util.writeRegularOutput(dSamToSigToLVals,strOutFile,sigNames)
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
                    if strVersion == 'rank_avg':
                        dSamToSigToLVals[strSam][strSig].append(dSamToGeneToRank[strSam][strCurSigGene])
                        lGlobalAvg.append(dSamToGeneToRank[strSam][strCurSigGene])
                    elif strVersion == 'rank_delta':
                        dSamToSigToLVals[strSam][strSig].append(dSamToGeneToRank[strSam][strCurSigGene] - dGeneToMeanRank[strCurSigGene])
                        lGlobalAvg.append(dSamToGeneToRank[strSam][strCurSigGene] - dGeneToMeanRank[strCurSigGene])
        util.writeDetailedOutput(dSigToGenes,dSamToGeneToRank,strOutFile+'.full.txt',sigNames)
        util.writeRegularOutput(dSamToSigToLVals,strOutFile,sigNames)


def writeNullModelHists(filename,sigNames,allValues,n,num_iter=100000,num_buckets=1000):
    """Writes each of the histograms to a pdf
    sigNames[i] should correspond with allValues[i]"""
    pdf = PdfPages(filename)
    for i in range(len(allValues)):
        nullmodel.getStatistics(allValues[i],n,pdf,sigNames[i],num_iter,num_buckets)
    pdf.close()

def writeNull(inputFile,version,nullFilename,n):
    """Writes the null distribution of inputFile to nullFilename
    version is either rank_avg, log, or val
    n is the number of genes we're averaging over
    """
    sams = util.readMatrix(inputFile)
    if version == 'rank_avg':
        #technically, this isn't quite accurate, since there could be some ties and this doesn't
        #   include them
        #consider redoing this
        names = ['All Samples']
        allVals = [list(range(1,num_genes))]
    else:
        names = [sam for sam in sorted(sams.keys())]
            allVals = [[sams[sam][gene] for gene in sams[sam]] for sam in names]
            if version == 'log':
                allVals = [[log(val) for val in line] for line in allVals]
    writeNullModelHists(nullFilename,names,allVals,n)


def generateHeatmap(inputFile,sigfile,abbrevs,n,version,zTransform,jobID,rowMetric,colMetric,invert,\
                    fixed,computeNull,isClient):
    """Main function which generates the heatmap
    inputFile - user-provided input
    sigfile - file with most important genes for each signature
        format: sigName\tgene1\tgene2...\nsigName2\t...
    abbrevs - file with abbreviations\tfull names
    n - number of genes to take from each signature
    version - what computation to do
        "log" - log-transformed values
        "val" - values
        "rank_avg" - rank average
        "rank_delta" - rank delta
    zTransform - "matrix" if going to be transformed, otherwise "none"
    jobID - ID of the current job, used for identifying and writing to file
    rowMetric, colMetric - "euclidean", "pearson", "none", used to cluster columns/rows
    invert - invert the axes of the heatmap
    fixed - fix the color values across multiple iterations
    computeNull - compute the null distribution (only works for rank average, value, log)
    isClient - always false when run on server (debug option)
    """
    #matplotlib imports - need to use a writeable directory
    if not options.client:
        os.environ["MPLCONFIGDIR"] = "/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space"
    import matplotlib
    matplotlib.use('Agg')
    import nullmodel
    from matplotlib.backends.backend_pdf import PdfPages

    #checks for errors and corrects whatever errors it can
    util.reformatFile(inputFile)
    util.checkForErrors(inputFile)

    #set up filenames
    compOutput = '/home/mike/workspace/PellegriniResearch/scripts/scratch/output.txt' if isClient\
       else '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_{0}/{0}.matrix.txt'.format(jobID)
    nullFilename = '/home/mike/workspace/PellegriniResearch/output/nulldist.pdf' if isClient\
       else '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/nulldist_' + jobID + '.pdf'
    RHeatmapOut ='/home/mike/workspace/PellegriniResearch/scripts/scratch/Rheatmap.pdf' if isClient\
        else '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/goTeles_tissueDeconvolution_{0}/{0}Rheatmap.pdf'.format(jobID)

    #computing null distribution
    if computeNull and version in ['rank_avg','log','val']:
        writeNull(inputFile,version,nullFilename,n)
    else:
        nullFilename = None

    abbrevsDict = loadAbbrevs(abbrevs)
    sigGenes = getSigGenes(sigfile,abbrevsDict.keys(),n) #TODO: make it filter by selected sigs
    writeValues(inputFile,sigGenes,compOutput,version,abbrevsDict) #TODO: write
    util.createHeatmap(compOutput,RHeatmapOut,version,zTransform,rowMetric,colMetric,jobID,invert,fixed,isClient,nullFilename)
        
#----------------------------------------------------------------------------
# main function call
#----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--gene_count", dest="gene_count", help="table of gene counts", metavar="TAB", default="")
    parser.add_option("-s", "--sig", dest="sig", help="path to pm ratio list", metavar="LIST", default="")
    parser.add_option("-x","--sigfile",dest="sigfile",help="path to precomputed signatures")
    parser.add_option("-g", "--group", dest="group", help="path to group files", metavar="PATH", default="")
    parser.add_option("-n", "--ngene_count", dest="ngene_count", help="number of genes to use", metavar="NUM", default=S_GENE_COUNT)
    parser.add_option("-v", "--version", dest="version", help="metric version (rank_avg,rank_delta,log)", metavar="VER", default=S_VERSION)
    parser.add_option("-z", "--zscore", dest="zscore", help="default no column z-score", metavar="ZSCORE", default=S_ZSCORE)
    parser.add_option("-j", "--job_id", dest="job_id", help="job id, for output path", default="-1")
    parser.add_option("-a", "--absolute", dest="bIsLog", action="store_false", default=True, help="no log transformation")
    parser.add_option("-l", "--log_transform", dest="bIsLog", action="store_true", help="log transform")
    parser.add_option("-r", "--row_metric", dest="row_metric", help="metric for clustering rows (samples)", default="pear_cor")
    parser.add_option("-c", "--col_metric", dest="col_metric", help="metric for clustering columns (signatures)", default="pear_cor")
    parser.add_option("-i","--invert",default=True,dest="invert",help="heatmap columns/rows swtiched")
    parser.add_option("-f", "--fixed", dest="fixed", default="none", help="use fixed color axis")
    parser.add_option("-u", "--null", dest="null", help="compute null model")
    parser.add_option("--client", dest="client", action="store_true", default=True, help="running client-side")
    parser.add_option("--server", dest="client", action="store_false", help="running server-side")
    (options, args) = parser.parse_args()
    
    #load categories
    dGroupToLTisDesc = loadGroupInfo(options.group)

    #list signature
    lTis = []
    for strGroup in sorted(dGroupToLTisDesc.keys()):
        for strTis,strDes in dGroupToLTisDesc[strGroup]:
            lTis.append(strTis)

    #load list of genes
    dGroupSigToGene = getSigGenes(options.sigfile,int(options.ngene_count))

    #process gene count matrix
    procGeneCountMatrix(options.gene_count,dGroupSigToGene,lTis,strOutMatrixTxt,options.version,options.bIsLog,util.loadAbbrevs(options.group))

    #generate heatmap pdf
        util.createHeatmap(strOutMatrixTxt,RHeatmapOut,options.version,options.zscore,options.row_metric,options.col_metric,options.job_id,False if options.invert=='none' else True, False if options.fixed == 'none' else True,options.client,nullFilename)
