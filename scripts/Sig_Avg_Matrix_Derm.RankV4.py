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

def getSigGenes(sigFile,selectedSigs,n):
    sigGenes = {}
    with open(sigFile) as sigToGenes:
        for line in sigToGenes:
            splitLine = line.strip("\n").split("\t")
            sig = splitLine[0]
            if sig in selectedSigs:
                sigGenes[sig] = splitLine[1:n+1]
    return sigGenes

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

def writeValues(sams,sigGenes,compOutput,version,abbrevsDict):
    """Writes the values after computation by version to compOutput
    inputFile - user-provided input file
    sigGenes - dictionary of signature : list of genes in signature
    compOutput - where to write output
    version - how to process input
    abbrevsDict - dictionary of abbreviation : full name for each signature in sigGenes"""
    util.writeDetailedOutput(sigGenes,sams,compOutput+'.full.txt',abbrevsDict)
    samSigVal = {}
    for sam in sams:
        samSigVal[sam] = {}
        for sig in sigGenes:
            geneVals = []
            for gene in sigGenes[sig]:
                if gene not in sams[sam]:
                    #the gene does not exist in our input, so it is reported as a N/A
                    continue
                geneVals.append(sams[sam][gene])
            samSigVal[sam][sig] = util.average(geneVals)
    if 'delta' in version:
        samSigVal = delta(samSigVal) #TODO: implement delta
    util.writeRegularOutput(samSigVal,compOutput,abbrevsDict)

def writeNullModelHists(filename,sigNames,allValues,n,num_iter=100000,num_buckets=1000):
    """Writes each of the histograms to a pdf
    sigNames[i] should correspond with allValues[i]"""
    pdf = PdfPages(filename)
    for i in range(len(allValues)):
        nullmodel.getStatistics(allValues[i],n,pdf,sigNames[i],num_iter,num_buckets)
    pdf.close()

def writeNull(sams,nullFilename,n):
    """Writes the null distribution of inputFile to nullFilename
    n is the number of genes we're averaging over
    """
    names = [sam for sam in sorted(sams.keys())]
    allVals = [[sams[sam][gene] for gene in sams[sam]] for sam in names]
    writeNullModelHists(nullFilename,names,allVals,n)

def generateHeatmap(inputFile,sigfile,abbrevs,n,version,zTransform,jobID,rowMetric,colMetric,invert,\
                    fixed,computeNull,isClient):
    """Main function which generates the heatmap
    inputFile - user-provided input
    sigfile - file with most important genes for each signature
        format: sigName\tgene1\tgene2...\nsigName2\t...
    abbrevs - file with abbreviations\tfull names
    n - number of genes to take from each signature
    version - what computation to do. looks for certain strings
        "avg" - use average values, or "delta" - difference from average
        "log" - log-transform values, or nothing
        "rank" - rank rather than raw values, or nothing
    zTransform - "matrix", "column", "row" if going to be transformed by those metrics, otherwise "none"
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

    #get sample values
    sams = readMatrix(inputFile,ordered=True)
    if 'rank' in version:
        sams = ranked(sams) #TODO: implement ranked
    if 'log' in version:
        sams = logall(sams) #TODO: implement logall

    #computing null distribution
    if computeNull and not zTransform and 'delta' not in version:
        writeNull(sams,version,nullFilename,n)
    else:
        nullFilename = None

    abbrevsDict = loadAbbrevs(abbrevs)
    sigGenes = getSigGenes(sigfile,abbrevsDict.keys(),n) 
    writeValues(sams,sigGenes,compOutput,version,abbrevsDict) #TODO: write
    util.createHeatmap(compOutput,RHeatmapOut,version,zTransform,rowMetric,colMetric,jobID,invert,fixed,isClient,nullFilename)
        
#----------------------------------------------------------------------------
# main function call
#----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("--input", dest="input", help="user-provided input")
    parser.add_option("--sigfile",dest="sigfile",help="path to precomputed signatures")
    parser.add_option("--abbrev", dest="abbrev", help="path to abbrevs file")
    parser.add_option("--n",dest='n', help="number of genes to use")
    parser.add_option("--version", dest="version", help="metric version with avg, log, rank, delta ", metavar="VER", default=S_VERSION)
    parser.add_option("--zTransform", dest="zTransform", help="how to transform matrix")
    parser.add_option("--job_id", dest="job_id", help="job id, for output path")
    parser.add_option("--row_metric", dest="row_metric", help="metric for clustering rows (samples)")
    parser.add_option("--col_metric", dest="col_metric", help="metric for clustering columns (signatures)")
    parser.add_option("--invert",default=False,dest="invert",action="store_true",help="heatmap columns/rows swtiched")
    parser.add_option("-f", "--fixed",default=False, dest="fixed",action="store_true", help="use fixed color axis")
    parser.add_option("-u", "--null",default=False,action="store_true", dest="null", help="compute null model")
    parser.add_option("--client", dest="client",default=False, action="store_true", help="running client-side")
    (options, args) = parser.parse_args()
    generateHeatmap(options.input,options.sigfile,options.abbrev,int(options.n),options.version,options.zTransform,options.job_id,options.row_metric,options.col_mteric,options.invert,options.fixed,options.null,options.client)
