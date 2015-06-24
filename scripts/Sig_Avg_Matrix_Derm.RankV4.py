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

def ranked(sams):
    """for each sample, replaces the gene's value with its rank"""
    rankedSams = {}
    for sam in sams:
        samvals = [sams[sam][gene] for gene in sams[sam]]
        ranks = util.getRanks(tuple(samvals))
        rankedSams[sam] = {gene:ranks[sams[sam][gene]] for gene in sams[sam]}
    return rankedSams

def logall(sams):
    """for each sample, replace the gene's value with the log of its value"""
    return {sam:{gene:math.log(sams[sam][gene]+1,10) for gene in sams[sam]} for sam in sams}

def delta(samSigVals):
    """for each signature, replace each samples values with the difference between their value
    and the average value across all samples"""
    for sig in samSigVals[samSigVals.keys()[0]]:
        av = util.average([samSigVals[sam][sig] for sam in samSigVals])
        for sam in samSigVals:
            samSigVals[sam][sig] = samSigVals[sam][sig] - av
    return samSigVals

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
    if not isClient:
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
    RHeatmapOut ='/home/mike/workspace/PellegriniResearch/output/Rheatmap.pdf' if isClient\
        else '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/goTeles_tissueDeconvolution_{0}/{0}Rheatmap.pdf'.format(jobID)

    #get sample values
    sams = util.readMatrix(inputFile,ordered=True)
    if 'rank' in version:
        sams = ranked(sams) #TODO: implement ranked
    if 'log' in version:
        sams = logall(sams) #TODO: implement logall

    #computing null distribution
    if computeNull and not zTransform and 'delta' not in version:
        writeNull(sams,version,nullFilename,n)
    else:
        nullFilename = None

    abbrevsDict = util.loadAbbrevs(abbrevs)
    sigGenes = getSigGenes(sigfile,abbrevsDict.keys(),n) 
    writeValues(sams,sigGenes,compOutput,version,abbrevsDict) #TODO: write
    print RHeatmapOut
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
    parser.add_option("--row_metric", dest="row_metric", help="metric for clustering rows (samples)")
    parser.add_option("--col_metric", dest="col_metric", help="metric for clustering columns (signatures)")
    parser.add_option("--invert",default=False,dest="invert",action="store_true",help="heatmap columns/rows swtiched")
    parser.add_option("--fixed",default=False, dest="fixed",action="store_true", help="use fixed color axis")
    parser.add_option("--null",default=False,action="store_true", dest="null", help="compute null model")
    (options, args) = parser.parse_args()
    generateHeatmap(options.input,options.sigfile,options.abbrev,int(options.n),options.version,options.zTransform,None,options.row_metric,options.col_metric,options.invert,options.fixed,options.null,True)
