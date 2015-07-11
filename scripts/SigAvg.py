"""File which does the actual analysis for the signature visualization
Written/modified by Mike Ambrose, mikeambrose@berkeley.edu
Feel free to email with any questions!
"""
import os
from optparse import OptionParser
import subprocess
import math
import iSigDBUtilities as util
#TODO: uncomment this when on server
#os.environ["MPLCONFIGDIR"] = "/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space"
import matplotlib
matplotlib.use('Agg')
import nullmodel
from matplotlib.backends.backend_pdf import PdfPages
from collections import OrderedDict
def getSigGenes(sigFile,selectedSigs,n):
    """Returns dictionary of signature : [top n genes]
    sigFile is the file with the shortened signatures
    n is the number of genes"""
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
    rankedSams = OrderedDict()
    for sam in sams:
        samvals = [sams[sam][gene] for gene in sams[sam]]
        ranks = util.getRanks(tuple(samvals))
        rankedSams[sam] = {gene:ranks[sams[sam][gene]] for gene in sams[sam]}
    return rankedSams

def logall(sams):
    """for each sample, replace the gene's value with the log of its value"""
    return OrderedDict({sam:{gene:math.log(sams[sam][gene]+1,10) for gene in sams[sam]} for sam in sams})

def delta(sams):
    """for each sample, replaces the gene's value with the difference between its value and the
    mean value across all samples"""
    genes = sams[sams.keys()[0]].keys()
    for gene in genes:
        av = util.average([sams[sam][gene] for sam in sams])
        for sam in sams:
            sams[sam][gene] = sams[sam][gene] - av
    return sams

def writeValues(sams,sigGenes,compOutput,version,abbrevsDict):
    """Writes the values after computation by version to compOutput
    inputFile - user-provided input file
    sigGenes - dictionary of signature : list of genes in signature
    compOutput - where to write output
    version - how to process input
    abbrevsDict - dictionary of abbreviation : full name for each signature in sigGenes"""
    samSigVal = OrderedDict()
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
    util.writeRegularOutput(samSigVal,compOutput,abbrevsDict)
    util.writeDetailedOutput(sigGenes,sams,compOutput+'.full.txt',abbrevsDict)

def writeNullModelHists(filename,sigNames,allValues,n,num_iter=100000,num_buckets=100):
    """Writes each of the histograms to a pdf
    sigNames[i] should correspond with allValues[i]"""
    pdf = PdfPages(filename) #version-dependent, but this seems to work for all versions
    for i in range(len(allValues)):
        nullmodel.getStatistics(allValues[i],n,pdf,sigNames[i],num_iter,num_buckets)
    pdf.close()

def writeInputFileHist(filename,sigNames,allValues,num_buckets=100):
    """Writes the distribution of each input signature
    sigNames is the names of the signatures, allValues is the set of values in the signature
    sigNames[i] corresponds to allValues[i]"""
    pdf = PdfPages(filename)
    for i in range(len(allValues)):
        matplotlib.pyplot.hist(allValues[i],num_buckets)
        matplotlib.pyplot.title(sigNames[i])
        pdf.savefig()
        matplotlib.pyplot.close()
    pdf.close()

def writeNull(sams,nullFilename,n,numIter):
    """Writes the null distribution of inputFile to nullFilename
    n is the number of genes we're averaging over
    """
    names = [sam for sam in sorted(sams.keys())]
    allVals = [[sams[sam][gene] for gene in sams[sam]] for sam in names]
    writeNullModelHists(nullFilename,names,allVals,n,numIter)

def writeInDist(sams,inDistFilename):
    """Writes the distribution of the input samples to inDistFilename"""
    names = [sam for sam in sorted(sams.keys())]
    allVals = [[sams[sam][gene] for gene in sams[sam]] for sam in names]
    writeInputFileHist(inDistFilename,names,allVals)

def generateHeatmap(inputFile,sigfile,abbrevs,n,version,zTransform,jobID,rowMetric,colMetric,invert,\
                    fixed,computeNull,isClient,nullIterations):
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
    nullIterations - number of iterations to run the simulation to find null distribution
    """
    #checks for errors and corrects whatever errors it can
    util.reformatFile(inputFile)
    util.checkForErrors(inputFile)

    #set up filenames
    compOutput = '/home/mike/workspace/PellegriniResearch/scripts/scratch/output.txt' if isClient\
       else '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_{0}/{0}.matrix.txt'.format(jobID)
    nullFilename = '/home/mike/workspace/PellegriniResearch/output/nulldist.pdf' if isClient\
       else '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/nulldist_' + jobID + '.pdf'
    inpHistFilename = '/home/mike/workspace/PellegriniResearch/output/inputdist.pdf' if isClient\
       else '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/inputdist_' + jobID + '.pdf'
    RHeatmapOut ='/home/mike/workspace/PellegriniResearch/output/Rheatmap.pdf' if isClient\
        else '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/{0}Rheatmap.pdf'.format(jobID)

    #get sample values
    sams = util.readMatrix(inputFile,ordered=True)

    if 'rank' in version:
        sams = ranked(sams)
    if 'log' in version:
        sams = logall(sams)
    if 'delta' in version:
        sams = delta(sams)

    #computing sample distribution
    writeInDist(sams,inpHistFilename)

    #computing null distribution
    if computeNull and zTransform != 'matrix' and 'delta' not in version:
        writeNull(sams,nullFilename,n,nullIterations)
    else:
        nullFilename = None

    abbrevsDict = util.loadAbbrevs(abbrevs)
    sigGenes = getSigGenes(sigfile,abbrevsDict.keys(),n) 
    writeValues(sams,sigGenes,compOutput,version,abbrevsDict)
    util.createHeatmap(compOutput,RHeatmapOut,version,zTransform,rowMetric,colMetric,jobID,invert,fixed,isClient,nullFilename,inpHistFilename)
        
#----------------------------------------------------------------------------
# main function call
#----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("--input", dest="input", help="user-provided input")
    parser.add_option("--sigfile",dest="sigfile",help="path to precomputed signatures")
    parser.add_option("--abbrev", dest="abbrev", help="path to abbrevs file")
    parser.add_option("--n",dest='n', help="number of genes to use")
    parser.add_option("--version", dest="version", help="metric version with avg, log, rank, delta ")
    parser.add_option("--zTransform", dest="zTransform", help="how to transform matrix")
    parser.add_option("--row_metric", dest="row_metric", help="metric for clustering rows (samples)")
    parser.add_option("--col_metric", dest="col_metric", help="metric for clustering columns (signatures)")
    parser.add_option("--invert",default=False,dest="invert",action="store_true",help="heatmap columns/rows swtiched")
    parser.add_option("--fixed",default=False, dest="fixed",action="store_true", help="use fixed color axis")
    parser.add_option("--null",default=False,action="store_true", dest="null", help="compute null model")
    parser.add_option("--nullIter",dest="nullIterations",help="how many iterations to use in null model")
    (options, args) = parser.parse_args()
    generateHeatmap(options.input,options.sigfile,options.abbrev,int(options.n),options.version,options.zTransform,None,options.row_metric,options.col_metric,options.invert,options.fixed,options.null,True,options.nullIterations if options.null else None)
