#!/usr/bin/env python

# Runs Spearman or Pearson on an input matrix
# Main hook for an outside program is runCorrelation
# Written by Mike Ambrose, mikeambrose@berkeley.edu
# Feel free to email with any questions!

import os
from optparse import OptionParser
import iSigDBUtilities as util
import math

def corrRank(sigs,sams,geneNames,strOutFile,sigNames,version):
    """Uses the spearman/pearson coefficient to evaluate the similarity between signatures and samples
    sigs is a dictionary where each entry has the n genes chosen for analysis
    sams is in the same format as sigDict, but has samples as keys instead of signatures
    geneNames is a list of gene names in the same order as they appear in sigs and sams
    sigNames is a dictionary of sig abbreviation : full name
    Similar to other functions, writes to strOutFile with the default output
    No value is returned"""
    #stored as a dictionary of (sig,sam):coef
    sigSamCoefficient = {}
    #dictionary of signature:gene:list of d values corresponding to (sig,gene,sam)
    for sig in sigs:
        sigSamCoefficient[sig] = {}
        for sam in sams:
            coef = getSpearmanVals(sigs[sig],sams[sam]) if version=="spearman"\
                                    else getPearsonVals(sigs[sig],sams[sam])
            sigSamCoefficient[(sig,sam)] = coef
    #write regular output
    samSigCoefficient = {}
    for sig in sigs:
        samSigCoefficient[sig] = {sam:sigSamCoefficient[(sig,sam)] for sam in sams}
    util.writeRegularOutput(samSigCoefficient,strOutFile,sigNames)


def getSpearmanVals(v1,v2):
    """Returns the spearman correlation coefficient between v1 and v2"""
    n = len(v1)
    v1Ranks,v2Ranks = util.getRanks(tuple(v1)),util.getRanks(tuple(v2))
    d = [v1Ranks[v1[i]] - v2Ranks[v2[i]] for i in range(n)]
    rho = 1 - 6*sum(x**2 for x in d)/float(n*(n**2-1))
    return rho

def getPearsonVals(v1,v2):
    """Returns the pearson correlation coefficient between v1 and v2"""
    n = len(v1)
    d = [v1[i] - v2[i] for i in range(n)]
    rho = (n * sum(v1[i]*v2[i] for i in range(n)) - sum(v1)*sum(v2))/\
            math.sqrt((n*sum(x**2 for x in v1) - sum(v1)**2)*(n*sum(x**2 for x in v2)-sum(v2)**2))
    return rho

def getCommonGenes(candGenes,allSets):
    """Returns the subset of candGenes which occur in every member of allSets"""
    commonGenes = set()
    for gene in candGenes:
        if all(gene in s for s in allSets):
            commonGenes.add(gene)
    return commonGenes

#TODO: rethink this with the value-matrix approach
def getSpearmanGenes(sams,sigs,compType,sigFile=None,n=50):
    """compType determines if all genes will be accessed or only the top genes
    if 'all' is passed in, computes the set of genes which are common to all samples and signatures
    if 'top' is passed in, takes the top n genes from each signature (as determined by the signature
    file) and computes the intersection of all those top genes, and then selects the subset of those
    which are present in all samples and signatures"""
    if compType=='all': #picks all genes in common
        commonGenes = set()
        candGenes = sams[sams.keys()[0]].keys()
        return list(getCommonGenes(candGenes,[sigs[sig] for sig in sigs]))
    elif compType=='top': #picks the top n
        candGenes = set()
        #generates a set of candidate genes by picking the top n genes for each sig
        with open(sigFile) as fullSigs:
            for line in fullSigs:
                line = line.split('\t')
                sig = line[0]
                if sig in sigs:
                    candGenes = candGenes.union(set(line[1:n+1]))
        allSets = [sams[sams.keys()[0]]] + [sigs[sig] for sig in sigs]
        return list(getCommonGenes(candGenes,allSets))
    elif compType=='mag': #picks all genes which are up/downregulated by a factor of n
        candGenes = set()
        for sig in sigs:
            for gene in sigs[sig]:
                if sigs[sig][gene] >= n:# or sigs[sig][gene] <= 1.0/n:
                    candGenes.add(gene)
        allSets = [sams[sams.keys()[0]]] + [sigs[sig] for sig in sigs]
        return list(getCommonGenes(candGenes,allSets))
    else:
        util.displayErrorMessage("Not a valid spearman gene selector " + str(compType))

def getSpearmanDict(inputDict,genes):
    """Contructs a matrix of key : values for each gene"""
    returnDict = {}
    for line in inputDict:
        returnDict[line] = []
        for gene in genes:
            returnDict[line].append(inputDict[line][gene])
    for line in inputDict:
        assert len(returnDict[line]) == len(genes)
    return returnDict

def runCorrelation(inputFile,version,invert,rowMetric,colMetric,geneMetric,geneVal,selMatrix,\
                    isClient,job_id,abbrevs):
    """Runs Spearman or Pearson correlation on the inputFile
    inputFile - user-uploaded file
    version - either "pearson" or "spearman"
    invert - invert the output heatmap
    rowMetric, colMetric - "euclidean", "pearson", or "none"; how output is clustered
    geneMetric - "all", "top", or "mag"; how underlying genes are chosen
    geneVal - corresponding value to go along with geneMetric
    sigMatrix - which matrix is selected
    isClient- debug flag, always False when run on server
    job_id - number generated to uniquely identify task
    abbrevs - file of abbreviation\tfull name
    """
    #checks for errors and corrects whichever errors it can
    util.reformatFile(inputFile)
    util.checkForErrors(inputFile)

    outFile = '/home/mike/workspace/PellegriniResearch/scripts/scratch/output.txt' if isClient\
            else '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_'\
                                +job_id+'/'+job_id+'.matrix.txt'
    #read matrix and sample files
    matrix = util.readMatrix(selMatrix,True,False)
    sams = util.readMatrix(inputFile,True,True)
    #find set of genes
    genes = getSpearmanGenes(sams,matrix,geneMetric,geneVal)
    if len(genes) == 0:
        util.displayErrorMessage("There are no genes in common between your samples and the matrix selected. Make sure the first column in your input is genes and that they have standard gene names")
    #restrict matrix/samples to only that set of genes
    spearmanSams, spearmanMatrix = getSpearmanDict(sams,genes), getSpearmanDict(matrix,genes)
    #get full names of signatures
    abbrevs = util.loadAbbrevs(abbrevs)
    #run correlation on them
    corrRank(spearmanMatrix,spearmanSams,genes,outFile,abbrevs,version)
    RHeatmapOut ='/home/mike/workspace/PellegriniResearch/scripts/scratch/Rheatmap.pdf' if isClient\
        else '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/goTeles_tissueDeconvolution_{0}/{0}Rheatmap.pdf'.format(options.job_id)
    #pass computation to R/make_heatmap
    util.createHeatmap(outFile,RHeatmapOut,version,"none",rowMetric,colMetric,job_id,invert,True,isClient,None)

#this file should always be imported by the server
#running from the command line is only for testing
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--input",dest="inputFile",help="user-uploaded file")
    parser.add_option("--version",dest="version",help="which correlation to run (spearman or pearson)")
    parser.add_option("--invert",dest="invert",default=False,action="store_true",help="invert the output")
    parser.add_option("--row",dest="row",help="row clustering (euclidean, pearson, none)")
    parser.add_option("--col",dest="col",help="column clustering (euclidean, pearson, none)")
    parser.add_option("--gene",dest="gene",help="gene metric (all, top, mag)")
    parser.add_option("--gval",dest="geneVal",help="value corresponding with gene metric")
    parser.add_option("--matrix",dest="selMatrix",help="matrix selected")
    parser.add_option("--abbrev",dest="abbrevs",help="abbreviation file")
    options, _ = parser.parse_args()
    runCorrelation(options.inputFile,options.version,options.invert,options.row,options.col,\
                    options.gene,options.geneVal,options.selMatrix,True,None,options.abbrevs)
