#!/usr/bin/env python

# Runs Spearman or Pearson on an input matrix
# Main hook for an outside program is runCorrelation
# Written by Mike Ambrose, mikeambrose@berkeley.edu
# Feel free to email with any questions!

import os
from optparse import OptionParser
import iSigDBUtilities as util
import math
from collections import OrderedDict
import subprocess
import iSigDBConstants

def corrRank(sigs,sams,strOutFile,version):
    """Uses the spearman/pearson coefficient to evaluate the similarity between signatures and samples
    sigs is a dictionary where each entry has the n genes chosen for analysis
    sams is in the same format as sigDict, but has samples as keys instead of signatures
    geneNames is a list of gene names in the same order as they appear in sigs and sams
    Similar to other functions, writes to strOutFile with the default output
    No value is returned"""
    #stored as a dictionary of (sig,sam):coef
    sigSamCoefficient = OrderedDict()
    #dictionary of signature:gene:list of d values corresponding to (sig,gene,sam)
    samSigCoefficient = OrderedDict()
    fn = getSpearmanVals if version=="spearman" else getPearsonVals
    for sam in sams:
        samSigCoefficient[sam] = {}
        for sig in sigs:
            samSigCoefficient[sam][sig] = fn(sigs[sig],sams[sam])
    util.writeRegularOutput(samSigCoefficient,strOutFile,{})

def decomp(sigs,sams,job_id,C,genes):
    """Runs the matrix decomposition on signatures/samples using DeconRNASeq
    sigs/sams/outFile same as in corrRank"""
    #first we write both matrices to a tab-delimited file
    for d in [sigs[sig] for sig in sigs] + [sams[sam] for sam in sams]:
        for gene in d.keys():
            if gene not in genes:
                del d[gene]
    util.writeRegularOutput(util.invertDict(sigs),C.DECOMP_SIGS)
    util.writeRegularOutput(util.invertDict(sams),C.DECOMP_SAMS)
    #next, we call the R script, which writes to the directory
    FNULL = open(os.devnull,'w')
    #if this is silently failing, it might be something with the rscript - remove the redirect to null
    subprocess.call([C.RSCRIPT,C.DECOMP_SCRIPT_PATH,C.DECOMP_SIGS,C.DECOMP_SAMS,C.COMP_OUTPUT],stdout=FNULL,stderr=FNULL)

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

def getSpearmanGenes(sams,matrix,compType,mName=None,n=50,high=True,low=False):
    """compType determines if all genes will be accessed or only the top genes
    if 'all' is passed in, computes the set of genes which are common to all samples and signatures
    if 'top' is passed in, takes the top n genes from each signature (as determined by the signature
    file) and computes the intersection of all those top genes, and then selects the subset of those
    which are present in all samples and signatures"""
    assert any([high,low])
    #baseDir = '/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/Matrices/topGenes/'
    baseDir = '/home/mike/workspace/PellegriniResearch/sigdir/MATRICES/topGenes/'
    samGenes = set(sams[sams.keys()[0]].keys())
    if compType=='all': #picks all genes in common
        return list(samGenes.intersection(set(matrix[matrix.keys()[0]].keys())))
    elif compType == 'cov':
        k = 1
        with open('{0}var_{1}'.format(baseDir,mName)) as topGenes:
            candGenes = set()
            for line in topGenes:
                if not line:    continue
                candGenes.add(line.split('\t')[0].upper())
                k += 1
                if n != None and k > n:
                    break
        return list(samGenes.intersection(candGenes))
    candGenes = set()
    #generates a set of candidate genes by picking the top n genes for each sig
    files = []
    if high:
        files.append(baseDir+'high_'+mName)
    if low:
        files.append(baseDir+'low_'+mName)
    if compType=='top': #picks the top n
        for fle in files:
            with open(fle) as sigGenes:
                for line in sigGenes:
                    if not line:
                        continue
                    line = line.upper().replace('\n','').split('\t')[1:]
                    candGenes = candGenes.union(set(x.split(',')[0].upper() for x in line[:int(n)]))
        return list(candGenes.intersection(samGenes))
    elif compType=='mag': #picks all genes which are up/downregulated by a factor of n
        for fle in files:
            with open(fle) as sigGenes:
                for line in sigGenes:
                    if line == "\n":
                        continue
                    lineDebug = line
                    line = line.upper().replace('\n','').split('\t')[1:]
                    for pair in line:
                        pair = pair.split(',')
                        try:
                            if (float(pair[1]) > n and high) or (float(pair[1]) < 1.0/n and low):
                                candGenes.add(pair[0])
                        except Exception:
                            util.displayErrorMessage(repr(lineDebug))
        return list(candGenes.intersection(samGenes))
    else:
        util.displayErrorMessage("Not a valid spearman gene selector " + str(compType))

def getSpearmanDict(inputDict,genes,ordered=False):
    """Contructs a matrix of key : values for each gene"""
    returnDict = {} if not ordered else OrderedDict()
    for line in inputDict:
        returnDict[line] = []
        for gene in genes:
            returnDict[line].append(inputDict[line][gene])
    for line in inputDict:
        assert len(returnDict[line]) == len(genes)
    return returnDict

def runCorrelation(inputFile,version,invert,mn,mx,rowMetric,colMetric,geneMetric,geneVal,selMatrix,\
                    isClient,jobID,low):
    """Runs Spearman or Pearson correlation on the inputFile
    inputFile - user-uploaded file
    version - either "pearson" or "spearman"
    invert - invert the output heatmap
    rowMetric, colMetric - "euclidean", "pearson", or "none"; how output is clustered
    geneMetric - "all", "top", or "mag"; how underlying genes are chosen
    geneVal - corresponding value to go along with geneMetric
    sigMatrix - which matrix is selected
    isClient- debug flag, always False when run on server
    jobID - number generated to uniquely identify task
    low - whether or not to include low genes
    """
    #checks for errors and corrects whichever errors it can
    util.reformatFile(inputFile)
    util.checkForErrors(inputFile)

    #get constants
    C = iSigDBConstants.Constants(jobID)

    #read matrix and sample files
    matrix = util.readMatrix(selMatrix,True,False)
    sams = util.readMatrix(inputFile,True,True)
    #find set of genes
    genes = getSpearmanGenes(sams,matrix,geneMetric,os.path.basename(selMatrix),geneVal,high=True,low=low)
    print "Using {0} genes".format(len(genes))
    if len(genes) == 0:
        util.displayErrorMessage("There are no genes in common between your samples and the matrix selected. Make sure the first column in your input is genes and that they have standard gene names")
    #restrict matrix/samples to only that set of genes
    spearmanSams, spearmanMatrix = getSpearmanDict(sams,genes,True), getSpearmanDict(matrix,genes)
    if version in ["spearman","pearson"]:
        #run correlation on them
        corrRank(spearmanMatrix,spearmanSams,C.COMP_OUTPUT,version)
    else:
        decomp(matrix,sams,jobID,C,genes)
    #pass computation to R/make_heatmap
    util.createHeatmap(C,version,"none",rowMetric,colMetric,invert,isClient,False,False,False,mn=mn,mx=mx)

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
    parser.add_option("--range",dest="range",help="color scale")
    parser.add_option("--low",dest="low",default=False,action="store_true",help="use low genes in addition to high genes")
    op, _ = parser.parse_args()
    mn,mx = (float(x) for x in op.range.split(','))
    runCorrelation(op.inputFile,op.version,op.invert,mn,mx,op.row,op.col,\
                    op.gene,float(op.geneVal),op.selMatrix,True,'',op.low)
