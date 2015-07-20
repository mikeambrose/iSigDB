"""To use this file, give the directory of all full signature files
and the number of genes that the maximum option needs (1000 at the moment)"""
from optparse import OptionParser
import numpy as np
import glob
import os
import MouseHumanDict as mhd

def readSigFiles(files):
    """Returns a dictionary of signature name : gene : value"""
    sigGeneVals = {}
    for fle in files:
        sigName = fle.split('/')[-1].split('--')[0]
        sigGeneVals[sigName] = {}
        with open(fle) as f:
            for line in f:
                gene,val = line.upper().replace('\n','').replace('\r','').split('\t')
                sigGeneVals[sigName][gene] = float(val)
    return sigGeneVals

def getTopNGenes(sigGeneVals,sig,n):
    """Get the top n genes in the signature"""
    allGeneVals = sigGeneVals[sig].items()
    allGeneVals.sort(key=lambda x:-x[1])
    return [x[0] for x in allGeneVals[:n]]

def geneStr(gene,mouseToHuman,humanToMouse):
    assert not (gene in humanToMouse and gene in mouseToHuman)
    if gene in mouseToHuman:
        return "{0}##{1}".format(mouseToHuman[gene],gene)
    elif gene in humanToMouse:
        return "{0}##{1}".format(gene,humanToMouse[gene])
    else:
        return gene

def writeSigGenes(sigDir,n,output):
    allSigs = glob.glob(sigDir+'/*--*')
    mouseToHuman,humanToMouse = mhd.getMouseHumanDicts(mhd.human,mhd.mouse,mhd.both)
    sigGeneVals = readSigFiles(allSigs)
    questionableHumanGenes = set()
    questionableMouseGenes = set()
    with open(output,'w') as out:
        for sig in sigGeneVals:
            topGenes = getTopNGenes(sigGeneVals,sig,n)
            #nonessential - checking how many genes for each signature are 'human' and 'mouse'
            mouseGenes = sum(g in mouseToHuman for g in topGenes)
            humanGenes = sum(g in humanToMouse for g in topGenes)
            if mouseGenes != 0 and humanGenes != 0:
                if humanGenes > mouseGenes:
                    for gene in mouseToHuman:
                        if gene in topGenes:
                            questionableHumanGenes.add(gene)
                if mouseGenes > humanGenes:
                    for gene in humanToMouse:
                        if gene in topGenes:
                            questionableMouseGenes.add(gene)
            #end nonessential section
            out.write("{sig}\t{genes}\n".format(sig=sig,genes="\t".join(\
                                    geneStr(g,mouseToHuman,humanToMouse) for g in topGenes)))
    print "'mouse' genes found in human signatures: ",questionableHumanGenes
    print "'human' genes found in mouse signatures: ",questionableMouseGenes

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-s", "--sig", dest="sig", help="signature directory")
    parser.add_option("-n", dest="n", help="max number of genes")
    parser.add_option("-o","--output",dest="output",help="output dictionary location")
    options,_ = parser.parse_args()
    writeSigGenes(options.sig,int(options.n),options.output)
