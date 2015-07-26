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

def getTopNGenes(geneVals,n):
    """Get the top n genes in the signature"""
    allGeneVals = geneVals.items()
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

def writeSigGenes(sigDir,n,humanOutput,mouseOutput,transDir):
    allSigs = glob.glob(sigDir+'/*--*')
    mouseToHuman,humanToMouse = mhd.getMouseHumanDicts(transDir+mhd.human,transDir+mhd.mouse,transDir+mhd.both)
    sigGeneVals = readSigFiles(allSigs)
    with open(humanOutput,'w') as hOut:
        with open(mouseOutput,'w') as mOut:
            for sig in sigGeneVals:
                #human output
                humanGeneCount = sum(gene in humanToMouse for gene in sigGeneVals[sig])
                mouseGeneCount = sum(gene in mouseToHuman for gene in sigGeneVals[sig])
                isHuman = humanGeneCount > mouseGeneCount
                print "Classifing signature {0} as {1}".format(sig,"human" if isHuman else "mouse")
                if isHuman:
                    topGenes = getTopNGenes(sigGeneVals[sig],n)
                else:
                    ambiguousGenes = set()
                    mGeneVals = {}
                    for gene in sigGeneVals[sig]:
                        if gene in mouseToHuman:
                            if mouseToHuman[gene] in mGeneVals and mouseToHuman[gene] not in ambiguousGenes:
                                print "conflicting/ambiguous gene {0} in {1} - removing".format(mouseToHuman[gene],sig,gene)
                                ambiguousGenes.add(mouseToHuman[gene])
                            mGeneVals[mouseToHuman[gene]] = sigGeneVals[sig][gene]
                    topGenes = getTopNGenes(mGeneVals,n)
                hOut.write("{sig}\t{genes}\n".format(sig=sig,genes="\t".join(topGenes)))
                #mouse output
                if not isHuman:
                    topGenes = getTopNGenes(sigGeneVals[sig],n)
                else:
                    hGeneVals = {}
                    for gene in sigGeneVals[sig]:
                        if gene in humanToMouse:
                            if humanToMouse[gene] in hGeneVals:
                                print "conflicting/ambiguous gene {0} in {1} - removing".format(mouseToHuman[gene],sig)
                                del hGeneVals[humanToMouse[gene]]
                            hGeneVals[humanToMouse[gene]] = sigGeneVals[sig][gene]
                    topGenes = getTopNGenes(hGeneVals,n)
                mOut.write("{sig}\t{genes}\n".format(sig=sig,genes="\t".join(topGenes)))

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-s", "--sig", dest="sig", help="signature directory")
    parser.add_option("-n", dest="n", help="max number of genes")
    parser.add_option("-d",dest="d",help="mouse-human file directory")
    parser.add_option("-x",dest="humanOutput",help="human output location")
    parser.add_option("-m",dest="mouseOutput",help="mouse output location")
    options,_ = parser.parse_args()
    writeSigGenes(options.sig,int(options.n),options.humanOutput,options.mouseOutput,options.d)
