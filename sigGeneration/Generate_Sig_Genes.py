"""To use this file, give the directory of all full signature files
and the number of genes that the maximum option needs (1000 at the moment)"""
from optparse import OptionParser
import numpy as np
import glob
import os
"""def readSigFiles(lFiles):
    dTisToFCToGene = {}
    for strFile in lFiles:
        head,tail = os.path.split(strFile)
        strTis =  tail.split('--')[0]
        dTisToFCToGene[strTis]={}
        f = open(strFile,'r')
        for strLine in f:
            strGene,strRatio = strLine.rstrip().split('\t')
            strGene = strGene.upper()
            if strRatio == '0':
                continue
            flogFC = float(strRatio)
            if not dTisToFCToGene[strTis].has_key(flogFC):
                dTisToFCToGene[strTis][flogFC] = []
            dTisToFCToGene[strTis][flogFC].append(strGene)
        f.close()
    return dTisToFCToGene"""

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

"""def getSigToGene(dTisToFCToGene,nGeneCount):
    dSigToGene = {}
    #load UP signatures
    for strTis in dTisToFCToGene.keys():
        dSigToGene[strTis] = []
        nCurGeneCount = 0
        for fFC in sorted(dTisToFCToGene[strTis].keys(),reverse=True):
            for strCurGene in dTisToFCToGene[strTis][fFC]:
                nCurGeneCount+=1
                dSigToGene[strTis].append(strCurGene)
                if nCurGeneCount >= nGeneCount:
                    break
            if nCurGeneCount >= nGeneCount:
                break
    return dSigToGene"""

def writeSigGenes(sigDir,n,output):
    allSigs = glob.glob(sigDir+'/*--*')
    sigGeneVals = readSigFiles(allSigs)
    with open(output,'w') as out:
        for sig in sigGeneVals:
            topGenes = getTopNGenes(sigGeneVals,sig,n)
            out.write("{sig}\t{genes}\n".format(sig=sig,genes="\t".join(topGenes)))

"""def writeSigGenes(strSigFile,nGeneCount,output):
    allSigs = glob.glob(strSigFile+'/*--*')
    dGroupTisToFCToGene = readSigFiles(lstrGroupsOnly)
    dGroupSigToGene = getSigToGene(dGroupTisToFCToGene,nGeneCount)
    with open(output,'w') as out:
    	for sig in dGroupSigToGene:
            out.write(sig+"\t")
            for gene in dGroupSigToGene[sig][:-1]:
                out.write(gene+"\t")
            out.write(dGroupSigToGene[sig][-1][0]+"\n")"""


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-s", "--sig", dest="sig", help="signature directory")
    parser.add_option("-n", dest="n", help="max number of genes")
    parser.add_option("-o","--output",dest="output",help="output dictionary location")
    options,_ = parser.parse_args()
    writeSigGenes(options.sig,int(options.n),options.output)
