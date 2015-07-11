"""To use this file, give the directory of all full signature files
and the number of genes that the maximum option needs (250 at the moment)"""
from optparse import OptionParser
import numpy as np
import glob
import os
def readSigFiles(lFiles):
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
    return dTisToFCToGene

def getSigToGene(dTisToFCToGene,nGeneCount):

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
    return dSigToGene

def writeSigGenes(strSigFile,nGeneCount,output):
    lstrGroupsOnly = glob.glob(strSigFile+'/*--*')
    dGroupTisToFCToGene = readSigFiles(lstrGroupsOnly)
    dGroupSigToGene = getSigToGene(dGroupTisToFCToGene,nGeneCount)
    with open(output,'w') as out:
    	for sig in dGroupSigToGene:
            out.write(sig+"\t")
            for gene in dGroupSigToGene[sig][:-1]:
                out.write(gene+"\t")
            out.write(dGroupSigToGene[sig][-1][0]+"\n")


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-s", "--sig", dest="sig", help="path to pm ratio list", metavar="LIST", default="")
    parser.add_option("-n", "--ngene_max", dest="ngene_max", help="number of genes to use", metavar="NUM")
    parser.add_option("-o","--output",dest="output",help="output of dictionary")
    (options, args) = parser.parse_args()
    dGroupSigToGene = writeSigGenes(options.sig,int(options.ngene_max),options.output)
