"""This file writes the top genes for each matrix in matrixAbbrevs
For each signature, it writes at least 1000 top genes and writes genes until the val is 2 
It also does the same for the bottom 1000 genes and the val is 1/2
two files are created - one for high genes and one for low
each line in the file is
sig\tgene,val\tgene,val\t...
"""
import glob
matrixAbbrevs = [("Immgen.txt","IMGN_"),("HumanBodyAtlas.txt","HBA_"),("MouseBodyAtlas.txt",""),("MacrophageActivation.txt","MA_")]#,("SkinDiseases.txt","")]
matrixPath = '/home/mike/workspace/PellegriniResearch/sigdir/MATRICES/{matrix}'
highGenePath = '/home/mike/workspace/PellegriniResearch/sigdir/MATRICES/topGenes/high_{matrix}'
lowGenePath = '/home/mike/workspace/PellegriniResearch/sigdir/MATRICES/topGenes/low_{matrix}'
allSignatures = glob.glob('/home/mike/workspace/PellegriniResearch/sigdir/SIGS/*')
allSignatureNames = [sig.split('/')[-1].split('--')[0] for sig in allSignatures]
splitLine = lambda line: (line.split('\t')[0],float(line.split('\t')[1]))

for matrix,abbrev in matrixAbbrevs:
    sigHighGenes = dict() #sig name : gene name : val
    sigLowGenes = dict() 
    with open(matrixPath.format(matrix=matrix)) as m:
        signatures = m.readline().replace('\n','').split('\t')[1:]
    for sig in signatures:
        print "looking at {sig} in {matrix}".format(sig=sig,matrix=matrix)
        sigHighGenes[sig] = {}
        sigLowGenes[sig] = {}
        i = allSignatureNames.index(abbrev+sig)
        geneVals = open(allSignatures[i]).read().split('\n')
        for j in range(len(geneVals)):
            line = geneVals[j]
            if not line:    continue
            gene, val = splitLine(line)
            if val < 2 and j >= 1000:
                break
            sigHighGenes[sig][gene] = val
        for k in range(-1,-len(geneVals)-1,-1):
            line = geneVals[k]
            if not line:    continue
            gene, val = splitLine(line)
            if val > 0.5 and k < -1000:
                break
            sigLowGenes[sig][gene] = val
    #write to files
    with open(highGenePath.format(matrix=matrix),'w') as highFile:
        for sig in sigHighGenes:
            high = sigHighGenes[sig].items()
            high.sort(key=lambda x:-x[1])
            highFile.write("{0}\t{1}\n".format(sig,'\t'.join([str(x)[1:-1].replace('\'','').replace(', ',',') for x in high])))
    with open(lowGenePath.format(matrix=matrix),'w') as lowFile:
        for sig in sigLowGenes:
            low = sigLowGenes[sig].items()
            low.sort(key = lambda x:x[1])
            lowFile.write("{0}\t{1}\n".format(sig,'\t'.join([str(x)[1:-1].replace('\'','').replace(', ',',') for x in low])))
