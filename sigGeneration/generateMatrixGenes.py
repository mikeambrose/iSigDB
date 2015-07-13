"""This file writes the top genes for each matrix in matrixAbbrevs
For each signature, it writes at least 1000 top genes and writes genes until the val is 2 
It also does the same for the bottom 1000 genes and the val is 1/2
two files are created - one for high genes and one for low
each line in the file is
sig\tgene,val\tgene,val\t...
"""
import glob
import math
matrixAbbrevs = [("Immgen.txt","IMGN_"),("HumanBodyAtlas.txt","HBA_"),("MouseBodyAtlas.txt",""),("MacrophageActivation.txt","MA_"),("SkinDiseases.txt","")]
matrixPath = '/home/mike/workspace/PellegriniResearch/sigdir/MATRICES/{matrix}'
highGenePath = '/home/mike/workspace/PellegriniResearch/sigdir/MATRICES/topGenes/high_{matrix}'
lowGenePath = '/home/mike/workspace/PellegriniResearch/sigdir/MATRICES/topGenes/low_{matrix}'
varGenePath = '/home/mike/workspace/PellegriniResearch/sigdir/MATRICES/topGenes/var_{matrix}'
allSignatures = glob.glob('/home/mike/workspace/PellegriniResearch/sigdir/SIGS/*')
allSignatureNames = [sig.split('/')[-1].split('--')[0] for sig in allSignatures]
splitLine = lambda line: (line.split('\t')[0],float(line.split('\t')[1]))
cov = lambda vals: stdev(vals,av(vals)) / av(vals)
av = lambda lst: sum(lst) / float(len(lst))
stdev = lambda lst,av: math.sqrt(sum((x-av)**2 for x in lst) / len(lst))

NUM_PM_GENES = 1000
THRESHOLD_PMS = 2,0.5
NUM_COV_GENES = 5000

for matrix,abbrev in matrixAbbrevs:

    print "Calculating high/low genes based on signature values for {0}".format(matrix)
    sigHighGenes = dict() #sig name : gene name : val
    sigLowGenes = dict() 
    with open(matrixPath.format(matrix=matrix)) as m:
        signatures = m.readline().replace('\n','').split('\t')[1:]
    for sig in signatures:
        print "looking at signature {0}".format(sig)
        sigHighGenes[sig] = {}
        sigLowGenes[sig] = {}
        try:
            i = allSignatureNames.index(abbrev+sig)
        except Exception:
            print "failure to find signature {sig} in {matrix}".format(sig=sig,matrix=matrix)
            exit()
        geneVals = open(allSignatures[i]).read().upper().split('\n')
        for j in range(len(geneVals)):
            line = geneVals[j]
            if not line:    continue
            gene, val = splitLine(line)
            if val < THRESHOLD_PMS[0] and j >= NUM_PM_GENES:
                break
            sigHighGenes[sig][gene] = val
        for k in range(-1,-len(geneVals)-1,-1):
            line = geneVals[k]
            if not line:    continue
            gene, val = splitLine(line)
            if val > THRESHOLD_PMS[1] and k < -NUM_PM_GENES:
                break
            sigLowGenes[sig][gene] = val
    
    print "calculating genes with high coefficient of variation for {0}".format(matrix)
    geneCovs = {}
    with open(matrixPath.format(matrix=matrix)) as m:
        m.readline() #skip header
        for line in m:
            if 'N/A' in line:   continue
            line = line.replace('\n','').upper().split('\t')
            gene, vals = line[0],[float(x) for x in line[1:]]
            geneCovs[gene] = cov(vals)
    geneCovs = geneCovs.items()
    geneCovs.sort(key=lambda x:-x[1])

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
    
    with open(varGenePath.format(matrix=matrix),'w') as covFile:
        for gene,c in geneCovs[:NUM_COV_GENES]:
            covFile.write("{gene}\t{c}\n".format(gene=gene,c=c))
