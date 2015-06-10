from math import log
import os.path, glob
def logSigs(inputDir,outputDir):
    allSigs = glob.glob(inputDir+'/*--*')
    allSigs = filter(lambda x:'pm' in x.lower() and 'log' not in x.lower(),allSigs)
    
    for sig in allSigs:
        geneVals = {}
        inp = open(sig).read().split('\n')
        for line in inp:
            if line == '':
                continue
            line = line.split('\t')
            gene = line[0]
            val = float(line[1])
            try:
                geneVals[gene] = abs(log(val,10)) if val != 0 else 0
            except:
                print gene, val, sig
                exit()
        with open(outputDir+'/'+os.path.basename(sig),'w') as out:
            geneToVal = sorted(geneVals.items(),key=lambda x:x[1])
            for pair in geneToVal:
                out.write(pair[0]+"\t"+str(pair[1])+"\n")

logSigs('SIGS','ABSLOGSIGS')
