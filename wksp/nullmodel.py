import random
import os
from optparse import OptionParser
import subprocess
def shuffleInput(inFile,outFile):
    f = open(inFile).read().split('\n')
    header = f[0]
    f = f[1:]
    genes,vals = [],[]
    for line in f:
        gene,val = line[:line.index('\t')],line[line.index('\t')+1:]
        genes.append(gene)
        vals.append(val)
    with open(outFile,'w') as f:
        f.write(header+"\n")
        while genes:
            i = random.randint(0,len(vals)-1)
            f.write(genes[0] + "\t" + vals[i] + "\n")
            del genes[0]
            del vals[i]
def getMinMaxPerSig(txtFile):
    f = open(txtFile).read().split('\n')
    sigs = f[0].split('\t')[1:]
    sigDict = {}
    for sig in sigs:
        sigDict[sig] = (float('+inf'),float('-inf'))
    for line in f[1:]:
        values = line.split('\t')[1:]
        for i in range(len(values)):
            sig, value = sigs[i],float(values[i])
            sigDict[sig] = (min(sigDict[sig][0],value), max(sigDict[sig][1],value))
    return sigDict
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-o", dest="output")
    parser.add_option("-t", dest="txt")
    parser.add_option("-n", dest="n")
    parser.add_option("-w", dest="w")
    (options, args) = parser.parse_args()
    sigMinsMaxes = dict()
    FNULL = open(os.devnull,'w')
    import time; t = time.time()
    for _ in range(int(options.n)):
        print "iteration",_
        shuffleInput(options.input,options.output)
        subprocess.call(["python","/home/mike/workspace/PellegriniResearch/Sig_Avg_Matrix_Derm.Rank.py", "-t",options.output, "-s","/home/mike/workspace/PellegriniResearch/SigGenes.txt", "--group","/home/mike/workspace/PellegriniResearch/SIGS/abbrevs.txt", "-n","50", "-v","rank_avg", "-z","none", "-j","0", "-l", "-r","none", "-c","none", "-i","checked"],stdout=FNULL,stderr=FNULL)
        sigDict = getMinMaxPerSig(options.txt)
        for sig in sigDict:
            if sig in sigMinsMaxes:
                sigMinsMaxes[sig][0].append(sigDict[sig][0])
                sigMinsMaxes[sig][1].append(sigDict[sig][1])
            else:
                sigMinsMaxes[sig] = ([sigDict[sig][0]],[sigDict[sig][1]])
    print time.time()-t,"seconds elapsed"
    with open(options.w,'w') as w:
        for sig in sigMinsMaxes:
            w.write(sig + "\t" + str(sigMinsMaxes[sig]) + "\n")
    import pdb;pdb.set_trace()
