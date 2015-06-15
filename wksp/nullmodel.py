import random
import os
from optparse import OptionParser
import subprocess
import matplotlib.pyplot as plt
from math import e
def getAllVals(inFile,rank):
    f = open(inFile).read().split('\n')
    samNames = f[0].split('\t')[1:]
    del f[0]
    inputs = []
    if not rank=="checked":
        for line in f:
            inputs.append([float(x) for x in line.split('\t')[1:]])
    else:
        n = 1
        for line in f:
            inputs.append([n])
            n += 1
        samNames = [samNames[0]]
    return inputs,samNames

def pickN(vals,n):
    subset = random.sample(vals,n)
    av = []
    for i in range(len(vals[0])):
        av.append(sum([sub[i] for sub in subset])/float(n))
    return av

def nullRank(N,n,x):
    """probability that the sum of n uniform(1,N) variables is greater than n*x"""
    t = 1.0/(N+1)
    #inner = (e - 1) / ((e**t - 1)*(e**(t*x))) / float(N)
    inner = (e-1)/(N*(e**t-1)*e**(x*t))
    return inner**n


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-n", dest="n")
    parser.add_option("-r", dest="rank")
    parser.add_option("-x", dest="num_iter")
    (options, args) = parser.parse_args()
    sigMinsMaxes = dict()
    FNULL = open(os.devnull,'w')
    import time; t = time.time()
    allVals,samNames = getAllVals(options.input,options.rank)
    samDists = {}
    for name in samNames:
        samDists[name] = []
    for _ in range(int(options.num_iter)):
        if _ % 100 == 0:
            print "iteration",_
        transformedNValues = pickN(allVals,int(options.n))
        for i in range(len(transformedNValues)):
            samDists[samNames[i]].append(transformedNValues[i])
    print time.time()-t,"seconds elapsed"
    import pdb;pdb.set_trace()
