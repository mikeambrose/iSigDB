import random
import os
from optparse import OptionParser
import subprocess
import matplotlib.pyplot as plt
from math import e,sqrt
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
    inputs = filter(lambda x:not all(p==0 for p in x), inputs)
    return inputs,samNames

def std(lst):
    mean = sum(lst)/float(len(lst))
    EXsquared = sum(x**2 for x in lst)/float(len(lst))
    return sqrt(EXsquared-mean**2)

def pickN(vals,n):
    subset = random.sample(vals,n)
    av = []
    stdev = []
    mn = []
    mx = []
    for i in range(len(vals[0])):
        vals = [sub[i] for sub in subset]
        av.append(sum(vals)/float(n))
        stdev.append(std(vals))
        mn.append(min(vals))
        mx.append(max(vals))
    return av,stdev,mn,mx

def nullRank(N,n,x):
    """bound on probability that the sum of n uniform(1,N) variables is greater than n*x
    derived from a Chernoff bound"""
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
    import time; t = time.time()
    allVals,samNames = getAllVals(options.input,options.rank)
    avs = {}
    stdevs = {}
    mxs = {}
    mns = {}
    for name in samNames:
        avs[name] = []
        stdevs[name] = []
        mxs[name] = []
        mns[name] = []
    for _ in range(int(options.num_iter)):
        if _ % 10000 == 0:
            print "iteration",_
        av, stdev, mn, mx = pickN(allVals,int(options.n))
        for i in range(len(av)):
            avs[samNames[i]].append(av[i])
            stdevs[samNames[i]].append(stdev[i])
            mxs[samNames[i]].append(mx[i])
            mns[samNames[i]].append(mn[i])
    print time.time()-t,"seconds elapsed"
    import pdb;pdb.set_trace()
