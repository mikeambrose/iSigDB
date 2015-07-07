import matplotlib.pyplot as plt
import math
import glob

def plotSig(s,l=False):
    vals = []
    with open(s) as sig:
        for line in sig:
            vals.append(float(line.split('\t')[1].replace('\n','')))
    if l:
        vals = [math.log(v,10) for v in vals]
    plt.hist(vals,50000)
    plt.show()


def plotSigs(ss,l=False):
    vals = []
    for s in ss:
        with open(s) as sig:
            for line in sig:
                vals.append(float(line.split('\t')[1].replace('\n','')))
    print sum(v<0.5 for v in vals), len(vals)
    if l:
        vals = [math.log(v,10) for v in vals]
    plt.hist(vals,50000)
    plt.show()
