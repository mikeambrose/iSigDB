import random
from optparse import OptionParser
import matplotlib.pyplot as plt
from math import sqrt,floor

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

def avN(vals,n):
    """Returns the average of a randomly selected n elements of vals"""
    subset = random.sample(vals,n)
    return sum(subset)/float(n)

def getStatistics(vals,n,imageLoc,num_iter=100000,num_buckets=1000):
    """Returns the approximate mean, standard deviation, and 95th, 99th, and 99.9th percentile values
    of the distribution created by averaging n values from vals
    Also writes a histogram of the distribution to imageLoc"""
    avs = []
    for _ in range(num_iter):
        avs.append(avN(vals,n))
    mean = sum(avs)/float(num_iter)
    std = sqrt(sum(x**2 for x in avs)/float(num_iter)-mean**2)
    avs.sort()
    p95 = avs[int(num_iter*0.95)]
    p99 = avs[int(num_iter*0.99)]
    p999 = avs[int(num_iter*0.999)]
    plt.hist(avs,num_buckets)
    plt.savefig(imageLoc,bbox_inches='tight')
    plt.close()
    return mean,std,p95,p99,p999

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-c", dest="column")
    parser.add_option("-n", dest="n")
    parser.add_option("-x", dest="num_iter")
    parser.add_option("-o", dest="imageLoc")
    (options, args) = parser.parse_args()
    import time; t = time.time()
    f = open(options.input).read().split('\n')[1:]
    vals = [float(line.split('\t')[int(options.column)+1]) for line in f]
    print getStatistics(vals,int(options.n),options.imageLoc,int(options.num_iter))
    print time.time()-t,"seconds elapsed"
