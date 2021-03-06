import random
from optparse import OptionParser
import matplotlib.pyplot as plt
from math import sqrt,floor

def avN(vals,n):
    """Returns the average of a randomly selected n elements of vals"""
    subset = random.sample(vals,n)
    return sum(subset)/float(n)

def getStatistics(vals,n,pdf,title,num_iter=100000,num_buckets=1000):
    """Returns the simulated values of the distribution of averaging n values from vals
    Also writes a histogram of the distribution to pdf"""
    avs = []
    for _ in range(num_iter):
        avs.append(avN(vals,n))
    mean = sum(avs)/float(num_iter)
    std = sqrt(sum(x**2 for x in avs)/float(num_iter)-mean**2)
    avs.sort()
    plt.hist(avs,num_buckets)
    plt.title(title)
    pdf.savefig()
    plt.close()
    return avs

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
