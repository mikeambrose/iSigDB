import glob,os,subprocess
from optparse import OptionParser

def spearmanRank(sigs,sams,geneNames,strOutFile,sigNames,pearson=False):
    """Uses the spearman rank coefficient to evaluate the similarity between signatures and samples
    sigs is a dictionary where each entry has the n genes chosen for analysis
    sams is in the same format as sigDict, but has samples as keys instead of signatures
    geneNames is a list of gene names in the same order as they appear in sigs and sams
    sigNames is a dictionary of sig abbreviation : full name
    Similar to other functions, writes to two files:
        strOutFile contains the basic information used by the R script
        strOutFile + ".matrix.txt" has the detailed information used in generating the detailed heatmaps
    No value is returned"""
    #stored as a dictionary of (sig,sam):coef
    sigSamCoefficient = {}
    #dictionary of signature:gene:list of d values corresponding to (sig,gene,sam)
    geneDetail = {}
    for sig in sigs:
        geneDetail[sig] = {}
        for gene in geneNames:
            geneDetail[sig][gene] = []
    for sig in sigs:
        for sam in sams:
            coef,rankDiffs = getSpearmanVals(sigs[sig],sams[sam]) if not pearson\
                                    else getPearsonVals(sigs[sig],sams[sam])
            for i in range(len(geneNames)):
                geneDetail[sig][geneNames[i]].append(rankDiffs[i])
            sigSamCoefficient[(sig,sam)] = coef
    #write detailed output
    with open(strOutFile+".matrix.txt","w") as detailedOut:
        for sig in sigs:
            detailedOut.write('# Signature: ' + (sigNames[sig] if sig in sigNames else sig)+'\n')
            #write sample names
            for sam in sams:
                detailedOut.write('\t'+sam)
            detailedOut.write('\n')
            #write gene detail
            for gene in geneNames:
                detailedOut.write(gene)
                samVals = geneDetail[sig][gene]
                assert len(samVals) == len(sams)
                for val in samVals:
                    detailedOut.write('\t'+str(val))
                detailedOut.write('\n')
            detailedOut.write('\n')
    #write regular output
    with open(strOutFile,"w") as basicOut:
        #write header
        basicOut.write("SAMPLE")
        for sig in sigs:
            #use the full name if it exists, otherwise use the abbreviation
            basicOut.write("\t"+(sigNames[sig] if sig in sigNames else sig))
        basicOut.write("\n")
        for sam in sams:
            basicOut.write(sam)
            for sig in sigs:
                basicOut.write("\t"+str(sigSamCoefficient[(sig,sam)]))
            basicOut.write("\n")


def getRanksMemoize(f):
    """fast memoization decorator for getRanks"""
    class memodict(dict):
        def __missing__(self,key):
            ret = self[key] = f(key)
            return ret
    return memodict().__getitem__

@getRanksMemoize
def getRanks(lst):
    """Returns a dictionary of entry : rank
    This is memoized and thus must be called with a tuple"""
    ranks = {}
    sorted_lst = sorted(lst)
    n = len(lst)
    i = 0
    while i < n:
        #detect if there is a run
        if i+1<n and sorted_lst[i+1] == sorted_lst[i]:
            run_length = 2
            j = i+1
            while j+1<n and sorted_lst[j+1] == sorted_lst[j]:
                run_length += 1
                j += 1
            rankVal = i+(run_length-1)/2.0
            for k in range(i,j+1):
                ranks[sorted_lst[k]] = rankVal
            i = i+run_length
        else:
            ranks[sorted_lst[i]] = i
            i += 1
    return ranks

def getSpearmanVals(v1,v2):
    """Returns the spearman correlation coefficient between v1 and v2
    as well as the list of rank differentials between v1 and v2"""
    n = len(v1)
    v1Ranks,v2Ranks = getRanks(tuple(v1)),getRanks(tuple(v2))
    d = [v1Ranks[v1[i]] - v2Ranks[v2[i]] for i in range(n)]
    #import pdb;pdb.set_trace()
    rho = 1 - 6*sum(x**2 for x in d)/float(n*(n**2-1))
    return rho,d
import math
def getPearsonVals(v1,v2):
    """Returns the pearson correlation coefficient between v1 and v2 as well as the list
    of differentials between v1 and v2"""
    n = len(v1)
    d = [v1[i] - v2[i] for i in range(n)]
    rho = (n * sum(v1[i]*v2[i] for i in range(n)) - sum(v1)*sum(v2))/\
            math.sqrt((n*sum(x**2 for x in v1) - sum(v1)**2)*(n*sum(x**2 for x in v2)-sum(v2)**2))
    return rho,d

def getLogPearsonVals(v1,v2):
    """runs getPearsonVals on log(v1) and log(v2)"""
    return getPearsonVals([math.log(x) for x in v1],[math.log(x) for x in v2])
    
def getSpearmanGenes(sams,sigs,compType,n=50):
    if compType=='all':
        commonGenes = set()
        sam = sams[sams.keys()[0]]
        for gene in sam:
            if all(gene in sigs[sig] for sig in sigs):
                commonGenes.add(gene)
        print len(commonGenes),'genes'
        return list(commonGenes)
    elif compType=='top':
        candGenes = set()
        #generates a set of candidate genes by picking the top n genes for each sig
        with open('SigGenes.txt') as fullSigs:
            for line in fullSigs:
                line = line.split('\t')
                sig = line[0]
                if sig in sigs:
                    candGenes = candGenes.union(set(line[1:n+1]))
        commonGenes = set()
        allSets = [sams[sams.keys()[0]]] + [sigs[sig] for sig in sigs]
        for gene in candGenes:
            if all(gene in s for s in allSets):
                commonGenes.add(gene)
        print len(commonGenes),'genes'
        return list(commonGenes)

def getSpearmanDict(inputDict,genes):
    """Contructs a matrix of key : values for each gene"""
    returnDict = {}
    for line in inputDict:
        returnDict[line] = []
        for gene in genes:
            returnDict[line].append(inputDict[line][gene])
    for line in inputDict:
        assert len(returnDict[line]) == len(genes)
    return returnDict

def getSamDict(f):
    """returns a dictionary of sample:gene:value for each gene in genes"""
    samDict = {}
    f = open(f).read().split('\n')
    sams = f[0].replace("\n","").replace("\r","").split('\t')[1:]
    for sam in sams:
        samDict[sam] = {}
    for line in f[1:-1]:
        line = line.split('\t')
        gene = line[0].upper()
        vals = [float(x) for x in line[1:]]
        for i in range(len(sams)):
            samDict[sams[i]][gene] = vals[i]
    return samDict

def getSigDict(dirSigs,selSigs):
    """returns a dictionary of signature:gene:value for each signature in selSigs"""
    #get all signature files
    allSigs = glob.glob(dirSigs+'/*--*')
    sigDict = {}
    for sigPath in allSigs:
        _,sig = os.path.split(sigPath)
        sig =  sig.split('--')[0]
        if sig not in selSigs:
            continue
        sigDict[sig] = {}
        sigFile = open(sigPath).read().split('\n')
        for line in sigFile:
            if not line:    continue
            line = line.split('\t')
            gene = line[0].upper()
            sigDict[sig][gene] = float(line[1])
    return sigDict

def getSelSigs(group):
    selSigs = []
    abbrevs = open(group).read().split('\n')
    for abbrev in abbrevs[:-1]:
        selSigs.append(abbrev.split('\t')[0])
    return selSigs

def loadSigDictionary(strPathToGroupFile):
    returnDict = {}
    f = open(strPathToGroupFile,'r')
    for line in f:
        abbrev, real = line.split('\t')
        returnDict[abbrev] = real[:-2]
    return returnDict

def createHeatMap(strMatrixFile,strOutFile,strVersion,strColumnZ,strRowMetric,strColMetric,invert,fixed):
    #print strMatrixFile
    #print strOutFile
    strTxtOutFile = 'routput.txt'
    strDend1OutFile = 'dend1.pdf'
    strDend2OutFile = 'dend2.pdf'
    subprocess.call(["Rscript",'heatsig.R',strMatrixFile,strOutFile,strColumnZ,strRowMetric,strColMetric,strTxtOutFile, strDend1OutFile, strDend2OutFile])
    import make_heatmap
    centerAroundZero = True
    maxVal,minVal = 1,-1
    make_heatmap.generateCanvas(strTxtOutFile,'heatmapoutput.html','Matrix Z-Score' if strColumnZ == 'matrix' else 'Value',
                                    invert,centerAroundZero,minVal,maxVal)
                                    

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--gene_count", dest="gene_count", help="table of gene counts", metavar="TAB", default="")
    parser.add_option("-s", "--sig", dest="sig", help="signature directory", metavar="LIST", default="")
    parser.add_option("-g", "--group", dest="group", help="path to group files", metavar="PATH", default="")
    parser.add_option("-r", "--row_metric", dest="row_metric", help="metric for clustering rows (samples)", default="pear_cor")
    parser.add_option("-c", "--col_metric", dest="col_metric", help="metric for clustering columns (signatures)", default="pear_cor")
    parser.add_option("-i", "--invert",default=True,dest="invert",help="heatmap columns inverted")
    parser.add_option("-f", "--fixed", dest="fixed", default="none", help="use fixed color axis")
    (options, args) = parser.parse_args()
    #check that it's in the form we expect and replace non-unix line endings
    strOutMatrixTxt = 'spearMatrix.txt'  
    selSigs = getSelSigs(options.group)
    #print selSigs[70:100]
    selSigs = selSigs[:70] + selSigs[95:]
    #selSigs = ['IL10_6h_05-15','IL10_24h_05-15','IL15_6h_05-15','IL15_24h_05-15','IFNB_2','IFNB_6','IFNB_24',
	#			'IFNG_2','IFNG_6','IFNG_24','TLP','LLP']
    #load dictionaries for sigs
    sigs = getSigDict(options.sig,selSigs)
    #load dictionary for samples
    sams = getSamDict(options.gene_count)
    genes = getSpearmanGenes(sams,sigs,'top',50)
    spearmanSams, spearmanSigs = getSpearmanDict(sams,genes),getSpearmanDict(sigs,genes)
    spearmanRank(spearmanSigs,spearmanSams,genes,strOutMatrixTxt,loadSigDictionary(options.group))
    #generate heatmap pdf
    createHeatMap(strOutMatrixTxt,'spearHeatmap.pdf','spearman','none',options.row_metric,options.col_metric,False if options.invert=='none' else True,False if options.fixed=='none' else True)
