"""
Various utilities called by different functions used for iSigDB
Made by Mike Ambrose, mikeambrose@berkeley.edu
"""
from collections import OrderedDict
import subprocess
import make_heatmapV4 as make_heatmap
import os
import random
import datetime
import re

def reformatFile(f):
    """Automatically changes f to be in the correct format
    Removes carriage returns (either with blank spaces if there are \r\n line endings
        or with newlines if there are just \r line endings)
    Replaces commas with semicolons
    Removes empty lines
    """
    s = open(f).read()
    s,ch1 = replaceNewlines(s)
    s,ch2 = replaceEmptyLines(s)
    s,ch3 = s.replace(',',';'),',' in s
    if any([ch1,ch2,ch3]):
        with open(f,'w') as fout:
            fout.write(s)

def replaceNewlines(s):
    """Replaces carriage returns in s
    If there are \r\n line endings, replaces with blank space - otherwise, replaces with \n
    Returns the modified string and whether or not the string was changed"""
    if '\r' not in s:
        return s, False
    elif '\r\n' in s:
        return s.replace('\r',''),True
    else:
        return s.replace('\r','\n'),True

def replaceEmptyLines(s):
    """Removes any empty lines in s
    Returns the modified string and whether or not the string was changed"""
    if '\n\n' not in s:
        return s,False
    while '\n\n' in s:
        s = s.replace('\n\n','\n')
    return s,True

def displayErrorMessage(message,htmlHeader=False):
    """Prints a brief html snippet which corresponds to an error message with context 'message'
    if htmlHeader is true, also prints the content-type header
    Exits after printing"""
    if htmlHeader:
        print "Content-type: text/html\n\n"
    print "<!DOCTYPE html>\n<title>iSigDB Error</title>\n<html>\n<body>"+message+"\n</body>\n</html>"""
    exit()

def checkForErrors(f):
    """Looks for errors which would cause us to terminate
    Specifically, checks for
        tab-delimiters with a consistent number of columns
        unique sample names
        Decimal values
    and calls displayErrorMessage if any of these are not met"""
    f = open(f).read().split('\n')
    if len(f) <= 1:
        displayErrorMessage("Formatting error: Only one line detected")
    names = f[0].split('\t')
    if len(set(names)) != len(names):
        displayErrorMessage("Some of your samples have the same name (maybe duplicates?)")
    numTabs = f[0].count('\t')
    if numTabs == 0:
        displayErrorMessage("Not a tab-separated file")
    for i in range(1,len(f)-1):
        line = f[i]
        if line.count('\t') != numTabs:
            displayErrorMessage("Inconsistent number of columns around line " + str(i+1))
        for x in line.split('\t')[1:]:
            try:
                float(x)
            except:
                displayErrorMessage("Non-decimal value around line " + str(i+1))

def writeDetailedOutput(sigGenes,samVals,outFile,fullNames,deltaSigs=None):
    """Writes detailed output to outFile
    sigGenes is a dictionary of signature : gene
    samVals is a dictionary of sample : gene : value
    fullNames is a dictionary of abbreviation to full name
    deltaSigs is a dictionary of signature name : average signature value
        used for detailed output with a delta"""
    with open(outFile,'w') as out:
        for sig in sigGenes:
            fullName = fullNames[sig] if sig in fullNames else sig
            out.write("# Signature: {0}\n".format(fullName))
            for sam in samVals:
                out.write('\t'+sam)
            out.write('\n')
            for gene in sigGenes[sig]:
                out.write(gene)
                for sam in samVals:
                    if gene not in samVals[sam]:
                        out.write('\tN/A')
                    else:
                        if deltaSigs and sig in deltaSigs:
                            out.write('\t'+str(samVals[sam][gene]-deltaSigs[sig]))
                        else:
                            out.write('\t'+str(samVals[sam][gene]))
                out.write('\n')
            out.write('\n')

average = lambda lst: sum(lst) / float(len(lst))

def invertDict(d):
    """Changes d from a:b:val to b:a:val"""
    newD = {}
    for a in d:
        for b in d[a]:
            if b not in newD:
                newD[b] = {}
            newD[b][a] = d[a][b]
    return newD    

def writeRegularOutput(samSigVals,outFile,fullNames={}):
    """Writes regular output to outFile
    samSigVals is a dictionary of sample : signature : value
        the value can be a list, in which case it is averaged, or a single value
        if value is an empty list, an error is thrown
    fullNames is a dictionary of abbreviation to full name"""
    sigNames = samSigVals[samSigVals.keys()[0]].keys()
    with open(outFile,'w') as out:
        out.write('SAMPLE')
        for sig in sigNames:
            fullName = fullNames[sig] if sig in fullNames else sig
            out.write('\t'+fullName)
        out.write('\n')
        for sam in samSigVals:
            out.write(sam)
            for sig in samSigVals[sam]:
                if type(samSigVals[sam][sig]) == type([]):
                    if not samSigVals[sam][sig]:
                        displayErrorMessage("No genes which intersect with signature "+sig)
                    out.write('\t'+str(average(samSigVals[sam][sig])))
                else:
                    out.write('\t'+str(samSigVals[sam][sig]))
            out.write('\n')

def createHeatmap(matrixFile,rPdfOutFile,version,zTransform,rowMetric,colMetric,jobID,invert,\
                    isClient,nullFilename,inpHistFilename=None,mn=None,mx=None):
    """Calls the R script to cluster and create the heatmap
        matrixFile is the location of the output
        rPdfOutFile is where the R heatmap will be output
        zTransform ("matrix","none") determines if the output is replaced with its z-scores
        rowMetric, colMetric ("euclidean","pearson","none") determine how R clusters the output
    Then control is passed to make_heatmap, which generates the HighCharts heatmap
        invert controls whether or not the output is inverted
        fixed controls whether or not the axes are fixed
        isClient is a debug flag, always passed as False on the server
        nullFilename is the location of the null distribution pdf
    """
    if not isClient:
        rTxtOutFile = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_'+jobID+'/'+jobID+'.matrixForHC.txt'
    else:
        rTxtOutFile = '/home/mike/workspace/PellegriniResearch/scripts/scratch/rOutput.txt'
    rscriptPath = "Rscript" if isClient\
                else "/UCSC/Pathways-Auxiliary/UCLApathways-R-3.1.1/R-3.1.1/bin/Rscript"
    heatsigPath = "/home/mike/workspace/PellegriniResearch/scripts/heatsigV4.R" if isClient\
                else '/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/heatsigV4.R'
    #call R and stuff all output
    if isClient:
        subprocess.call([rscriptPath,heatsigPath,matrixFile,rPdfOutFile,zTransform,\
                    rowMetric,colMetric,rTxtOutFile])
    else:
        FNULL = open(os.devnull, 'w')
        subprocess.call([rscriptPath,heatsigPath,matrixFile,rPdfOutFile,zTransform,\
                        rowMetric,colMetric,rTxtOutFile],stdout=FNULL,stderr=FNULL)
    #only center around zero  for certain input types
    centerAroundZero = (zTransform=="matrix") or ("delta" in version)
    #if fixed is selected, choose the fixed values
    maxVal,minVal = None,None
    if mn != None:
        minVal = mn
    if mx != None:
        maxVal = mx
    #debug, location of html output from make_heatmap
    out = '/home/mike/workspace/PellegriniResearch/output/HighChartsHeatmap.html' if isClient\
            else None
    #whether or not to include the detailed output
    includeDetailed = version not in ['pearson','spearman']
    #if we have a null filename, replace it with the html-accessable one
    if nullFilename:
        nullFilename="http://pathways-pellegrini.mcdb.ucla.edu//submit/img/" +\
                            os.path.basename(nullFilename)
    if inpHistFilename:
        inpHistFilename="http://pathways-pellegrini.mcdb.ucla.edu//submit/img/" +\
                            os.path.basename(inpHistFilename)
    #pass control to make_heatmap
    make_heatmap.generateCanvas(rTxtOutFile, out,'Matrix Z-Score' if zTransform == 'matrix' else 'Value',invert,centerAroundZero,minVal,maxVal,rPdfOutFile,includeDetailed,nullFilename,inpHistFilename)

def readMatrix(f,filterAllZero=True,ordered=False):
    """accepts file of the form:
    A           0.1     0.2     ...
    B           1.2     -0.4    ...
    ...
    and returns a dictionary of the form
    {COL1 : {A : 0.1, B : 1.2}, COL2 : {A : 0.2, B : -0.4}}
    in other words, returns a dictionary of sample:gene:value
    if it encounters a N/A value, it will skip the line
    if filterAllZero is false, removes any line which is all zeroes
    if ordered is True, uses an OrderedDict() to maintain the order
    otherwise will use a regular dictionary (better for performance, especially with many samples)"""
    samDict = OrderedDict() if ordered else dict()
    f = open(f).read().split('\n')
    #names of the samples
    sams = f[0].replace("\n","").replace("\r","").split('\t')[1:]
    for sam in sams:
        samDict[sam] = {}
    for line in f[1:]:
        if not line:
            continue
        line = line.split('\t')
        if any(val == 'N/A' for val in line) or \
           (all(float(val)==0 for val in line[1:]) and filterAllZero):
            continue
        #process when there are multiple genes in one line with // operator
        genes = re.split(" ?//+ ?",line[0].upper())
        genes = [gene.strip() for gene in genes]
        vals = [float(x) for x in line[1:]]
        for gene in genes:
            for i in range(len(sams)):
                samDict[sams[i]][gene] = vals[i]
    return samDict

def getSigDict(dirSigs,selSigs):
    """dirSigs is the directory with all signature files
    selSigs is the signatures selected by the user
    returns a dictionary of signature:gene:value for each signature in selSigs"""
    #get all signature files
    allSigs = glob.glob(dirSigs+'/*--*')
    sigDict = {}
    for sigPath in allSigs:
        _,sig = os.path.split(sigPath)
        #removing non-pms
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

def loadAbbrevs(abbrevs):
    returnDict = {}
    with open(abbrevs) as a:
        for line in a:
            abbrev, real = line.replace('\n','').replace('\r','').split('\t')
            returnDict[abbrev] = real
    return returnDict

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

zeroExtend = lambda s,k: s if len(s)==k else zeroExtend("0"+s,k) #pads s with 0s until length k
def getJobID():
    """Generates the job id
    {year}{month}{day}{hour}{minute}{second}xxxx
    where xxxx is a random four-digit number"""
    randomSeed = zeroExtend(str(random.randint(0,9999)),4)
    currentTime = datetime.datetime.now()
    return "{0:04d}{1:02d}{2:02d}{3:02d}{4:02d}{5:02d}".format(currentTime.year,currentTime.month,currentTime.day,\
                            currentTime.hour,currentTime.minute,currentTime.second)+randomSeed

def copyFile(f,loc):
    """Copies the opened file f to loc"""
    with open(loc,'wb') as out:
        out.write(f.read())

