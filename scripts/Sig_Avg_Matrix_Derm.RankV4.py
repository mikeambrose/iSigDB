#==============================================================================
# Date Modified: 2015.04.24
#==============================================================================
import os
from optparse import OptionParser
import numpy as np
import glob
import os
import gzip
import subprocess
import math
from collections import OrderedDict
import iSigDBUtilities as util
S_GENE_COUNT = '50'
S_VERSION = 'rank_delta'
S_ZSCORE = 'column'

def getSigGenes(strSigFile,nGeneCount):
    dSigToGene = {}
    with open(strSigFile) as sigToGenes:
        for line in sigToGenes:
            splitLine = line.strip("\n").split("\t")
            sig = splitLine[0]
            dSigToGene[splitLine[0]] = splitLine[1:nGeneCount+1]
    return dSigToGene

def readSigFiles(lFiles):
    dTisToFCToGene = {}
    dTisToGeneToFC = {}
    for strFile in lFiles:
        head,tail = os.path.split(strFile)
        strTis =  tail.split('--')[0]
        dTisToFCToGene[strTis]={}
        dTisToGeneToFC[strTis] = {}
        for strLine in open(strFile,'r'):
            strGene,strRatio = strLine.rstrip().split('\t')
            strGene = strGene.upper()
            if strRatio == '0':
                continue
            flogFC = float(strRatio)
            dTisToGeneToFC[strTis][strGene] = flogFC
            if dTisToFCToGene[strTis].has_key(flogFC) == False:
                dTisToFCToGene[strTis][flogFC] = []
            dTisToFCToGene[strTis][flogFC].append(strGene)
    return dTisToFCToGene,dTisToGeneToFC

def getSigToGene(dTisToFCToGene,nGeneCount):
    dSigToGene = {}

    'load DN signatures'
    '''
    for strTis in dTisToFCToGene.keys():
        dSigToGene[strTis+'_DN'] = []
        nCurGeneCount = 0
        for fFC in sorted(dTisToFCToGene[strTis].keys()):
            for strCurGene in dTisToFCToGene[strTis][fFC]:
                nCurGeneCount+=1
                dSigToGene[strTis+'_DN'].append(strCurGene)
            if nCurGeneCount >= nGeneCount:
                break
    '''

    'load UP signatures'
    for strTis in dTisToFCToGene.keys():
        dSigToGene[strTis] = []
        nCurGeneCount = 0
        for fFC in sorted(dTisToFCToGene[strTis].keys(),reverse=True):
            for strCurGene in dTisToFCToGene[strTis][fFC]:
                nCurGeneCount+=1
                dSigToGene[strTis].append(strCurGene)
            if nCurGeneCount >= nGeneCount:
                break
    return dSigToGene

def getSamToGeneToRank(dSamToCountToGene):
    dSamToGeneToRank = OrderedDict()
    dGeneToMeanRank = {}

    #process rank
    for strCurSam in dSamToCountToGene:
        nRank = 1
        currentRankVal = None
        currentRunLength = 0
        dSamToGeneToRank[strCurSam] = {}
        for fCount in sorted(dSamToCountToGene[strCurSam].keys()):
            meanRank = nRank + len(dSamToCountToGene[strCurSam][fCount])/2.0
            for strCurGene in sorted(dSamToCountToGene[strCurSam][fCount]):
                dSamToGeneToRank[strCurSam][strCurGene] = meanRank
                nRank+=1

    #process mean gene rank
    for strCurGene in sorted(dSamToGeneToRank[dSamToGeneToRank.keys()[0]].keys()):
        lRowRanks = []
        for strCurSam in sorted(dSamToGeneToRank.keys()):
            lRowRanks.append(dSamToGeneToRank[strCurSam][strCurGene])
            dGeneToMeanRank[strCurGene] = np.average(lRowRanks)

    #QC print rankings

    #fout = open('./test_rank.log','w')
    #strHeader = 'GENE'
    #for strCurSam in sorted(dSamToGeneToRank.keys()):
    #    strHeader+='\t'+strCurSam
    #fout.write(strHeader+'\tAVG\n')
    #for strCurGene in sorted(dSamToGeneToRank[dSamToGeneToRank.keys()[0]].keys()):
    #    strRow = strCurGene
    #    for strCurSam in sorted(dSamToGeneToRank.keys()):
    #        strRow+='\t%d'%dSamToGeneToRank[strCurSam][strCurGene]
    #    fout.write(strRow+'\t%f\n'%dGeneToMeanRank[strCurGene])
    #fout.close()

    return dSamToGeneToRank,dGeneToMeanRank

def procGeneCountMatrix(strGeneCount,dSigToGenes,lSigs,strOutFile,strVersion, bIsLog,sigNames):
    dColIDToColLbl = {}
    dSamToGeneToCount = OrderedDict()
    dSamToSigToLVals = OrderedDict()
    dSamToCountToGene = OrderedDict()
    lGlobalAvg = []

    #load Sample To Gene To Count
    fopen = open(strGeneCount,'r')

    for strLine in fopen:
        lcols = strLine.rstrip().split('\t')
        if len(dColIDToColLbl.keys())==0:
            for i in range(1,len(lcols)):
                dColIDToColLbl[i] = lcols[i]
                dSamToGeneToCount[lcols[i]] = {}
                dSamToSigToLVals[lcols[i]] = {}
                dSamToCountToGene[lcols[i]] = {}
                for strSig in dSigToGenes.keys():
                    dSamToSigToLVals[lcols[i]][strSig] = []
        else:
            strCurGene = lcols[0]
            strCurGene = strCurGene.upper()
            #for averaging
            """if strCurGene in numGenes:
                numGenes[strCurGene] += 1
            else:
                numGenes[strCurGene] = 1"""
            if all(float(lcols[i])==0 for i in range(1,len(lcols))):
                continue
            for i in range(1,len(lcols)):
                if strCurGene not in dSamToGeneToCount[dColIDToColLbl[i]]:
                    dSamToGeneToCount[dColIDToColLbl[i]][strCurGene] = float(lcols[i])
                    if dSamToCountToGene[dColIDToColLbl[i]].has_key(float(lcols[i])) == False:
                        dSamToCountToGene[dColIDToColLbl[i]][float(lcols[i])] = []
                    dSamToCountToGene[dColIDToColLbl[i]][float(lcols[i])].append(strCurGene)
                else:
                    dSamToGeneToCount[dColIDToColLbl[i]][strCurGene] = max(float(lcols[i]),dSamToGeneToCount[dColIDToColLbl[i]][strCurGene])
    fopen.close()

    #qc print gene
    #===================================================================================================================
    #write avg log gene count
    #===================================================================================================================
    if strVersion == 'log':
        for strSig in dSigToGenes.keys():
            for strCurSigGene in dSigToGenes[strSig]:
                for strSam in dSamToGeneToCount.keys():
                    if(dSamToGeneToCount[strSam].has_key(strCurSigGene) == False):
                        continue
                    if bIsLog:
                        dSamToSigToLVals[strSam][strSig].append(np.log10(dSamToGeneToCount[strSam][strCurSigGene]+1))
                    else:
                        dSamToSigToLVals[strSam][strSig].append(dSamToGeneToCount[strSam][strCurSigGene])
        #write outputs
        util.writeDetailedOutput(dSigToGenes,dSamToGeneToCount,strOutFile+'.full.txt',sigNames)
        util.writeRegularOutput(dSamToSigToLVals,strOutFile,sigNames)
    #===================================================================================================================
    #write delta rank metric
    #===================================================================================================================
    else:
        dSamToGeneToRank,dGeneToMeanRank = getSamToGeneToRank(dSamToCountToGene)

        #dSamToSigToLVals DELTA MEAN MATRIX
        for strSig in dSigToGenes.keys():
            for strCurSigGene in dSigToGenes[strSig]:
                for strSam in dSamToGeneToRank.keys():
                    if(dSamToGeneToRank[strSam].has_key(strCurSigGene) == False):
                        continue
                    if strVersion == 'rank_avg':
                        dSamToSigToLVals[strSam][strSig].append(dSamToGeneToRank[strSam][strCurSigGene])
                        lGlobalAvg.append(dSamToGeneToRank[strSam][strCurSigGene])
                    elif strVersion == 'rank_delta':
                        dSamToSigToLVals[strSam][strSig].append(dSamToGeneToRank[strSam][strCurSigGene] - dGeneToMeanRank[strCurSigGene])
                        lGlobalAvg.append(dSamToGeneToRank[strSam][strCurSigGene] - dGeneToMeanRank[strCurSigGene])
        util.writeDetailedOutput(dSigToGenes,dSamToGeneToRank,strOutFile+'.full.txt',sigNames)
        util.writeRegularOutput(dSamToSigToLVals,strOutFile,sigNames)

#------SPEARMAN/PEARSON CORRELATION
def corrRank(sigs,sams,geneNames,strOutFile,sigNames,version):
    """Uses the spearman/pearson coefficient to evaluate the similarity between signatures and samples
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
            coef,rankDiffs = getSpearmanVals(sigs[sig],sams[sam]) if version=="spearman"\
                                    else getPearsonVals(sigs[sig],sams[sam])
            for i in range(len(geneNames)):
                geneDetail[sig][geneNames[i]].append(rankDiffs[i])
            sigSamCoefficient[(sig,sam)] = coef
    #write detailed output
    with open(strOutFile+".full.txt","w") as detailedOut:
        for sig in sigs:
            detailedOut.write('# Signature: ' + (sigNames[sig] if sig in sigNames else sig)+'\n')
            #write sample names
            for sam in sams:
                etailedOut.write('\t'+sam)
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
    rho = 1 - 6*sum(x**2 for x in d)/float(n*(n**2-1))
    return rho,d

def getPearsonVals(v1,v2):
    """Returns the pearson correlation coefficient between v1 and v2 as well as the list
    of differentials between v1 and v2"""
    n = len(v1)
    d = [v1[i] - v2[i] for i in range(n)]
    rho = (n * sum(v1[i]*v2[i] for i in range(n)) - sum(v1)*sum(v2))/\
            math.sqrt((n*sum(x**2 for x in v1) - sum(v1)**2)*(n*sum(x**2 for x in v2)-sum(v2)**2))
    return rho,d

def getCommonGenes(candGenes,allSets):
    """Returns the subset of candGenes which occur in every member of allSets"""
    commonGenes = set()
    for gene in candGenes:
        if all(gene in s for s in allSets):
            commonGenes.add(gene)
    return commonGenes

def getSpearmanGenes(sams,sigs,compType,sigFile=None,n=50):
    """compType determines if all genes will be accessed or only the top genes
    if 'all' is passed in, computes the set of genes which are common to all samples and signatures
    if 'top' is passed in, takes the top n genes from each signature (as determined by the signature
    file) and computes the intersection of all those top genes, and then selects the subset of those
    which are present in all samples and signatures"""
    if compType=='all': #picks all genes in common
        commonGenes = set()
        candGenes = sams[sams.keys()[0]].keys()
        return list(getCommonGenes(candGenes,[sigs[sig] for sig in sigs]))
    elif compType=='top': #picks the top n
        candGenes = set()
        #generates a set of candidate genes by picking the top n genes for each sig
        with open(sigFile) as fullSigs:
            for line in fullSigs:
                line = line.split('\t')
                sig = line[0]
                if sig in sigs:
                    candGenes = candGenes.union(set(line[1:n+1]))
        allSets = [sams[sams.keys()[0]]] + [sigs[sig] for sig in sigs]
        return list(getCommonGenes(candGenes,allSets))
    elif compType=='mag': #picks all genes which are up/downregulated by a factor of n
        candGenes = set()
        for sig in sigs:
            for gene in sigs[sig]:
                if sigs[sig][gene] >= n:# or sigs[sig][gene] <= 1.0/n:
                    candGenes.add(gene)
        allSets = [sams[sams.keys()[0]]] + [sigs[sig] for sig in sigs]
        return list(getCommonGenes(candGenes,allSets))
    else:
        util.displayErrorMessage("Not a valid spearman gene selector " + str(compType))

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
    samDict = OrderedDict()
    f = open(f).read().split('\n')
    sams = f[0].replace("\n","").replace("\r","").split('\t')[1:]
    for sam in sams:
        samDict[sam] = {}
    for line in f[1:-1]:
        line = line.split('\t')
        if all(float(line[i])==0 for i in range(1,len(line))):
            continue
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
        #removing non-pms
        if 'PM' not in sig.upper():
            continue
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


#------END SPEARMAN/PEARSON CORRELATION
def writeNullModelHists(filename,sigNames,allValues,n,num_iter=100000,num_buckets=1000):
    """Writes each of the histograms to a pdf
    sigNames[i] should correspond with allValues[i]"""
    pdf = PdfPages(filename)
    for i in range(len(allValues)):
        nullmodel.getStatistics(allValues[i],n,pdf,sigNames[i],num_iter,num_buckets)
    pdf.close()

def loadGroupInfo(strPathToGroupFile):
    dGroupToLTisDesc = {}
    strGroup = strPathToGroupFile.split('.')[0]
    dGroupToLTisDesc[strGroup] = []
    for strLine in open(strPathToGroupFile,'r'):
        if strLine == '':
            continue
        lcols = strLine.rstrip().split('\t')
        strTissue = lcols[0]
        strDesc = lcols[1]
        dGroupToLTisDesc[strGroup].append((strTissue,strDesc))
    return dGroupToLTisDesc

def loadSigDictionary(strPathToGroupFile):
    returnDict = {}
    f = open(strPathToGroupFile,'r')
    for line in f:
        abbrev, real = line.split('\t')
        returnDict[abbrev] = real.replace("\n", "")
    return returnDict

def createHeatMap(strMatrixFile,strOutFile,strVersion,strColumnZ,strRowMetric,strColMetric,jobID,invert,fixed,isClient,nullFilename):
    if not isClient:
        strTxtOutFile = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_'+jobID+'/'+jobID+'.matrixForHC.txt'
        strHeatmapOutFile = strOutFile+'Rheatmap.pdf'
    else:
        strTxtOutFile = '/home/mike/workspace/PellegriniResearch/scripts/scratch/rOutput.txt'
        strHeatmapOutFile = '/home/mike/workspace/PellegriniResearch/output/RHeatmap.pdf'
    rscriptPath = "Rscript" if isClient else "/UCSC/Pathways-Auxiliary/UCLApathways-R-3.1.1/R-3.1.1/bin/Rscript"
    heatsigPath = "/home/mike/workspace/PellegriniResearch/scripts/heatsigV4.R" if isClient else '/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/heatsigV4.R'
    FNULL = open(os.devnull, 'w')
    subprocess.call([rscriptPath,heatsigPath,strMatrixFile,strHeatmapOutFile,strColumnZ,strRowMetric,strColMetric,strTxtOutFile],stdout=FNULL,stderr=FNULL)
    centerAroundZero = (strColumnZ == 'matrix') or (strVersion in ["rank_delta","pearson","spearman"])
    maxVal,minVal = None,None
    if fixed:
        if strColumnZ == 'matrix':
            minVal = -5
            maxVal = 5
        else:
            if strVersion == "rank_delta":
                maxVal = 10000
                minVal = -10000
            elif strVersion == "rank_avg":
                maxVal = 30000
                minVal = 0
            elif strVersion == "log": #TODO: fix version names
                maxVal = 10
                minVal = 0
            elif strVersion in ['pearson','spearman']:
                maxVal = 1
                minVal = -1
    import make_heatmapV4 as make_heatmap
    out = '/home/mike/workspace/PellegriniResearch/output/HighChartsHeatmap.html' if isClient else None
    includeDetailed = strVersion not in ['pearson','spearman']
    if nullFilename:
        nullFilename="http://pathways-pellegrini.mcdb.ucla.edu//submit/img/" + os.path.basename(nullFilename)
    make_heatmap.generateCanvas(strTxtOutFile, out,'Matrix Z-Score' if strColumnZ == 'matrix' else 'Value',invert,centerAroundZero,minVal,maxVal,strOutFile,includeDetailed,nullFilename)

#----------------------------------------------------------------------------
# main function call
#----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--gene_count", dest="gene_count", help="table of gene counts", metavar="TAB", default="")
    parser.add_option("-s", "--sig", dest="sig", help="path to pm ratio list", metavar="LIST", default="")
    parser.add_option("-x","--sigfile",dest="sigfile",help="path to precomputed signatures")
    parser.add_option("-g", "--group", dest="group", help="path to group files", metavar="PATH", default="")
    parser.add_option("-n", "--ngene_count", dest="ngene_count", help="number of genes to use", metavar="NUM", default=S_GENE_COUNT)
    parser.add_option("-v", "--version", dest="version", help="metric version (rank_avg,rank_delta,log)", metavar="VER", default=S_VERSION)
    parser.add_option("-z", "--zscore", dest="zscore", help="default no column z-score", metavar="ZSCORE", default=S_ZSCORE)
    parser.add_option("-j", "--job_id", dest="job_id", help="job id, for output path", default="-1")
    parser.add_option("-a", "--absolute", dest="bIsLog", action="store_false", default=True, help="no log transformation")
    parser.add_option("-l", "--log_transform", dest="bIsLog", action="store_true", help="log transform")
    parser.add_option("-r", "--row_metric", dest="row_metric", help="metric for clustering rows (samples)", default="pear_cor")
    parser.add_option("-c", "--col_metric", dest="col_metric", help="metric for clustering columns (signatures)", default="pear_cor")
    parser.add_option("-i","--invert",default=True,dest="invert",help="heatmap columns/rows swtiched")
    parser.add_option("-f", "--fixed", dest="fixed", default="none", help="use fixed color axis")
    parser.add_option("-o", "--gene_option", dest="gene_option", help="spearman/pearson gene computation option")
    parser.add_option("-u", "--null", dest="null", help="compute null model")
    parser.add_option("--client", dest="client", action="store_true", default=True, help="running client-side")
    parser.add_option("--server", dest="client", action="store_false", help="running server-side")
    (options, args) = parser.parse_args()

    #matplotlib imports
    if not options.client:
        os.environ["MPLCONFIGDIR"] = "/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space"
    import matplotlib
    matplotlib.use('Agg')
    import nullmodel
    from matplotlib.backends.backend_pdf import PdfPages

    util.reformatFile(options.gene_count)
    util.checkForErrors(options.gene_count)

    if options.client:
        strOutMatrixTxt = '/home/mike/workspace/PellegriniResearch/scripts/scratch/output.txt'
    else:
        strOutMatrixTxt = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_'+options.job_id+'/'+options.job_id+'.matrix.txt'
    #null model
    nullFilename = '/home/mike/workspace/PellegriniResearch/output/nulldist.pdf' if options.client\
                        else '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/nulldist_' + options.job_id + '.pdf'
    if options.version == 'log' and not options.bIsLog and options.null != 'none':
        sams = getSamDict(options.gene_count)
        names = []
        allVals = []
        for sam in sorted(sams.keys()):
            names.append(sam)
            allVals.append([sams[sam][gene] for gene in sams[sam]])
        writeNullModelHists(nullFilename,names,allVals,int(options.ngene_count))
    elif options.version == 'log' and options.bIsLog and options.null != 'none':
        sams = getSamDict(options.gene_count)
        names = []
        allVals = []
        for sam in sams:
            names.append(sam)
            allVals.append([math.log(sams[sam][gene]+1,10) for gene in sams[sam]])
        writeNullModelHists(nullFilename,names,allVals,int(options.ngene_count))
    elif options.version == 'rank_avg' and options.null != 'none':
        names = ['All Samples']
        d = getSamDict(options.gene_count)
        num_genes = len(d[d.keys()[0]])
        allVals = [list(range(1,num_genes))]
        writeNullModelHists(nullFilename,names,allVals,int(options.ngene_count))
    else:
        nullFilename = None

    if options.version not in ['spearman','pearson']:
        #load categories
        dGroupToLTisDesc = loadGroupInfo(options.group)

        #list signature
        lTis = []
        for strGroup in sorted(dGroupToLTisDesc.keys()):
            for strTis,strDes in dGroupToLTisDesc[strGroup]:
                lTis.append(strTis)

        #load list of genes
        dGroupSigToGene = getSigGenes(options.sigfile,int(options.ngene_count))

        #process gene count matrix
        procGeneCountMatrix(options.gene_count,dGroupSigToGene,lTis,strOutMatrixTxt,options.version,options.bIsLog,loadSigDictionary(options.group))
    else: #spearman or pearson
        selSigs = getSelSigs(options.group)
        sigs = getSigDict(options.sig,selSigs)
        sams = getSamDict(options.gene_count)
        #filter sigs with few genes
        i = 0
        while i < len(sigs):
            sig = sigs.keys()[i]
            if len(sigs[sig]) < 15000:
                del sigs[sig]
            else:
                i += 1
        genes = getSpearmanGenes(sams,sigs,options.gene_option,options.sigfile,int(options.ngene_count))
        if len(genes) == 0:
            util.displayErrorMessage("There are no genes which are in common between your samples and all signatures. Make sure the gene name in your sample is in the first column")
        spearmanSams, spearmanSigs = getSpearmanDict(sams,genes), getSpearmanDict(sigs,genes)
        corrRank(spearmanSigs,spearmanSams,genes,strOutMatrixTxt,loadSigDictionary(options.group),options.version)
        print "Based on " + str(len(genes)) + " genes"
    #generate heatmap pdf
    RHeatmapOut = '/home/mike/workspace/PellegriniResearch/scripts/scratch/' if options.client else '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/goTeles_tissueDeconvolution_'+options.job_id
    createHeatMap(strOutMatrixTxt,RHeatmapOut,options.version,options.zscore,options.row_metric,options.col_metric,options.job_id,False if options.invert=='none' else True, False if options.fixed == 'none' else True,options.client,nullFilename)
