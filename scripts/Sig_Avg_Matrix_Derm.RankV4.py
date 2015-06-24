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
    procGeneCountMatrix(options.gene_count,dGroupSigToGene,lTis,strOutMatrixTxt,options.version,options.bIsLog,util.loadAbbrevs(options.group))

    #generate heatmap pdf
    RHeatmapOut ='/home/mike/workspace/PellegriniResearch/scripts/scratch/Rheatmap.pdf' if options.client\
        else '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/goTeles_tissueDeconvolution_{0}/{0}Rheatmap.pdf'.format(options.job_id)
    util.createHeatMap(strOutMatrixTxt,RHeatmapOut,options.version,options.zscore,options.row_metric,options.col_metric,options.job_id,False if options.invert=='none' else True, False if options.fixed == 'none' else True,options.client,nullFilename)
