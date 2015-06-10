#==============================================================================
# Date Modified: 2015.04.24
#==============================================================================
from optparse import OptionParser
import numpy as np
import glob
import os
import gzip
import subprocess
S_GENE_COUNT = '50'
S_VERSION = 'rank_delta'
S_ZSCORE = 'column'

def replaceLineEndings(f):
    s = open(f).read()
    unix = '\n' in s
    dos = '\r' in s
    if unix and dos:
        s = s.replace('\r\n','\n')
        if '\r' in s: #we didn't successfully rid ourselves
            errorMessage("File has non-standard line endings")
    elif dos:
        s = s.replace('\r','\n')
    if dos:
        with open(f,'w') as fout:
            fout.write(s)

def errorMessage(s):
    print "<!doctype html> <html> <body>" + s + "</body> </html>"
    exit()

def checkForErrors(f):
    """Checks that f is vaguely in the format we expect
    Checks that it appears to be tab-delimited with a consistent number of columns
    and that all entries are decimals.
    If not, prints out html corresponding to an error message and exits execution"""
    f = open(f).read().split('\n')
    if len(f) <= 1:
        errorMessage("Formatting error: Only one line detected")
    numTabs = f[0].count('\t')
    if numTabs == 0:
        errorMessage("Not a tab-separated file")
    for i in range(1,len(f)-1):
        line = f[i]
        if line.count('\t') != numTabs:
            errorMessage("Inconsistent number of columns around line " + str(i+1))
        for x in line.split('\t')[1:]:
            try:
                float(x)
            except:
                errorMessage("Non-decimal alue around line " + str(i+1))

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
    dSamToGeneToRank = {}
    dGeneToMeanRank = {}

    #process rank
    for strCurSam in dSamToCountToGene:
        nRank = 1
        dSamToGeneToRank[strCurSam] = {}
        for fCount in sorted(dSamToCountToGene[strCurSam].keys()):
            for strCurGene in sorted(dSamToCountToGene[strCurSam][fCount]):
                dSamToGeneToRank[strCurSam][strCurGene] = nRank
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
    dSamToGeneToCount = {}
    dSamToSigToLVals = {}
    dSamToCountToGene = {}
    lGlobalAvg = []

    #load Sample To Gene To Count
    fopen = None
    if strGeneCount.endswith('.gz'):
        fopen = gzip.open(strGeneCount,'rb')
    else:
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
            for i in range(1,len(lcols)):
                dSamToGeneToCount[dColIDToColLbl[i]][strCurGene] = float(lcols[i])
                if dSamToCountToGene[dColIDToColLbl[i]].has_key(float(lcols[i])) == False:
                    dSamToCountToGene[dColIDToColLbl[i]][float(lcols[i])] = []
                dSamToCountToGene[dColIDToColLbl[i]][float(lcols[i])].append(strCurGene)
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

        fout2 = open(strOutFile+'.full.txt', 'w')
        for strSig in sorted(dSigToGenes.keys()):
            fout2.write('# Signature: %s\n' % sigNames[strSig] if strSig in sigNames else strSig)
            for strSam in sorted(dSamToGeneToCount.keys()):
                fout2.write('\t%s' % strSam)
            fout2.write('\n')
            for strCurSigGene in dSigToGenes[strSig]:
                fout2.write('%s' % strCurSigGene)
                for strSam in sorted(dSamToGeneToCount.keys()):
                    if strCurSigGene in dSamToGeneToCount[strSam]:
                        if not bIsLog:
                            fout2.write('\t%f' % (dSamToGeneToCount[strSam][strCurSigGene]))
                        else:
                            fout2.write('\t%f' % (np.log10(dSamToGeneToCount[strSam][strCurSigGene]+1)))
                    else:
                        fout2.write('\t%N/A')
                fout2.write('\n')
            fout2.write('\n')

        #print header
        fout = open(strOutFile,'w')
        strHeader = 'SAMPLE'
        for i in range(len(lSigs)):
            strHeader+='\t'+(sigNames[lSigs[i]] if lSigs[i] in sigNames else lSigs[i])
        fout.write(strHeader+'\n')
        #print sample Row Averages
        for strSam in sorted(dSamToGeneToCount.keys()):
            strRow = strSam
            for j in range(len(lSigs)):
                strCurSig = lSigs[j]
                strRow += '\t%f'%np.average(dSamToSigToLVals[strSam][strCurSig])
            fout.write(strRow+'\n')
        fout.close()

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
                    if strVersion.startswith('rank_avg'):
                        dSamToSigToLVals[strSam][strSig].append(dSamToGeneToRank[strSam][strCurSigGene])
                        lGlobalAvg.append(dSamToGeneToRank[strSam][strCurSigGene])
                    elif strVersion.startswith('rank_delta'):
                        dSamToSigToLVals[strSam][strSig].append(dSamToGeneToRank[strSam][strCurSigGene] - dGeneToMeanRank[strCurSigGene])
                        lGlobalAvg.append(dSamToGeneToRank[strSam][strCurSigGene] - dGeneToMeanRank[strCurSigGene])
        #lLoadGlobalAvg
        '''
        for strSam in sorted(dSamToGeneToCount.keys()):
            strRow = strSam
            for j in range(len(lSigs)):
                strCurSig = lSigs[j]
                lGlobalAvg.append(np.average(dSamToSigToLVals[strSam][strCurSig]))
        '''

        # TODO: extra output here
        fout2 = open(strOutFile+'.full.txt', 'w')
        for strSig in sorted(dSigToGenes.keys()):
            fout2.write('# Signature: %s\n' % (sigNames[strSig] if strSig in sigNames else strSig))
            for strSam in sorted(dSamToGeneToRank.keys()):
                fout2.write('\t%s' % strSam)
            fout2.write('\n')
            for strCurSigGene in dSigToGenes[strSig]:
                fout2.write('%s' % strCurSigGene)
                for strSam in sorted(dSamToGeneToRank.keys()):
                    if(dSamToGeneToRank[strSam].has_key(strCurSigGene) == False):
                        fout2.write('\tN/A')
                    else:
                        fout2.write('\t%f' % dSamToGeneToRank[strSam][strCurSigGene])
                fout2.write('\n')
            fout2.write('\n')

        #find signature gene
        fout = open(strOutFile,'w')
        strHeader = 'SAMPLE'
        for i in range(len(lSigs)):
            strHeader+='\t'+(sigNames[lSigs[i]] if lSigs[i] in sigNames else lSigs[i])
        fout.write(strHeader+'\n')
        #print sample Row Averages
        for strSam in sorted(dSamToGeneToCount.keys()):
            strRow = strSam
            for j in range(len(lSigs)):
                strCurSig = lSigs[j]
                if strVersion.endswith('_z'):
                    fMetric = (np.average(dSamToSigToLVals[strSam][strCurSig]) - np.average(lGlobalAvg))/np.std(lGlobalAvg)
                else:
                    fMetric = np.average(dSamToSigToLVals[strSam][strCurSig])
                strRow += '\t%f'%fMetric
            fout.write(strRow+'\n')
        fout.close()

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

def createHeatMap(strMatrixFile,strOutFile,strVersion,strColumnZ,strRowMetric,strColMetric,jobID,invert,fixed,isClient):
    #print strMatrixFile
    #print strOutFile
    if not isClient:
        strTxtOutFile = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_'+jobID+'/'+jobID+'.matrixForHC.txt'
        strHeatmapOutFile = strOutFile+'Rheatmap.pdf'
    else:
        strTxtOutFile = '/home/mike/workspace/PellegriniResearch/scripts/scratch/rOutput.txt'
        strHeatmapOutFile = '/home/mike/workspace/PellegriniResearch/output/RHeatmap.pdf'
    rscriptPath = "Rscript" if isClient else "/UCSC/Pathways-Auxiliary/UCLApathways-R-3.1.1/R-3.1.1/bin/Rscript"
    heatsigPath = "/home/mike/workspace/PellegriniResearch/scripts/heatsigV4.R" if isClient else '/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/heatsigV4.R'
    subprocess.call([rscriptPath,heatsigPath,strMatrixFile,strHeatmapOutFile,strColumnZ,strRowMetric,strColMetric,strTxtOutFile])
    centerAroundZero = (strColumnZ == 'matrix') or (strVersion == "rank_delta")
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
            elif strVersion == "log":
                maxVal = 10
                minVal = 0
            else:
                maxVal = 30
                minVal = 0
    import make_heatmapV4 as make_heatmap
    out = '/home/mike/workspace/PellegriniResearch/output/HighChartsHeatmap.html' if isClient else None
    make_heatmap.generateCanvas(strTxtOutFile, out,'Matrix Z-Score' if strColumnZ == 'matrix' else 'Value',invert,centerAroundZero,minVal,maxVal,strOutFile)

#----------------------------------------------------------------------------
# main function call
#----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--gene_count", dest="gene_count", help="table of gene counts", metavar="TAB", default="")
    parser.add_option("-s", "--sig", dest="sig", help="path to pm ratio list", metavar="LIST", default="")
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
    parser.add_option("--client", dest="client", action="store_true", help="running client-side")
    parser.add_option("--server", dest="client", action="store_false", default=True, help="running server-side")
    (options, args) = parser.parse_args()

    replaceLineEndings(options.gene_count)
    checkForErrors(options.gene_count)
    #load categories
    dGroupToLTisDesc = loadGroupInfo(options.group)

    #list signatures
    lTis = []
    for strGroup in sorted(dGroupToLTisDesc.keys()):
        for strTis,strDes in dGroupToLTisDesc[strGroup]:
            lTis.append(strTis)

    #load list of genes
    dGroupSigToGene = getSigGenes(options.sig,int(options.ngene_count))
    if options.client:
        strOutMatrixTxt = '/home/mike/workspace/PellegriniResearch/scripts/scratch/output.txt'
    else:
        strOutMatrixTxt = '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_'+options.job_id+'/'+options.job_id+'.matrix.txt'

    #process gene count matrix
    procGeneCountMatrix(options.gene_count,dGroupSigToGene,lTis,strOutMatrixTxt,options.version,options.bIsLog,loadSigDictionary(options.group))

    #generate heatmap pdf
    RHeatmapOut = '/home/mike/workspace/PellegriniResearch/scripts/scratch/' if options.client else '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/goTeles_tissueDeconvolution_'+options.job_id
    createHeatMap(strOutMatrixTxt,RHeatmapOut,options.version,options.zscore,options.row_metric,options.col_metric,options.job_id,False if options.invert=='none' else True, False if options.fixed == 'none' else True,options.client)
