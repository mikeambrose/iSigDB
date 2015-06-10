from __future__ import division
def column(lst,col):
    return [row[col] for row in lst]
def median(lst):
    cpy = lst[::]
    cpy.sort()
    if len(cpy) % 2 == 0:
        halfway = len(cpy)//2
        return (cpy[halfway] + cpy[halfway-1])/2
    else:
        return cpy[len(cpy)//2]
floatAll = lambda x: [float(p) for p in x]

def reformatData(inFile,outFiles,impColumns):
    assert len(outFiles) == len(impColumns)
    f = open(inFile).read().split('\n')
    while not f[-1]:    del f[-1]
    data = [x.split(',') for x in f]
    geneNames = column(data,1)[1:]
    for i in range(len(impColumns)):
        vals = column(data,impColumns[i])
        name = vals[0]
        del vals[0]
        vals = floatAll(vals)
        totalMedian = median(vals)
        print totalMedian,i
        assert len(geneNames) == len(vals)
        geneDict = dict()
        for j in range(len(geneNames)):
            if geneNames[j] in geneDict:
                geneDict[geneNames[j]].append(vals[j])
            else:
                geneDict[geneNames[j]] = [vals[j]]
        for gene in geneDict:
            geneDict[gene] = median(geneDict[gene])/totalMedian
        f_out = open(outFiles[i],'w')
        geneTuples = geneDict.items()
        geneTuples.sort(key = lambda x:-x[1])
        for t in geneTuples:
            f_out.write(t[0] + "\t" + str(t[1]) + "\n")
        f_out.close()

def reformatData_compMedian(inFile,outFiles,baseMedianCol,impColumns):
    assert len(outFiles) == len(impColumns)
    f = open(inFile).read().split('\n')
    while not f[-1]:    del f[-1]
    data = [x.split(',') for x in f]
    geneNames = column(data,1)[1:]
    basevals = floatAll(column(data,baseMedianCol)[1:])
    baseDict = dict()
    baseMedian = median(basevals)
    for j in range(len(geneNames)):
        if geneNames[j] in baseDict:
            baseDict[geneNames[j]].append(basevals[j])
        else:
            baseDict[geneNames[j]] = [basevals[j]]
    for gene in baseDict:
        baseDict[gene] = median(baseDict[gene])/baseMedian 

    for i in range(len(impColumns)):
        vals = column(data,impColumns[i])
        name = vals[0]
        del vals[0]
        vals = floatAll(vals)
        totalMedian = median(vals)
        assert len(geneNames) == len(vals)
        geneDict = dict()
        for j in range(len(geneNames)):
            if geneNames[j] in geneDict:
                geneDict[geneNames[j]].append(vals[j])
            else:
                geneDict[geneNames[j]] = [vals[j]]
        for gene in geneDict:
            if baseDict[gene] == 0 and median(geneDict[gene]) == 0:
                geneDict[gene] = 0
            elif baseDict[gene] == 0 and median(geneDict[gene]) != 0:
                print "weird case"
                geneDict[gene] = 1
            else:
                geneDict[gene] = median(geneDict[gene])/totalMedian / baseDict[gene]
        f_out = open(outFiles[i],'w')
        geneTuples = geneDict.items()
        geneTuples.sort(key = lambda x:-x[1])
        for t in geneTuples:
            f_out.write(t[0] + "\t" + str(t[1]) + "\n")
        f_out.close()

reformatData('Hum_avg_mLEP_FoldChange_matteo.csv',['mLEP-' + x + '-human--pms_sorted.txt' for x in ['t0','t2','t6','t24']],[5,8,11,14])
reformatData('Mouse_avg_mLEP_FoldChange_matteo.csv',['mLEP-' + x + '-mouse--pms_sorted.txt' for x in ['t0','t2','t6','t24']],[6,10,14,18])
#reformatData_compMedian('Hum_avg_mLEP_FoldChange_matteo.csv',['mLEP-' + x + '-human-medt0--pms_sorted.txt' for x in ['t0','t2','t6','t24']],5,[5,8,11,14])
#reformatData_compMedian('Mouse_avg_mLEP_FoldChange_matteo.csv',['mLEP-' + x + '-mouse-medt0--pms_sorted.txt' for x in ['t0','t2','t6','t24']],6,[6,10,14,18])
