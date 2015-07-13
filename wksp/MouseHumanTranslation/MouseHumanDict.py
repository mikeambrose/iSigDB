def readEnsemblDict(filename):
    ensGene = {}
    with open(filename) as f:
        f.readline()
        for line in f:
            ens,gene = line.upper().replace('\n','').split('\t') 
            if not gene:    continue
            if ens in ensGene:
                ensGene[ens].append(gene)
            else:
                ensGene[ens] = [gene]
    return ensGene

human = readEnsemblDict('BioMart_EnsemblGenes80_HomoSapiensGenesGRCh38.p2_EnsemblGeneID_HGNCsymbol.txt')
mouse = readEnsemblDict('BioMart_EnsemblGenes80_MusMusculusGenesGRCm38.p3_EnsemblGeneID_MGIsymbol.txt')
humanToMouse = {}
mouseToHuman = {}
with open('BioMart_EnsemblGenes80_MusMusculusGenesGRCm38.p3_EnsemblGeneID_HumanEnsemblGeneID.txt')\
    as f:
    f.readline()
    for line in f:
        mEns,hEns = line.upper().replace('\n','').split('\t')
        if not hEns or hEns not in human or mEns not in mouse:  continue
        for mGene in mouse[mEns]:
            if mGene not in mouseToHuman:
                mouseToHuman[mGene] = human[hEns]
            else:
                mouseToHuman[mGene].extend(human[hEns])
        for hGene in human[hEns]:
            if hGene not in humanToMouse:
                humanToMouse[hGene] = mouse[mEns]
            else:
                humanToMouse[hGene].extend(mouse[mEns])

for d in [humanToMouse,mouseToHuman]:
    for gene in d.keys():
        while gene in d[gene]:
            d[gene].remove(gene)
        if len(d[gene]) == 0:
            del d[gene]
        else:
            d[gene] = list(set(d[gene]))
