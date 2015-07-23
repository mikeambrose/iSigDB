def readEnsemblDict(filename):
    ensGene = {}
    allGenes = set()
    with open(filename) as f:
        f.readline()
        for line in f:
            ens,gene = line.upper().replace('\n','').split('\t') 
            if not gene:    continue
            allGenes.add(gene)
            if ens in ensGene:
                ensGene[ens].append(gene)
            else:
                ensGene[ens] = [gene]
    return ensGene,allGenes

def getMouseHumanDicts(human,mouse,both):
    """Finds the set of all mappings which are one-to-one, excluding genes which map
        only to themselves
    'human','mouse','both' are their respective Ensembl filenames
        (constants contain the expected filenames)"""
    human,humanGenes = readEnsemblDict(human)
    mouse,mouseGenes = readEnsemblDict(mouse)
    humanToMouse = {}
    mouseToHuman = {}
    with open(both) as f:
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

    #impose one-to-one translation
    for d in [humanToMouse,mouseToHuman]:
        for gene in d.keys():
            if len(d[gene]) == 0 or len(d[gene]) > 1:
                del d[gene]
            else:
                d[gene] = list(set(d[gene]))
    #remove gene->gene mappings
    for d in [humanToMouse,mouseToHuman]:
        for gene in d.keys():
            if d[gene][0] == gene:
                del d[gene]
            else:
                d[gene] = d[gene][0]
    #if a gene is defined as both a mouse and human gene, remove it
    for gene in humanToMouse.keys():
        if gene in mouseGenes:
            del humanToMouse[gene]
    for gene in mouseToHuman.keys():
        if gene in humanGenes:
            del mouseToHuman[gene]
    return mouseToHuman,humanToMouse



human = 'BioMart_EnsemblGenes80_HomoSapiensGenesGRCh38.p2_EnsemblGeneID_HGNCsymbol.txt'
mouse = 'BioMart_EnsemblGenes80_MusMusculusGenesGRCm38.p3_EnsemblGeneID_MGIsymbol.txt'
both = 'BioMart_EnsemblGenes80_MusMusculusGenesGRCm38.p3_EnsemblGeneID_HumanEnsemblGeneID.txt'
