baseConstants = {
    "R_TXT": '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_{jobID}/{jobID}.matrixForHC.txt',
    "RSCRIPT": '/UCSC/Pathways-Auxiliary/UCLApathways-R-3.1.1/R-3.1.1/bin/Rscript',
    "HEATSIG": '/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/heatsigV4.R',
    "R_DOWNLOAD": '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/data/{jobID}.data.txt',
    "HEATMAP_TEMPLATE": '/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/heatmapTemplate.html',
    "COMP_OUTPUT": '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_{jobID}/{jobID}.matrix.txt',
    "NULL_PDF": '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/nulldist_{jobID}.pdf',
    "INPUT_DIST_PDF": '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/inputdist_{jobID}.pdf',
    "R_HEATMAP": '/UCSC/Apache-2.2.11/htdocs-UCLApathways-pellegrini/submit/img/{jobID}Rheatmap.pdf',
    "DECOMP_SIGS": '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_{jobID}/{jobID}.sigs.txt',
    "DECOMP_SAMS": '/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_{jobID}/{jobID}.sams.txt',
    "DECOMP_SCRIPT_PATH": '/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/decomp.R',
    "ACCESIBLE_LINK": 'http://pathways-pellegrini.mcdb.ucla.edu/submit/',
    "MATRIX_GENE_DIR": '/UCSC/Pathways-Auxiliary/UCLApathways-Larry-Execs/SigByRank/Matrices/topGenes/'
}
class Constants:
    def __init__(self,jobID):
        self.constants = baseConstants
        for c in self.constants:
            if "{jobID}" in self.constants[c]:
                self.constants[c] = self.constants[c].replace("{jobID}",jobID)
    def __getattr__(self,name):
        if name not in self.constants:
            print "Invalid constant name {0}".format(name)
            exit()
        return self.constants[name]
