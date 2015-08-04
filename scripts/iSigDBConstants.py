baseConstants = {
    "R_TXT": '/home/mike/workspace/PellegriniResearch/scripts/scratch/rOutput.txt',
    "RSCRIPT": 'Rscript',
    "HEATSIG": '/home/mike/workspace/PellegriniResearch/scripts/heatsigV4.R',
    "R_DOWNLOAD": "/home/mike/workspace/PellegriniResearch/scripts/scratch/rDownload.txt",
    "HEATMAP_TEMPLATE": "/home/mike/workspace/PellegriniResearch/scripts/heatmapTemplate.html",
    "COMP_OUTPUT": '/home/mike/workspace/PellegriniResearch/scripts/scratch/output.txt',
    "NULL_PDF": '/home/mike/workspace/PellegriniResearch/output/nulldist.pdf',
    "INPUT_DIST_PDF": '/home/mike/workspace/PellegriniResearch/output/inputdist.pdf',
    "R_HEATMAP": '/home/mike/workspace/PellegriniResearch/output/Rheatmap.pdf',
    "DECOMP_SIGS": '/home/mike/workspace/PellegriniResearch/scripts/scratch/sigs.txt',
    "DECOMP_SAMS": '/home/mike/workspace/PellegriniResearch/scripts/scratch/sams.txt',
    "DECOMP_SCRIPT_PATH": '/home/mike/workspace/PellegriniResearch/scripts/decomp.R',
    "ACCESIBLE_LINK": 'http://pathways-pellegrini.mcdb.ucla.edu/submit/',
    "MATRIX_GENE_DIR": '/home/mike/workspace/PellegriniResearch/sigdir/MATRICES/topGenes/'  
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
