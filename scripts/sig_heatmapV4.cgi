import cgi
import cgitb
import detailedHeatmap
cgitb.enable()
form = cgi.FieldStorage()
signature = form["signatureName"].value
scaleArg = "scale" in form
invertArg = "invert" in form
rowArg = form["row_metric"].value
colArg = form["col_metric"].value
seed = form["seed"].value
fileSource = "/UCSC/Pathways-Auxiliary/UCLApathways-Scratch-Space/goTeles_tissueDeconvolutionV2_" + seed + '/' + seed + ".matrix.txt.full.txt"
detailedHeatmap.generateCanvas(fileSource,None,signature,invertArg,rowArg,colArg,seed,False)

