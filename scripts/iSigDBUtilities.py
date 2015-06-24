"""
Various utilities called by different functions used for iSigDB
Made by Mike Ambrose, mikeambrose@berkeley.edu
"""
def reformatFile(f):
    """Automatically changes f to be in the correct format
    Removes carriage returns (either with blank spaces if there are \r\n line endings
        or with newlines if there are just \r line endings)
    Replaces commas with semicolons
    Removes empty lines
    Splits genes with // into two lines
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

def displayErrorMessage(message):
    """Prints a brief html snippet which corresponds to an error message with context 'message'
    Exits after printing"""
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

def writeDetailedOutput(sigGenes,samVals,outFile,fullNames):
    """Writes detailed output to outFile
    sigGenes is a dictionary of signature : gene
    samVals is a dictionary of sample : gene : value
    fullNames is a dictionary of abbreviation to full name"""
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
                        out.write('\t'+str(samVals[sam][gene]))
                out.write('\n')
            out.write('\n')

average = lambda lst: sum(lst) / float(len(lst))

def writeRegularOutput(samSigVals,outFile,fullNames):
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
                    out.write(str(average(samSigVals[sam][sig])))
                else:
                    out.write(str(samSigVals[sam][sig]))
            out.write('\n')