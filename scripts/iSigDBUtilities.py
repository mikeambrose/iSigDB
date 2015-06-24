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
    if any(ch1,ch2,ch3):
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

