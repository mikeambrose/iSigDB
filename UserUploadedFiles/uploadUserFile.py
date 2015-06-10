#!/usr/bin/env python
"""python script called to handle file upload
written by Mike Ambrose, mikeambrose@berkeley.edu
"""
import os
from datetime import datetime
import subprocess
import hashlib
import random
from optparse import OptionParser

MAX_FILES_PER_IP_PER_DAY = 5
MAX_FILE_SIZE = 1e8 #this is measured in bytes
MAX_FILE_SIZE_HUMAN_READABLE = "100mb"
MAX_LEN_FILENAME = 100
MIN_SPACE_REMAINING_ON_DISK = 1e7 #10gb (this is what is reported by df, kbs)
def uploadFile(userFile,fileName,userIP,fileSource='files/',assocFile='associations.txt',logFile='logs.txt',\
               hashFile='hashes.txt',outputHTMLFile=None,debug=False):
    """This method is called to upload the file given by userFile
     First, we check a few things:
      the client has not uploaded too many files (MAX_FILES_PER_IP_PER_HOUR)
      the file size is not too large (MAX_FILE_SIZE)
      the file is unique
      the user-given file name is unique
      there is enough space left on the disk
      fileName is acceptable (only ascii, length at most MAX_LEN_FILENAME
    If all these are satisfied, we generate a random name for the file
        This protects against a variety of attacks in which an attacker names a file in a way which
        causes an otherwise trustworthy program to execute it
    We then associate the user-given fileName with the randomly generated location of the file
        this association is added to assocFile, a tab-delimited file of user-made file name : real filepath
    We add the incident to the logFile, a tab-delimited file which tracks IP, file name, time uploaded for the
      last day of uploads
    We write userFile to our randomly generated location
    The function then generates HTML which describes whether or not the upload was successful
    That HTML is printed to standard output if outputHTMLFile=None, otherwise written to the given file
    This function returns whether or not the upload is successful
    """
    fileName = fileName.replace('<','&lt').replace('>','&gt')
    if not checkLogfile(userIP,logFile):
        return reportError("You have uploaded too many files today (limit of " + str(MAX_FILES_PER_IP_PER_DAY)\
                     + " per day)",outputHTMLFile,debug)
    filesize = getFilesize(userFile)
    if debug:
        print "filesize:",filesize
    if filesize > MAX_FILE_SIZE:
        return reportError("Your file is too large (limit of " + MAX_FILE_SIZE_HUMAN_READABLE + ")",outputHTMLFile,debug)
    
    hashVal = hashUserFile(userFile)
    if debug:
        print "hashVal:",hashVal
    if not hashVal: #non-ascii value
        return reportError("Your file has non-ascii characters",outputHTMLFile,debug)
    fileUnique,dupeFileName = checkFileUnique(hashVal,hashFile)
    if not fileUnique:
        return reportError("This file has already been uploaded under the name" + dupeFileName,outputHTMLFile,debug)
    if not checkNameUnique(fileName,assocFile):
        return reportError("This filename has already been used",outputHTMLFile,debug) 
    if not checkDiskSpace(logFile):
        return reportError("There is no more space on the disk. Please email us and let us know",outputHTMLFile,debug)
    if len(fileName) > MAX_LEN_FILENAME:
        return reportError("Your name or your filename is too long",outputHTMLFile,debug)
    if not checkAscii(fileName):
        return reportError("No special characters are allowed in the filename",outputHTMLFile,debug)
    if not checkFileFormat(userFile):
        return reportError("Your file does not have the correct format. Make sure that it is tab-delimited, has a constant number of columns, and every entry is a decimal number.",outputHTMLFile,debug)
    #we generate the random filename with the same method as we generate random seeds
    #concatenate a random 4-digit number to datetime
    randSeed = str(random.randint(0,9999))
    while len(randSeed) != 4:
        randSeed = "0" + randSeed
    randName = datetime.now().strftime('%Y%m%d%H%M%S')+randSeed
    writeNickname(randName,fileName,assocFile)
    writeHash(hashVal,fileName,hashFile)
    logIncident(userIP,fileName,filesize,logFile)
    writeFile(userFile,fileSource+randName)
    return reportSuccess(outputHTMLFile,fileName)


def checkLogfile(userIP,logFile):
    """Checks if the user represented by IP has uploaded too many files
    if the number of instances of IP in logFile is more than MAX_FILES_PER_IP_PER_DAY, returns False
    otherwise returns True
    """
    with open(logFile) as logs:
        logs.readline() #header line
        count = 0
        for line in logs:
            if line.split('\t')[0] == userIP:
                count += 1
                if count >= 5:
                    return False
        return True

def getFilesize(userFile):
    """Returns the filesize of userFile using seek"""
    filePos = userFile.tell()
    userFile.seek(0,os.SEEK_END)
    size = userFile.tell()
    userFile.seek(filePos,os.SEEK_SET)
    return size

def checkAscii(fileName):
    """Returns whether or not every character in fileName is ascii"""
    return all(ord(c) < 128 for c in fileName)

def checkDiskSpace(logFile):
    """Returns whether or not there is at least MIN_SPACE_REMAINING_ON_DISK space left on disk
    Assumes that logFile is stored on disk and uses unix command df"""
    df = subprocess.Popen(["df",logFile], stdout=subprocess.PIPE)
    output = df.communicate()[0]
    _,_,_,space,_,_ = output.split('\n')[1].split()
    return int(space) > MIN_SPACE_REMAINING_ON_DISK

def checkFileUnique(hashVal,hashFile):
    """Returns a tuple of
        (whether or not hashVal is not found in hashFile,
         the name of the file if it is already in our database)
    """
    with open(hashFile) as hashes:
        hashes.readline() #header line
        for line in hashes:
            if line.replace("\n","").split('\t')[1] == hashVal:
                return (False, line.split('\t')[0])
        return (True,'')

def checkNameUnique(fileName,assocFile):
    """Returns whether the nickname fileName is already in use"""
    with open(assocFile) as assocs:
        assocs.readline()
        for line in assocs:
            if line.split('\t')[0] == fileName:
                return False
    return True

def checkFileFormat(userFile):
    """Checks that userFile is vaguely in the format we expect
    Checks that it appears to be tab-delimited with a consistent number of columns
    and that all entries are decimals
    """
    filePos = userFile.tell()
    firstLine = userFile.readline()
    numTabs = firstLine.count('\t')
    if numTabs == 0:
        return False
    for line in userFile:
        if line.count('\t') != numTabs:
            userFile.seek(filePos,os.SEEK_SET)
            return False
        for f in line.split('\t')[1:]:
            try:
                float(f)
            except:
                userFile.seek(filePos,os.SEEK_SET)
                return False
    userFile.seek(filePos,os.SEEK_SET)
    return True

def writeNickname(randName,fileName,assocFile):
    """Writes the randName : fileName association in assocFile"""
    with open(assocFile,'a') as assocs:
        assocs.write(fileName+"\t"+randName+"\n")

def writeHash(hashVal,fileName,hashFile):
    """Writes the fileName : hashVal association in hashFile"""
    with open(hashFile,'a') as hashes:
        hashes.write(fileName+"\t"+hashVal+"\n")

def logIncident(userIP,fileName,fileSize,logFile):
    """Logs the upload in logFile"""
    with open(logFile,'a') as logs:
        logs.write(userIP+"\t"+fileName+"\t"+str(fileSize)+"\n")

def writeFile(userFile,fileLoc):
    """Writes userFile to fileLoc
    (it's times like these that make me question my dedication to writing a comment for every function)"""
    with open(fileLoc,'wb') as outputFile:
        outputFile.write(userFile.read())

def hashUserFile(userFile,blocksize=65536):
    """Returns the MD5 hash of the input file
    blocksize controls how much of the file is read at once
    If any line is not all ascii, fails to complete the hash and returns False
    """
    filePos = userFile.tell()
    hasher = hashlib.md5()
    buf = userFile.read(blocksize)
    while len(buf) > 0:
        if not checkAscii(buf):
            return False
        hasher.update(buf)
        buf = userFile.read(blocksize)
    userFile.seek(filePos,os.SEEK_SET)
    return hasher.hexdigest()

def reportError(errorMsg,outputHTMLFile,debug):
    """Creates the HTML page corresponding to the given error message
    Outputs to a file given in outputHTMLFile or to standard output if outputHTMLFile=None"""
    if debug:
        print errorMsg
    html = """<!doctype html>
<html>
    <head>
        <title>DermDB -- File Upload Error</title>
    </head>
    <body>
        An error occurred: """+errorMsg+"""
    </body>
</html>"""
    if outputHTMLFile:
        with open(outputHTMLFile,'w') as htmlOut:
            htmlOut.write(html)
    else:
        print "Content-type: text/html\n\n" + html
    return 1

def reportSuccess(outputHTMLFile,fileName):
    """Creates the HTML page corresponding to a successful upload
    Outputs to a file given in outputHTMLFile or to standard output if outputHTMLFile=None"""
    html = """<!doctype html>
<html>
    <head>
        <title>DermDB -- File Upload Successful</title>
    </head>
    <body>
        File uploaded successfully under the name """ + fileName + """
    </body>
</html>"""
    if outputHTMLFile:
        with open(outputHTMLFile,'w') as htmlOut:
            htmlOut.write(html)
    else:
        print "Content-type: text/html\n\n"+html
    return 0

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i","--input",dest="userFileLoc",help="File which user is uploading")
    parser.add_option("-n","--name",dest="fileName",help="user's name for file")
    parser.add_option("-p","--ip",dest="userIP",help="user's IP address")
    parser.add_option("-s","--source",dest="fileSource",help="where files are stored")
    parser.add_option("-a","--assoc",dest="assocFile",help="location of associations file")
    parser.add_option("-l","--log",dest="logFile",help="location of log file")
    parser.add_option("-x","--hash",dest="hashFile",help="location of hash file")
    parser.add_option("-o","--output",dest="outputHTMLFile",help="output HTML destination")
    (options,args) = parser.parse_args()
    userFile = open(options.userFileLoc)
    uploadFile(userFile,options.fileName,options.userIP,options.fileSource,options.assocFile,options.logFile,\
               options.hashFile,options.outputHTMLFile)
