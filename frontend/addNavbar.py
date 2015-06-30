from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i","--input",dest="input",help="file to add navbar to")
parser.add_option("-n","--navbar",dest="navbar",help="location of navbar html")
parser.add_option("-o","--output",dest="output",help="where to output file")
options,args = parser.parse_args()

navbar = open(options.navbar).read()
content = open(options.input).read()

#find the head content and add it
headStart = content.index('<head>')+len('<head>')
headEnd = content.index('</head>')
headContents = content[headStart:headEnd]

#get rid of illegal elements
illegal = [('<a href="#" class="tooltip">','</a>'), ("<a href='#' class='tooltip'>","</a>")]
for elem in illegal:
    print "removing element {0}".format(elem)
    start,end = elem
    while start in content:
        startInd = content.index(start)
        print "removing instance at {0}".format(startInd)
        endInd = startInd+content[startInd:].index(end)
        print "which goes to {0}".format(endInd)
        content = content[:startInd] + content[endInd:]

#find the body content
bodyStart = content.index('<body>')+len('<body>')
bodyEnd = content.index('</body>')
bodyContents = content[bodyStart:bodyEnd]

#find everything after the body
postBodyStart = content.index('</body>') + len('<body>')
postBodyEnd = content.index('</html>')
postBodyContents = content[postBodyStart:postBodyEnd]

navbar = navbar.format(head=headContents,body=bodyContents,postbody=postBodyContents)

with open(options.output,'w') as out:
  out.write(navbar)
