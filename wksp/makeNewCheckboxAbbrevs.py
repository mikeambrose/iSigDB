#!/usr/env/bin python

abbrevs = open('/home/mike/workspace/PellegriniResearch/sigdir/abbrevs.txt').read().split('\n')
long_to_short = {}
for line in abbrevs:
    if not line:
        continue
    line = line.split('\t')
    long_to_short[line[1]] = line[0]

new_names_file = open('/home/mike/workspace/PellegriniResearch/wksp/newNameToOld.txt').read().split('\n')
new_names = []
old_names = []
in_new = True
for line in new_names_file:
    if not line:
        continue
    if '~' in line:
        in_new = False
    else:
        if in_new:
            new_names.append(line)
        else:
            old_names.append(line)
assert len(old_names) == len(new_names)
new_to_old = {}
for i in range(len(new_names)):
    new_to_old[new_names[i]] = old_names[i]
with open('/home/mike/workspace/PellegriniResearch/sigdir/JSTreeAbbrevs.txt','w') as out:
    for new in new_to_old:
        out.write("{1}\t{0}\n".format(new,long_to_short[new_to_old[new]]))
