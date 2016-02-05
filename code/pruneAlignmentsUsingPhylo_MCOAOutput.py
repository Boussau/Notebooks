#! /usr/bin/env python
#input: 
# GeneFamiliesToDiscard
# GenesToDiscard
# directory where all alignments are placed.



#Select sequences from a Fasta file according to a list of names of sequences to remove (listsp list)
def copyAlnAndRemoveSeqs(file, listsp, out):
    try:
        fout=open(out, 'w')
    except IOError:
        print ("Unknown file: ",out)
        sys.exit()
    try:
        f=open(file, 'r')
    except IOError:
        print ("Unknown file: ",f)
        sys.exit()
    for l in f:
        if '>' in l:
            if (l.replace('>','').strip() not in listsp):
                tokeep=True
                fout.write(l)
            else:
                tokeep=False
        elif tokeep:
            fout.write(l)
    f.close()   
    fout.close()
    return


import math
import sys
print (sys.argv)
import os.path


argv=sys.argv[1:]
fam=argv[0]
genes=argv[1]
dir = argv[2]
out=argv[0]+".sel"


try:
    fa=open(fam, 'r')
except IOError:
    print ("Unknown file: ",fam)
    sys.exit()

familiesToDiscard = list()
for l in fa:
    familiesToDiscard.append(int(l.strip()))
fa.close()

try:
    fg=open(genes, 'r')
except IOError:
    print ("Unknown file: ",genes)
    sys.exit()

genesToDiscard = dict()
for l in fg:
    if 'Species' in l:
        pass
    elif (genesToDiscard.__contains__(l.split()[1]) ):
        genesToDiscard[l.split()[1]].append(l.split()[0])
    else:
        genesToDiscard[l.split()[1]] = list()
        genesToDiscard[l.split()[1]].append(l.split()[0])
fg.close()

print (familiesToDiscard)
print (genesToDiscard)


num = 0
for file in os.listdir(dir):
    #    print file
    num = num + 1
    if num in familiesToDiscard: 
        print ("Discarding family "+file)
        pass
    else:
        if ( genesToDiscard.__contains__(str(num)) ):
            print ("Discarding genes in " + file)
            copyAlnAndRemoveSeqs(os.path.join(dir,file), genesToDiscard[str(num)], os.path.join(dir,file)+"_cleaned")
        else:
            empty = list()
            print ("Keeping all genes in " + file)
            copyAlnAndRemoveSeqs(os.path.join(dir,file), empty, os.path.join(dir,file)+"_cleaned")







#       
#
#try:
#    fout=open(out, 'w')
#except IOError, e:
#    print "Unknown file: ",out
#    sys.exit()
#    
#try:
#    fref=open(ref, 'r')
#except IOError, e:
#    print "Unknown file: ",ref
#    sys.exit()
#
#listsp=list()    
#for l in fref:
##    listsp.append(l.replace(" ",""))
#    listsp.append(l.strip())
#
#tokeep=False    
#for l in f:
#    if '>' in l:
##        print l
##    	if (l.replace('>','').replace(" ","") in listsp):
#    	if (l.replace('>','').strip() in listsp):
#            tokeep=True
#            fout.write(l)
#    	else:
#            tokeep=False
#            
#    elif tokeep:
#    	fout.write(l)
#f.close()   
#fout.close()
#fref.close()
#    
