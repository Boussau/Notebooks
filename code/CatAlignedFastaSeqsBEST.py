#! /usr/bin/env python
#Opens all the files in the input directory and concatenates their sequences in respect to their names

import math
import sys
import os.path
#print sys.argv

argv=sys.argv[1:]
dir=argv[0]
outfile=argv[1]

hash=dict()
l_old=""
toAvoid = list()

currentLength = 0

#Get all the names found in the sequences
for file in os.listdir(dir):
    try:
        f=open(os.path.join(dir,file), 'r')
    except IOError:
        print ("Unknown file: "+file)
        sys.exit()
    numseq = 0
    for l in f:
        if '>' in l:
            numseq = numseq + 1
            if l_old != "":
                hash[l_old]=""
            l_old=l
    hash[l_old]=""
    if numseq <= 1:
        toAvoid.append(file)
    print("In "+file+" we have "+str(numseq)+" sequences.") 

#Builds a new dictionary to know whether we found all the species desired sequences
TestHash=dict.fromkeys(hash, 0)

totalNum = 0
hashMissingGenes = dict.fromkeys(hash, 0)

for file in os.listdir(dir):
    #print (file)
    totalNum +=1
    test=True
    for f in toAvoid:
        if file == f:
            test = False
            break
    if test:
        try:
            f=open(os.path.join(dir,file), 'r')
        except IOError:
            print ("Unknown file: "+file)
            sys.exit()
        name=""
        seq=""
        length=0
        for l in f:
            if '>' in l:
                if  name != "":
                    hash[name]=hash[name]+seq
                    seq=seq.strip().replace(' ','')
                    length=seq.__len__()
                    TestHash[name]=1
                name=l
                seq=""
            else:
                seq=seq+l.strip().replace(' ','')
        if length != 0:
            seq=seq.strip().replace(' ','')
            hash[name]=hash[name]+seq
            TestHash[name]=1
            f.close()
            for x in TestHash.keys():
                if TestHash[x]==1:
                    TestHash[x]=0
                else:
                    hash[x]=hash[x]+'-'*length
                    hashMissingGenes[x]+=1
            print("Gene Family "+ file + " from "+ str(currentLength) + " to " + str(currentLength + length-1) + ".")
            currentLength = currentLength + length

print("Total number of gene families: "+str(totalNum))
print("Number of missing genes per species:\n")
for x in hashMissingGenes.keys():
    print("\t"+x.strip().replace(">","")+" : " + str(hashMissingGenes[x]))

            
try:
    fout=open(outfile, 'w')
except IOError:
    print ("Unknown file: "+outfile)
    sys.exit()         

for eachseq in hash:
    if (eachseq.strip()=="") :
    	pass
    else :
    	fout.write (eachseq)
    	fout.write (hash[eachseq]+'\n')
fout.close()
