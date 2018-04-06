import sys
from ete3 import Tree
import re

if (len(sys.argv)<3):
    print("Usage:\tpython convertPhyloBayesChronogramForFigtree.py infile outfile\n")

in_file = sys.argv[1]
out_file = sys.argv[2]

with open(in_file,  "rU") as handle, open(out_file, "w") as output_handle:
    line = ""
    for l in handle:
        line += l.strip()
    if (line.count(";") != 1):
        print("There should be one and only one tree in the input file "+in_file)
        print("Exiting.")
        exit(1)
    line2 = re.sub('(\d*)(\.*)(\d*)_(\d*)(\.*)(\d*)', "[&height_95%_HPD={\g<1>\g<2>\g<3>,\g<4>\g<5>\g<6>}]", line)
    line2 = re.sub('\](\:)(\d*)(\.*)(\d*)', ", length=\g<2>\g<3>\g<4>]:\g<2>\g<3>\g<4>", line2)

    print (line2)
    line3 = re.sub('(\d*)(\.)*(\d*)_(\d*)(\.)*(\d*)', "", line)

    output_handle.write("#NEXUS\nbegin taxa;")

    t = Tree( line3 )
    names = t.get_leaf_names ()
    output_handle.write("dimensions ntax="+str(len(names)) +";\ntaxlabels\n")
    for n in names:
        output_handle.write(n+"\n")
    output_handle.write(";\nend;\n\nbegin trees;\ntree TREE1 = [&R] [&R=true]\n")
    output_handle.write(line2)
    output_handle.write("end;")
#)735.09_300.332
#[&height_95%_HPD={11.150000000000087,11.150000000000095}]
#length=0.5455417729722768
