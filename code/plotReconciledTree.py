from ete3 import Tree, TreeStyle, NodeStyle
import sys



if (len(sys.argv) < 2):
    print("\n\tUsage: python plotReconciledTree.py treeFile outputFile")
    print("\tWith treeFile a Newick tree with internal node annotations,")
    print("\t outputFile is the name of the output pdf file you want, without the pdf extension.\n")
    exit(-1)


###########################################
###########################################
################### MAIN ##################
###########################################
###########################################

treefile=sys.argv[1]
outputfile=sys.argv[2]

try:
    f=open(treefile, 'r')
except IOError:
    print ("Unknown file: ", treefile)
    sys.exit()


conse = ""
for l in f:
    conse = conse+l
f.close()

t = Tree(conse, format=1 )
t.write( features=[], outfile="temp")




ts = TreeStyle()
#ts.scale =  20 # pixels per branch length unit
ts.min_leaf_separation= 10
ts.show_scale = True
ts.rotation = 0 #90

# Base style for duplicatioon nodes
duplication_width = 1
nstyle_duplication = NodeStyle()
nstyle_duplication["size"] = 2.0
nstyle_duplication["fgcolor"] = "#FF0000"
#nstyle_duplication["vt_line_color"] = "#ff0000"
#nstyle_duplication["hz_line_color"] = "#ff0000"
nstyle_duplication["vt_line_width"] = duplication_width
nstyle_duplication["hz_line_width"] = duplication_width
nstyle_duplication["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
nstyle_duplication["hz_line_type"] = 0

vt_line_width =  1
hz_line_width =  1


# Base style for nodes
nstyle = NodeStyle()
nstyle["size"] = 0.0
nstyle["fgcolor"] = "#000000"
#nstyle["vt_line_color"] = "#ff0000"
#nstyle["hz_line_color"] = "#ff0000"
nstyle["vt_line_width"] = vt_line_width #8
nstyle["hz_line_width"] = hz_line_width
nstyle["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
nstyle["hz_line_type"] = 0


t.sort_descendants()
for n in t.traverse():
    if n.Ev == "D":
        n.set_style(nstyle_duplication)
    else:
        n.set_style(nstyle)
coord = t.render(outputfile + ".pdf", tree_style=ts)#, w=80, h=60, units="mm")
