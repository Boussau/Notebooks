from ete3 import Tree, TreeStyle, NodeStyle
from PIL import Image, ImageDraw
import sys

if (len(sys.argv) < 4):
    print("\n\tUsage: python plotTransfersOnTree.py speciesTreeFile transferFile outputFile")
    print("\tWith speciesTreeFile a Newick tree with internal node names,")
    print("\tand transferFile a file containing transfers in the following format:")
    print("\tdonor1        receptor1        weight1")
    print("\tdonor2        receptor2        weight2")
    print("\t...")
    print("\tand outputFile is the name of the output png file you want, without the png extension.\n")
    exit(-1)

# Finds the coordinates of the node
def findCoordinatesOfNode(nodeName, nameToIdMap):
    nodeId = nameToIdMap[nodeName]
    for c in coord['nodes']:
        if c[len(c)-2]==nodeId:
            return (c[0],c[1])


# Plots one single transfer on a tree
def plotOrientedTransferOnTree (xyDonor, xyReceptor, imagedraw):
    imagedraw.line((xyDonor[0],xyDonor[1], xyReceptor[0],xyReceptor[1]), fill=(255,0,0,128), width=4)
    width = 4
    imagedraw.ellipse([xyReceptor[0]-width,xyReceptor[1]-width,xyReceptor[0]+width,xyReceptor[1]+width], fill=(255,255,0,255))

# Computes transparencies
def computeTransparencies (amounts):
        maxAmount = max(amounts)
        minAmount = min (amounts)
        # We want minAmount to be at 20, and maxAmount at 255, for the transarency
        minTransp = 10
        maxTransp = 255
        # We look for a linear regression a+b*x
        b =  (maxTransp - minTransp) / (maxAmount-minAmount)
        a = maxTransp - b * maxAmount
        transparencies = list()
        for i in range(len(amounts)):
            transparencies.append(int(a + b * amounts[i]))
        return transparencies

# Plots a list of transfers, with transparency function of the number of genes transferred
def plotSeveralOrientedTransfersOnTree (xyDonors, xyReceptors, amounts, imagedraw):
    transparencies = computeTransparencies (amounts)
    width = 5
    for i in range(len(xyDonors)-1, -1, -1):
        xyDonor = xyDonors[i]
        xyReceptor = xyReceptors[i]
        imagedraw.line((xyDonor[0],xyDonor[1], xyReceptor[0],xyReceptor[1]), fill=(255,0,0,transparencies[i]), width=width)
        imagedraw.ellipse([xyReceptor[0]-width,xyReceptor[1]-width,xyReceptor[0]+width,xyReceptor[1]+width], fill=(255,255,0,255))




treefile=sys.argv[1]
transferfile=sys.argv[2]
outputfile=sys.argv[3]

try:
    f=open(treefile, 'r')
except IOError:
    print ("Unknown file: ", treefile)
    sys.exit()

lines = ""
for l in f:
    lines = lines+l
f.close()

t = Tree(lines, format=1 )

 #"((HalmarSJ:59,(BdeexoJS:7,(BdebacW:6,(Bdebacst:4,Bdebac:4)64:2)65:1)66:52)853:1,((((SorcelSoe7:5,SorcelSo:5)87:11,HalochDS:16)97:1,(((((MyxxanDK:8,MyxfulHW:8)86:2,(MyxstiDS:9,Myxful12:9)85:1)96:1,CorcorDS:11)103:1,Arcgep:12)106:3,(AnaspFw:13,(Anadeh2Ca8:2,(Anadeh2C:1,AnaspK:1)61:1)62:11)63:2)109:2)112:41,(((PelcarDS:31,Geosub:31)81:1,((PelproDS:22,GeolovSZ:22)79:8,((Geopic:24,((GeosulPC:3,GeosulKN:3)82:20,GeometGS:23)94:1)101:5,((GeouraRf:27,GeodalFR:27)78:1,(GeospM1:18,(GeospM2:14,GeobemBe:14)77:4)80:10)93:1)104:1)107:2)110:25,(((SynaciSB:54,DestieDS:54)75:1,((((DesoleHx:35,(DestolTo:34,DesautHR:34)70:1)72:1,DesalkAK:36)89:14,((((DessulDS:37,DespsyLS:37)73:5,DesproDS:42)90:1,DesalkAH:43)98:6,((DesretDS:44,DesbacDS:44)71:4,((((Desvulstn6:26,((Desvulst:21,DesvulRC:21)76:4,DesvulDP:25)92:1)100:12,((LawintPH:19,LawintN3:19)84:14,Desdessu:33)95:5)102:1,DesalaG2:39)105:8,(DesmagRS:46,(Desafrst:45,(DessalDS:41,(DesdesND:40,DesaesAs:40)69:1)74:4)91:1)99:1)108:1)111:1)113:1)114:3,(SynfumMP:52,(DesbaaDS:51,DesaceDS:51)68:1)88:1)115:2)116:1,(HipmarDS:20,DesaceA6:20)67:36)117:1)118:1)119:2);", format=1 )


ts = TreeStyle()
ts.min_leaf_separation= 0
#ts.scale = 88
nstyle = NodeStyle()
nstyle["size"] = 0.0001
for n in t.traverse():
    n.set_style(nstyle)

#t.render("tree.png", tree_style=ts)

coord = t.render(outputfile+".png", tree_style=ts, w=183, units="mm")


#coord = t.render("tree.png")
#print(coord)

# Map between node name and node id
nodeNameToId = dict()
for node in t.traverse("postorder"):
  #print (node.name, node._nid)#, node.dist, node.support)
  nodeNameToId[node.name] = node._nid



im = Image.open(outputfile+".png")

draw = ImageDraw.Draw(im)



# Reading the transfers
try:
    f=open(transferfile, 'r')
except IOError:
    print ("Unknown file: ", transferfile)
    sys.exit()

xy1=(0,0)
xy2=(0,0)

donorsXY=list()
receptorsXY=list()
amounts = list()


for l in f:
    li = l.split()
    donor = li[0]
    receptor = li[1]
    donorsXY.append(findCoordinatesOfNode(donor, nodeNameToId))
    receptorsXY.append(findCoordinatesOfNode(receptor, nodeNameToId))
    amounts.append(float(li[2]))

plotSeveralOrientedTransfersOnTree (donorsXY, receptorsXY, amounts, draw)

for l in f:
    li = l.split()
    donor = li[0]
    receptor = li[1]
    donorId = nodeNameToId[donor]
    receptorId = nodeNameToId[receptor]
    xy1 = findCoordinatesOfNode(donor, nodeNameToId)
    xy2 = findCoordinatesOfNode(receptor, nodeNameToId)
    plotOrientedTransferOnTree (xy1, xy2, draw)

f.close()

im.save(outputfile+"WithTransfers.png", "png")


#
#
# donor = "92"
# receptor="DesalaG2"
#
# donorId = nodeNameToId[donor]
# receptorId = nodeNameToId[receptor]
#
# #print("FACES: "+str(len(coord['faces'])))
# #print("NODES: "+str(len(coord['nodes'])))
#
# #faces: positions of the leaves
# #node_areas: top-left to bottom right rectangles containing clades
#
#
#
#
# print(donorId)
# print(receptorId)
# print(xy1)
# print(xy2)
#
#
#
# #xy1 = (181.0, 501.0)
# #xy2 = (225, 511)
# #draw.line((xy1[0],xy1[1], xy2[0],xy2[1]), fill="#ff0000", width=2)
#
# #xy1 = (64.0, 501.0)
# #xy2 = (225, 511)
#
# #draw.line((xy1[0],xy1[1], xy2[0],xy2[1]), fill=128, width=2)
#
# #im.show()
