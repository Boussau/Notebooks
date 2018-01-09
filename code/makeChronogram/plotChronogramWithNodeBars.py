from ete3 import Tree, TreeStyle, NodeStyle, AttrFace
import math
from PIL import Image, ImageDraw, ImageFont
import sys

# Finds the coordinates of the node
def findCoordinatesOfNodeWithId(nodeId):
    for c in coord['nodes']:
        if c[len(c)-2]==nodeId:
            return (c[0],c[1])


def findCoordinatesOfNode(nodeName, nameToIdMap):
    nodeId = nameToIdMap[nodeName]
    return (findCoordinatesOfNodeWithId(nodeId) )


# To get the 95% HPD: mynode.age_quant_5_95
def get95PcHPD(node):
    lowup = "{0.0_0.0}"
    try:
        lowup = node.age_quant_5_95
    except AttributeError:
        print ("\tNode "+str(node._nid) +" has no age_quant_5_95 property")
    low = float(lowup.split('_')[0].replace("{",""))
    up = float(lowup.split('_')[1].replace("}",""))
    return (low,up)


def getNodeAge(node):
    age=0.0
    try:
        age = float(node.age_median)
    except AttributeError:
        print ("\tNode "+str(node._nid) +" has no age_median property")
    return (age)


def plotAnEllipse(x, y, imagedraw):
    width = 4
    imagedraw.ellipse([x-width,y-width,x+width,y+width], fill=(255,255,0,255))
    return


def plotConfidenceInterval(x0, x1, y, drawTrans):
    #imagedraw.line((x0,y, x1,y), fill=(0,0,255,128), width=4)
    drawTrans.line((x0, y, x1, y), fill=(0,0,255,128), width=4)
    #drawTrans.line((x0,y, x1,y), fill=(100), width=4)
    return


# Here we have to do some maths to convert ages into X coordinates.
# minAge corresponds to maxX
# maxAge corresponds to minX
# for age t, we have coordinate x = maxX - (t - minAge)/(maxAge-minAge)*(maxX-minX)
def computeXCoordinate(minX, maxX, minAge, maxAge, t):
    x = maxX - (t - minAge)/(maxAge-minAge)*(maxX-minX)
    return x


def getTickCoordinatesForTimeline(minX, maxX, minAge, maxAge, mint_timeline, maxt_timeline, every, offsetX):
    current = mint_timeline
    coords = dict()
    while current <= maxt_timeline:
        coord = computeXCoordinate(minX, maxX, minAge, maxAge, current)+ offsetX
        coords[coord] = current
        current = current + every
    return coords


def plotTimeline(minX, maxX, y, minAge, maxAge, ticks, draw):
    draw.line((minX,y, maxX,y), fill=(0,0,0,255), width=1)
    draw.line((minX,y, minX,y+5), fill=(0,0,0,255), width=1)
    draw.line((maxX,y, maxX,y+5), fill=(0,0,0,255), width=1)
    font = ImageFont.truetype('Verdana.ttf', 10)
    #font = ImageFont.load_default()
    for coord in ticks.keys():
        draw.text((coord,y-12),'{:.0f}'.format(ticks[coord]),(0,0,0),font=font)
        draw.line((coord,y, coord,y+5), fill=(0,0,0,255), width=1)
#    draw.text((minX,y-12),'{:.0f}'.format(minAge),(0,0,0),font=font)
#    draw.text((maxX,y-12),'{:.0f}'.format(maxAge),(0,0,0),font=font)
    return



if (len(sys.argv) < 7):
    print("\n\tUsage: python plotChronogramWithNodeBars.py speciesTreeFile outputFile mint_timeline maxt_timeline every codefile [xscale]\n")
    print("\tWith speciesTreeFile a NHX tree with internal node HPD interval (annotation age_quant_5_95),\n")
    print("\toutputFile is the name of the output png file you want, without the png extension,\n")
    print("\tmint_timeline is the minimum age to display on the timeline,\n")
    print("\tmaxt_timeline is the maximum age to display on the timeline,\n")
    print("\tevery is how often ticks will be put on the timeline,\n")
    print("\tcodefile is file giving the correspondence between code names and full species names,\n")
    print("\t[xscale] is the optional scaling factor for branch lengths.\n")
    print("\t Instead, "+str(len(sys.argv))+" arguments were given.\n")
    exit(-1)

treefile=sys.argv[1]
outputfile=sys.argv[2]
mint_timeline = float(sys.argv[3])
maxt_timeline = float(sys.argv[4])
every = float(sys.argv[5])
codefile = sys.argv[6]
scale = 0.3
if len(sys.argv) >7:
    scale = float(sys.argv[7])

widthImage = 2*253

try:
    fco=open(codefile, 'r')
except IOError:
    print ("Unknown file: ", codefile)
    sys.exit()

codeToName = dict()
nameToCode = dict()
for l in fco:
    li = l.strip().split("\t")
    codeToName[li[0]] = li[1]
    nameToCode[li[1]] = li[0]


try:
    f=open(treefile, 'r')
except IOError:
    print ("Unknown file: ", treefile)
    sys.exit()

lines = ""
for l in f:
    lines = lines+l
f.close()

last_comments = lines.rfind("[")
#print(lines[0:last_comments]+";")
t = Tree( lines[0:last_comments]+";" )
#print(t)
# Parse the node information for the root

root_info = lines[last_comments:]
median = float(root_info.split("age_median=")[1].split(":")[0])
mean = float(root_info.split("age_mean=")[1].split("]")[0])
sd = float(root_info.split("age_sd=")[1].split(":")[0])
min = float(root_info.split("age_range={")[1].split("_")[0])
max = float(root_info.split("age_range={")[1].split("_")[1].split("}")[0])
ciMin = float(root_info.split("age_quant_5_95={")[1].split("_")[0])
ciMax = float(root_info.split("age_quant_5_95={")[1].split("_")[1].split("}")[0])
id = float(root_info.split("id=")[1].split(":")[0])
t.add_features(support=1.0, age_median=median, age_mean=mean, age_sd=sd, age_range="{"+str(min)+"_"+str(max)+"}", age_quant_5_95="{"+str(ciMin)+"_"+str(ciMax)+"}", id=id)

ts = TreeStyle()
ts.min_leaf_separation= 0
ts.show_scale=False
ts.show_leaf_name = False
ts.scale = scale #0.3 #0.1# 0.3  # 10 pixels per branch length unit
nstyle = NodeStyle()
nstyle["size"] = 0.0001
for n in t.traverse():
    n.set_style(nstyle)

for leaf in t:
   leaf.add_face(AttrFace("name", fstyle="italic"), column=0, position='branch-right')


#This first rendering is used to get _nid attributes for each node
coord = t.render(outputfile+".png", tree_style=ts, w=widthImage, units="mm")

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

for node in t.traverse("postorder"):
    if node.name == "" or isfloat(node.name):
        node.name = "Internal_"+str(node._nid)
    else:
        node.name = codeToName[node.name]

coord = t.render(outputfile+".png", tree_style=ts, w=widthImage, units="mm")

# Map between node name and node id
nodeNameToId = dict()
nodeIdToXY = dict()
nodeIdToAge = dict()
nodeIdToLowup = dict()
maxX = 0.0
minX = 1000000.0
minAge = 0.0
maxAge = 0.0
maxAgeConf = 0.0
maxY = 0.0
minY = 1000000.0
for node in t.traverse("postorder"):
    nodeNameToId[node.name] = node._nid
    xy = findCoordinatesOfNodeWithId(node._nid)
    if xy[0] > maxX:
        maxX = xy[0]
    if xy[0] < minX:
        minX = xy[0]
    if xy[1] > maxY:
        maxY = xy[1]
    if xy[1] < minY:
        minY = xy[1]
    nodeIdToXY[node._nid] = xy
    #print("XY: "+str(node._nid)+" : " +str(xy))
    lowup = get95PcHPD(node)
    nodeIdToLowup[node._nid] = lowup
    #print("LOWUP: "+str(node._nid) +" : "+ str(lowup))
    age = getNodeAge(node)
    nodeIdToAge[node._nid] = age
    if age > maxAge:
        maxAge = age
    if age < minAge:
        minAge = age
    if lowup[1] > maxAgeConf:
        maxAgeConf = lowup[1]


#print("maxAge: "+str(maxAge))
#print("maxAgeConf: "+str(maxAgeConf))

coord = t.render(outputfile+".png", tree_style=ts, w=widthImage, units="mm")

im = Image.open(outputfile+".png")
siz = im.size

#enlarger = math.ceil(computeXCoordinate(minX, maxX, minAge, maxAge, int(maxAgeConf - maxAge)))
enlarger = 600
sizLarger = [siz[0]+enlarger,siz[1]+30]


# For transparency
#mask=Image.new("RGBA", sizLarger, color=(0,0,0,0))
#drawTrans=ImageDraw.Draw(mask)

offsetX = 0.0
for node in t.traverse("postorder"):
    if not node.is_leaf():
        nid = node._nid
        xy = nodeIdToXY[nid]
        #plotAnEllipse(xy[0], xy[1], draw)
        lowup = nodeIdToLowup[node._nid]
        lowX = computeXCoordinate(minX, maxX, minAge, maxAge, lowup[0])
        upX = computeXCoordinate(minX, maxX, minAge, maxAge, lowup[1])
        if upX < offsetX:
            offsetX = upX
        #plotConfidenceInterval(lowX , upX , xy[1], drawTrans)
offsetX = - math.ceil(offsetX)
#print("offset: "+ str(offsetX))

#offsetX = 300 #math.ceil(computeXCoordinate(minX, maxX, minAge, maxAgeConf, int(maxAgeConf - maxAge)))
#print("offset: "+ str(offsetX))


# Plotting the timeline
y_timeline = sizLarger[1]-10
minXTimeline = computeXCoordinate(minX, maxX, minAge, maxAge, mint_timeline) + offsetX+int(enlarger/2) #+ enlarger
maxXTimeline = computeXCoordinate(minX, maxX, minAge, maxAge, maxt_timeline) + offsetX+int(enlarger/2) #+ enlarger
ticks = getTickCoordinatesForTimeline(minX, maxX, minAge, maxAge, mint_timeline, maxt_timeline, every, offsetX+int(enlarger/2))



canvas = Image.new("RGB", (sizLarger[0]-offsetX, sizLarger[1]), color="white") #(1,1,1))
canvas.paste(im, (offsetX+int(enlarger/2), 0))#, 0.5)
draw = ImageDraw.Draw(canvas)
plotTimeline(minXTimeline, maxXTimeline, y_timeline, mint_timeline, maxt_timeline, ticks, draw)

# For transparency
mask=Image.new("RGBA", (sizLarger[0]-offsetX, sizLarger[1]), color=(0,0,0,0))
drawTrans=ImageDraw.Draw(mask)


for node in t.traverse("postorder"):
    if not node.is_leaf():
        nid = node._nid
        xy = nodeIdToXY[nid]
        #plotAnEllipse(xy[0], xy[1], draw)
        lowup = nodeIdToLowup[node._nid]
        lowX = computeXCoordinate(minX, maxX, minAge, maxAge, lowup[0]) + offsetX+int(enlarger/2)
        upX = computeXCoordinate(minX, maxX, minAge, maxAge, lowup[1]) + offsetX+int(enlarger/2)
        #plotAnEllipse(0+int(enlarger/2),  xy[1], drawTrans)
        plotConfidenceInterval(lowX , upX , xy[1], drawTrans)


canvas.paste(mask, (0, 0), mask=mask)  #, 0.5)

canvas.save(outputfile+"WithIntervals.png", "png")
