

## The two following functions return dictionaries which map node ids to underlying tips, and vice versa
## First, if the nodes do not have names, we put names into the name and the support attributes:
def createNameToLeavesLink( t ):
    node2leaves = t.get_cached_content()
    nodeId2LeafList = dict()
    leafList2NodeId = dict()
    nodeId = 0
    for k in node2leaves.keys():
        if len(node2leaves[k]) == 1: #leaf node
            pass
        else:
            nodelist = list()
            for n in node2leaves[k]:
                nodelist.append( n.name )
            nodelist.sort()
            nodeId2LeafList[str(nodeId)] = nodelist
            leafList2NodeId[tuple(nodelist)] = str(nodeId) # Need to transform the mutable list into an immutable tuple
            k.name = str(float(nodeId))
            k.support = str(nodeId)
            nodeId = nodeId + 1
    return nodeId2LeafList, leafList2NodeId


## If the nodes have names:
def getNameToLeavesLink( t ):
    node2leaves = t.get_cached_content()
    nodeId2LeafList = dict()
    leafList2NodeId = dict()
    for k in node2leaves.keys():
        if len(node2leaves[k]) == 1: #leaf node
            pass
        else:
            nodelist = list()
            for n in node2leaves[k]:
                nodelist.append( n.name )
            nodelist.sort()
            nodeId2LeafList[k.support] = nodelist
            leafList2NodeId[tuple(nodelist)] = int(k.support) # Need to transform the mutable list into an immutable tuple
    return nodeId2LeafList, leafList2NodeId


## Renames the nodes of a tree according to the names of the underlying tips
def renumberNodes( treeToAnnotate, leafList2NodeId ):
    #print treeToAnnotate.get_ascii(attributes=[ "name"], show_internal=False)
    node2leaves = treeToAnnotate.get_cached_content()
    for k in node2leaves.keys():
        if len(node2leaves[k]) == 1: #leaf node
            pass
        else:
            nodelist = list()
            for n in node2leaves[k]:
                nodelist.append( n.name )
            nodelist.sort()
            k.support = leafList2NodeId[tuple(nodelist)]


## Gets node heights in a tree, and returns 2 dictionaries node->heights and id->heights
def getNodeHeights( t ):
    node2Height = dict()
    id2Height = dict()
    for node in t.traverse("postorder"):
        if node not in node2Height:
            node2Height[node] = 0.0
            id2Height[node.name] = 0.0
        if node.up:
            node2Height[node.up] = node2Height[node] + node.dist
            id2Height[str(node.up.support)] = node2Height[node] + node.dist
    #print (str(node.support) + " : " + str(node2Height[node]))
    return node2Height,id2Height


## Given a series of node age constraints, counts those that are satisfied by the tree
def checkConstraints(id2Height, constraints):
    numberOfMetConstraints = 0
    for c in constraints:
        try:
            if (id2Height[str(float(c[0]))] > id2Height[str(float(c[1]))]):
                numberOfMetConstraints = numberOfMetConstraints + 1
        except KeyError:
            print ("Constraint involving nodes " + str(float(c[0])) + " and " + str(float(c[1])) + " cannot be checked.")
    return numberOfMetConstraints


## A function to compute weighted means and medians from a dataframe with columns kItem and scores
def computeWeightedMeanAndMedian(kItemdf):
    dftempSorted = kItemdf.sort_values(by='kItem', ascending=[0])
    kItem = list(dftempSorted['kItem'])
    weights = list(dftempSorted['scores'])
    sum = 0.0
    mean = 0.0
    i = 0
    while sum <= 0.5:
        sum = sum + weights[i]
        mean += weights[i] * kItem[i]
        i = i +1
    median = kItem[i]
    for j in range(i+1, len(kItem)):
        mean += weights[j] * kItem[j]
    return mean, median


## Function to compute a 95% credible interval based on a dataframe with columns kItem and scores
def compute95PcCredibilityIntervalWeightedDF(kItemdf):
    dftempSorted = kItemdf.sort_values(by='kItem', ascending=[1])
    kItem = list(dftempSorted['kItem'])
    weights = list(dftempSorted['scores'])
    sumi = 0.0
    i = 0
    while sumi <= 0.05:
        sumi = sumi + weights[i]
        i = i +1
    mini = kItem[i]
    sumi = 0.0
    i = 0
    while sumi <= 0.05:
        sumi = sumi + weights[len(kItem) - i-1]
        i = i +1
    maxi = kItem[len(kItem) - i-1]
    return mini, maxi


def addRootNodeInformationToNhx(rootNodeString, treTxt):
    t = treTxt.replace(";",rootNodeString)+"\n"
    return t

## Function to write to file a Nexus tree with credibility intervals on its branches. Can be read by Figtree.
## WARNING: There shouldn't be any "0_" in the species names...
def outputsWeightedChronogram (tre, id2Heights, outFileRadical, weights):
    unweightedAges = list()
    rootNodeString=""
    for node in tre.traverse("postorder"):
        id = 0
        if node.name != "":
            id = node.name
        else:
            id = str(node.support)
            node.name=id
        node.id=id
        print("Node: "+str(id)+" : "+ str(node.dist))
        v = list()
        for i in id2Heights:
            v.append(i[str(id)])
        median = np.median(v)
        sd = np.std(v)
        mean = np.mean(v)
        min = np.min(v)
        max = np.max(v)
        dftempNoWeight = pd.DataFrame({ 'kItem':v, 'scores': weights }, dtype='float')
        ciMin = -1.0
        cMax = -1.0
        ciMin, ciMax = compute95PcCredibilityIntervalWeightedDF(dftempNoWeight)
        node.add_features(support=1.0, age_median=median, age_mean=mean, age_sd=sd, age_range="{"+str(min)+","+str(max)+"}", age_quant_5_95="{"+str(ciMin)+","+str(ciMax)+"}", id=id)
        if node.is_root():
            rootNodeString = "[&&NHX:age_quant_5_95={"+str(ciMin)+"_"+str(ciMax)+"}:name="+node.name+":id="+node.id+":age_sd="+str(sd)+":dist="+str(0.0)+":age_range={"+str(min)+"_"+str(max)+"}:support=1.0:age_median="+str(median)+":age_mean="+str(mean)+"]:0.0 ;"
        children = node.get_children()
        for c in children:
            c.dist = node.age_median - c.age_median
        unweightedAges.append(median)
        #print("Node: "+str(id)+" : "+ str(node.dist))


    treTxt = tre.write(features=[])
    treTxt = addRootNodeInformationToNhx(rootNodeString, treTxt)
    try:
        f=open(outFileRadical+".nhx", 'w')
    except IOError:
        print ("Unknown file: "+outFileRadical+".nhx")
        sys.exit()
    f.write(treTxt)
    f.close()

    treFigTree = treTxt.replace("&&NHX:", "&").replace(":age", ",age").replace(":name", ",name").replace(":id", ",id").replace(":support", ",support").replace(":dist", ",dist").replace("0_", "0,").replace("1_", "1,").replace("2_", "2,").replace("3_", "3,").replace("4_", "4,").replace("5_", "5,").replace("6_", "6,").replace("7_", "7,").replace("8_", "8,").replace("9_", "9,")

    fout = outFileRadical+".nxs"

    try:
        f=open(fout, 'w')
    except IOError:
        print ("Unknown file: "+fout)
        sys.exit()

    f.write("#NEXUS\n")
    f.write("BEGIN TAXA;\n")
    f.write("DIMENSIONS NTAX="+str(len(tre))+";\n")
    f.write("TAXLABELS\n")
    for l in tre:
        f.write(l.name+"\n")

    f.write(";\n")
    f.write("END;\n")

    f.write("BEGIN TREES;\n")

    f.write("TREE 0 = [&R] " + treFigTree)

    f.write("\nEND;\n")

    f.close()
