
def readTreeFromFile(file):
    try:
        f=open(file, 'r')
    except IOError:
        print ("Unknown file: "+file)
        sys.exit()

    line = ""
    for l in f:
        line += l.strip()

    f.close()
    t = Tree( line )
    return t



def readMAPChronogramFromRBOutput (file):
    try:
        f=open(file, 'r')
    except IOError:
        print ("Unknown file: "+file)
        sys.exit()
    line = ""
    treeStrings = list()
    for l in f:
        if "tree TREE1 = [&R]" in l:
            line = l.replace("tree TREE1 = [&R]", "")
            tree = re.sub('\[&index=\d+([,\w=\d,%\.\{\}])*\]', "", line)#[&index=102,posterior=1.000000,ccp=1.000000,height_95%_HPD={0.025722,0.071446}]
            #print(tree)
            return Tree(tree)


def readMAPChronogramFromRBOutputAndExtract95Hpd (file):
    try:
        f=open(file, 'r')
    except IOError:
        print ("Unknown file: "+file)
        sys.exit()
    line = ""
    treeStrings = list()
    for l in f:
        if "tree TREE1 = [&R]" in l:
            line = l.replace("tree TREE1 = [&R]", "")
            # We are going to create a tree string with node indices
            # At the same time, we'll store 95% HPD along with the node indices.
            # We'll return the tree and a map between the node indices and the 95% HPD.
            line2 = re.sub('\[&index=(\d+)]', '', line)
            tree = re.sub('\[&index=(\d+)([,\w=\d,%\.\{\}])+\]', r'\g<1>', line2)#[&index=102,posterior=1.000000,ccp=1.000000,height_95%_HPD={0.025722,0.071446}]
            brackets = re.findall('\[&index=[,\w=\d,%\.\{\}]*\]', line2)
            idToHPD = dict()
            for i in brackets:
                index = re.findall('\[&index=(\d+)', i)
                hpd = re.findall('age_95%_HPD={(\d+\.*\d*,\d+\.*\d*)}', i)
                idToHPD[index[0]] = list()
                idToHPD[index[0]].append(float(hpd[0].split(",")[0]))
                idToHPD[index[0]].append(float(hpd[0].split(",")[1]))
            return Tree(tree), idToHPD

### The two following functions are used to number nodes in a tree according to another tree. This is necessary to compare node ages between two trees.

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

def renumberNodes( treeToAnnotate, leafList2NodeId ):
    #print treeToAnnotate.get_ascii(attributes=[ "name"], show_internal=False)
    node2leaves = treeToAnnotate.get_cached_content()
    for k in node2leaves.keys():
        if len(node2leaves[k]) == 1: #leaf node
            pass
        else:
            nodelist = list()
            for n in node2leaves[k]: # for all the leaves in the subtree
                nodelist.append( n.name )
            nodelist.sort()
            k.support = leafList2NodeId[tuple(nodelist)]


def renumberNodesAndUpdate95HPDAccordingly( treeToAnnotate, leafList2NodeId, idToHPD ):
    newIdToHPD = dict()
    #print treeToAnnotate.get_ascii(attributes=[ "name"], show_internal=False)
    node2leaves = treeToAnnotate.get_cached_content()
    for k in node2leaves.keys():
        if len(node2leaves[k]) == 1: #leaf node
            pass
        else:
            nodelist = list()
            for n in node2leaves[k]: # for all the leaves in the subtree
                nodelist.append( n.name )
            nodelist.sort()
            newName = leafList2NodeId[tuple(nodelist)]
            newIdToHPD[newName] = idToHPD[str(int(k.support))]
            k.support = str(newName)
    return(newIdToHPD)

def getNodeHeights( t ):
    node2Height = dict()
    id2Height = dict()
    for node in t.traverse("postorder"):
        if node not in node2Height:
            node2Height[node] = 0.0
            id2Height[node.support] = 0.0
        if node.up:
            if node.up.name =='':
                leaves = node.up.get_leaves()
                name=""
                for l in leaves:
                    name += l.name
                node.up.name=name
            node2Height[node.up] = node2Height[node] + node.dist
            id2Height[str(int(node.up.support))] = node2Height[node] + node.dist
      # print node.name + " : " + str(node2Height[node])
    #return node2Height,id2Height
    return id2Height


def plotAndComputeCorrelation(x,y,namex, namey, logyn, logxn, limx=None, limy=None):
    print("Pearson correlation coefficient and p-value: "+ str(scipy.stats.pearsonr(x, y)))
    #Plotting:
    fig, ax = plt.subplots(figsize=(20, 10))
    ax.plot(x, y, 'bo')
    # draw diagonal line
    xy = x+y
    maxxy=max(xy)
    print(maxxy)
    ax.plot([0, 2*maxxy], [0, 2*maxxy ], color='k', linestyle='--',  lw=2)
    plt.axis([0, 1.1*maxxy, 0, 1.1*maxxy])
    plt.xlabel(namex, fontsize=15)
    plt.ylabel(namey, fontsize=15)
#plt.legend(['data'], loc='upper left')
    if logyn:
        ax.yscale('log')
    if logxn:
        ax.xscale('log')
    if not limx == None:
        ax.xlim(limx)
    if not limy == None:
        ax.ylim(limy)
    fig.show()


def plotAndComputeSeveralCorrelations(x, y1, y2, y3, y4, namex, namey, namey1, namey2, namey3, namey4, logyn, logxn, limx=None, limy=None):
    print(namey1 + " : Pearson correlation coefficient and p-value: "+ str(scipy.stats.pearsonr(x, y1)))
    print(namey2 + " : Pearson correlation coefficient and p-value: "+ str(scipy.stats.pearsonr(x, y2)))
    print(namey3 + " : Pearson correlation coefficient and p-value: "+ str(scipy.stats.pearsonr(x, y3)))
    print(namey4 + " : Pearson correlation coefficient and p-value: "+ str(scipy.stats.pearsonr(x, y4)))

    #Plotting:
    fig, ax = plt.subplots(figsize=(20, 15))
    cl, = ax.plot(x, y1, 'bh')
    cal, = ax.plot(x, y2, 'gD')
    con, = ax.plot(x, y3, 'rv')
    calcon, = ax.plot(x, y4, 'co')

    # draw diagonal line
    xy = x+y1
    maxxy=max(xy)
    ax.plot([0, 2*maxxy], [0, 2*maxxy ], color='k', linestyle='--',  lw=2)
    plt.axis([0, 1.1*maxxy, 0, 1.1*maxxy])
    plt.xlabel(namex, fontsize=15)
    plt.ylabel(namey, fontsize=15)
    plt.legend([cl, cal, con,calcon], [namey1, namey2, namey3, namey4], loc='upper left')
    if logyn:
        ax.yscale('log')
    if logxn:
        ax.xscale('log')
    if not limx == None:
        ax.set_xlim(limx)
    if not limy == None:
        ax.set_ylim(limy)
    #fig.show()
