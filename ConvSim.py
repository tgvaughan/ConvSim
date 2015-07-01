#!/usr/bin/env python

from argparse import ArgumentParser, FileType
from sys import argv, exit, stdout, stderr
from random import Random


def getNthElement(alist, el, n):
    """Return index of nth occurrence of element "el" in alist."""

    idx = None
    k = 0
    for i,thisel in enumerate(alist):
        if thisel==el:
            if k<n:
                k += 1
                continue
            else:
                idx = i
                break

    return idx


class Node:

    def __init__(self):
        self.parents = []
        self.children = []
        self.height = 0.0
        self.ancestral = [(0.0,1.0)]
        self.ancestralParents = []
        self.ancestralChildren = []
        self.label = ""

    def addChild(self, child):
        """This is only sensible when child is a lineage node."""
        self.children.append(child)
        child.parents.append(self)
        child.ancestralParents.append(None)

    def setLabel(self, label):
        self.label = str(label)

    def equals(self, node):
        """Tests whether two BINARY trees are equivalent in the sense
        that their node heights and topologies are equivalent."""

        if self.height != node.height:
            return False

        if len(self.children) != len(node.children):
            return False

        if (len(self.children)>0):
            if len(self.children) != 2:
                raise Exception("Can only compare binary trees!")

            leftEqual =  self.children[0].equals(node.children[0])
            rightEqual = self.children[1].equals(node.children[1])
            if leftEqual and rightEqual:
                return True

            leftEqual = self.children[0].equals(node.children[1])
            rightEqual = self.children[1].equals(node.children[0])
            if leftEqual and rightEqual:
                return True

            return False

        else:
            return True

    def deleteLineage(self):
        """Remove this lineage node from the graph, joining any
        children to any parents and placing self.ancestral in the
        ancestralParents dictionary of its child."""

        if len(self.children) != 1 or len(self.parents) != 1:
            raise Exception("deleteLineage() called for non-lineage node.")

        child = self.children[0]
        parent = self.parents[0]

        idx = child.parents.index(self)
        child.parents[idx] = parent
        child.ancestralParents[idx] = self.ancestral

        parent.children.remove(self)
        parent.children.append(child)

    def populateAncestralChildren(self):
        """Add ancestralChildren lists to nodes in graph."""

        self.ancestralChildren = [None]*len(self.children)

        uniqueChildren = list(set(self.children))

        for child in uniqueChildren:
            count = sum([x==child for x in self.children])

            for n in range(count):
                cidx = getNthElement(self.children, child, n)
                pidx = getNthElement(child.parents, self, n)
                
                self.ancestralChildren[cidx] = child.ancestralParents[pidx]

            child.populateAncestralChildren()

    def findHybrids(self, seenNodes, hybridIDs):
        """Identify hybrid nodes in subgraph."""
        
        if self in seenNodes:
            if self not in hybridIDs.keys():
                newID = len(hybridIDs.keys()) + 1
                hybridIDs[self] = newID
            return

        seenNodes.append(self)

        for child in self.children:
            child.findHybrids(seenNodes, hybridIDs)


    def getNewick(self, prevNode=None, seenNodes=None, hybridIDs=None, childSeenCount=0, marginal=False):
        """Return extended Newick string representing subgraph."""

        if prevNode == None:
            hybridIDs = {}
            seenNodes = []
            self.findHybrids([], hybridIDs)

        if self not in seenNodes:
            previouslyUnseen = True
            seenNodes.append(self)
        else:
            previouslyUnseen = False

        newickStr = ""
        if previouslyUnseen and len(self.children)>0:
            newickStr += "("
            childrenSeen = {}
            for i,child in enumerate(self.children):

                if child not in childrenSeen:
                    childrenSeen[child] = 0
                else:
                    childrenSeen[child] += 1

                if i>0:
                    newickStr += ","

                newickStr += child.getNewick(self,
                                             seenNodes,
                                             hybridIDs,
                                             childrenSeen[child],
                                             marginal=marginal)
            newickStr += ")"

        newickStr += self.label

        if self in hybridIDs:
            newickStr += "#" + str(hybridIDs[self])

        # Add annotation describing ancestral material:
        if not marginal:
            if prevNode == None:
                ancestral = self.ancestral
            else:
                idx = 0
                k = 0
                while self.parents[idx] != prevNode or k<childSeenCount:
                    idx = self.parents.index(prevNode, idx+1)
                    k += 1
                ancestral = self.ancestralParents[idx]


            newickStr += '[&ancestral={'
            for i,interval in enumerate(ancestral):
                if i>0:
                    newickStr += ','
                newickStr +=  str(interval[0]) + ',' + str(interval[1])
            newickStr += '}]'

        if prevNode != None:
            branchLength = prevNode.height - self.height
        else:
            branchLength = 0.0

        newickStr += ":" + str(branchLength)

        if prevNode == None:
            newickStr += ";"

        return newickStr


def ancestralUnion(ancestral1, ancestral2):
    """Return the union of the two sets of ancestral material"""

    #print "AU: {} {}".format(ancestral1, ancestral2)

    newIntervals = []

    for interval in ancestral1:
        newIntervals.append(interval)

    for interval in ancestral2:
        newIntervals.append(interval)

    newIntervals.sort()
    #print newIntervals

    while True:
        overlapsFound = False
        newIntervalsP = []
        i = 0
        while i<len(newIntervals):
            if i<len(newIntervals)-1 and newIntervals[i][1]>=newIntervals[i+1][0]:
                newIntervalsP.append((newIntervals[i][0],max(newIntervals[i][1],newIntervals[i+1][1])))
                i += 2
                overlapsFound = True
            else:
                newIntervalsP.append(newIntervals[i])
                i += 1

        newIntervals = newIntervalsP
        if not overlapsFound:
            break

    #print "AU: {}".format(newIntervals)

    return newIntervals


def ancestralPartition(ancestral, convInterval):
    """Return a pair of mutually exclusive sets of ancestral material
    with the division being made at the chosen interval."""

    ancestral1 = []
    ancestral2 = []

    for interval in ancestral:
        if interval[0]<convInterval[0]:
            if interval[1]<convInterval[0]:
                ancestral1.append(interval)
            else:
                ancestral1.append((interval[0],convInterval[0]))
                ancestral2.append((convInterval[0],min(interval[1],convInterval[1])))
                if interval[1]>convInterval[1]:
                    ancestral1.append((convInterval[1],interval[1]))
        else:
            if interval[1]>convInterval[1]:
                if interval[0]>convInterval[1]:
                    ancestral1.append(interval)
                else:
                    ancestral1.append((convInterval[1],interval[1]))
                    ancestral2.append((interval[0],convInterval[1]))
            else:
                ancestral2.append(interval)

    ancestral1.sort()
    ancestral2.sort()

    return((ancestral1,ancestral2))


def simulate(n, rho, theta, delta, debug=False):
    """Perform simulation.

    Be aware that the nodes in activeLineages are
    quite distinct from the nodes that eventually wind up on the graph.
    The former are _represented_ by node objects, but these objects have
    no definite height.  Coalescence leaves anihiliates two lineages leaving
    one lineage and one node remaining, while conversion results in a duplication
    of lineages but leaves only a single node (with two parents) on the graph."""

    rand = Random()

    activeLineages = []
    for i in range(n):
        node = Node()
        node.setLabel(i)
        lineage = Node() # Lineages are separate to nodes on the graph!
        lineage.addChild(node)
        activeLineages.append(lineage)

    t = 0.0

    # Simulation loop
    while len(activeLineages)>1:
        n = len(activeLineages)

        # Coalescence propensity
        cProp = theta*0.5*n*(n-1)
        
        # Recombination/conversion propensity
        rProp = rho*n

        # Choose time of next event
        totProp = cProp + rProp
        t += rand.expovariate(totProp)

        # Select type of event
        if rand.uniform(0,totProp)<cProp:
            # Coalescence

            # Select random pair of lineages:
            lineages = rand.sample(activeLineages, 2)
            
            # Coalesce nodes:
            parent = Node()
            parent.height = t

            parent.addChild(lineages[0])
            parent.addChild(lineages[1])
            parent.ancestral = ancestralUnion(lineages[0].ancestral, lineages[1].ancestral)

            # Replace coalesced nodes by parent node in active lineages:

            activeLineages.remove(lineages[0])
            activeLineages.remove(lineages[1])
            lineages[0].deleteLineage()
            lineages[1].deleteLineage()

            parentLineage = Node()
            parentLineage.addChild(parent)
            parentLineage.ancestral = parent.ancestral
            activeLineages.append(parentLineage)

        else:
            # Recombination/conversion

            # Select lineage at random
            lineage = rand.sample(activeLineages, 1)[0]

            # Select start and end of converted region:
            boundary1 = rand.uniform(0,1)
            if (rand.uniform(0,1)<0.5):
                boundary2 = min(1,boundary1 + rand.expovariate(1/delta))
            else:
                boundary2 = max(0,boundary1 - rand.expovariate(1/delta))
            boundary1, boundary2 = sorted([boundary1,boundary2])

            # Partition ancestral material:
            newAncestrals = ancestralPartition(lineage.ancestral, (boundary1,boundary2))

            # Continue only if conversion has effect:
            if len(newAncestrals[0])>0 and len(newAncestrals[1])>0:

                if debug:
                    print "t={}: Conversion: {} => {} {}".format(t, lineage.ancestral, newAncestrals[0], newAncestrals[1])

                # Set original node height:
                lineage.height = t
                
                # Generate parents:
                parent1 = Node()
                parent1.addChild(lineage)
                parent1.ancestral = newAncestrals[0]

                parent2 = Node()
                parent2.addChild(lineage)
                parent2.ancestral = newAncestrals[1]

                # Now that the lineage node is on the graph, ensure its child has
                # a corresponding ancestralParents entry:
                idx = lineage.children[0].parents.index(lineage)
                lineage.children[0].ancestralParents[idx] = lineage.ancestral

                # Replace original node with parents in active lineages:
                activeLineages.remove(lineage)
                activeLineages.append(parent1)
                activeLineages.append(parent2)

    root = activeLineages[0].children[0]

    # Populate ancestralChildren arrays:
    root.populateAncestralChildren()

    return root


#### MARGINAL TREE SET CONSTRUCTION ###

def getMarginalTrees(root):
    """Obtain marginal trees from ARG."""

    margRoot = Node()
    margRoot.height = root.height
    margRoot.label = root.label
    nextLocus = getMarginalTree(root, margRoot, 0.0, 1.0)
    margRoots = [margRoot]
    eventLoci = [0.0, nextLocus]

    while nextLocus<1.0:
        margRoot = Node()
        margRoot.height = root.height
        margRoot.label = root.label
        nextLocus = getMarginalTree(root, margRoot, nextLocus, 1.0)
        margRoots.append(margRoot)
        eventLoci.append(nextLocus)

    for i in range(len(margRoots)):
        margRoots[i] = cleanMarginalTree(margRoots[i])

    margIDs = getGenealogyIDs(margRoots)

    return (margRoots, eventLoci, margIDs)

def intervalContainingLocus(ancestral, locus):
    for interval in ancestral:
        if locus>=interval[0] and locus<interval[1]:
            return interval
    return None

def getMarginalTree(node, margNode, curLocus, nextLocus):
    """Construct marginal subtree rooted at margNode corresponding
    to alignment bween curLocus and the next locus at which a
    lineage bifurcation occurs.  The location of this bifurcation is
    returned."""

    for i,child in enumerate(node.children):
        interval = intervalContainingLocus(node.ancestralChildren[i], curLocus)
        if interval == None:
            continue

        nextLocus = min(interval[1], nextLocus)

        margChild = Node()
        margChild.height = child.height
        margChild.label = child.label
        margNode.addChild(margChild)

        nextLocus = getMarginalTree(child, margChild, curLocus, nextLocus)

    return nextLocus

def cleanMarginalTree(node):
    """Remove single-child+parent nodes from clade."""

    if len(node.children)==1:
        
        child = node.children[0]
        child.parents = node.parents

        if len(node.parents)>0:
            parent = node.parents[0]
            idx = parent.children.index(node)
            parent.children[idx] = child
            return cleanMarginalTree(child)
        else:
            return cleanMarginalTree(child)

    else:
        for i in range(len(node.children)):
            cleanMarginalTree(node.children[i])
        return node

def getGenealogyIDs(margRoots):
    """Assign unique integer ID to each unique marginal genealogy."""

    margIDs = []
    uniqueGenealogies = {}
    nextID = 0
    for margRoot in margRoots:
        found = False
        for genID in uniqueGenealogies.keys():
            if margRoot.equals(uniqueGenealogies[genID]):
                margIDs.append(genID)
                found = True
                break
        if not found:
            uniqueGenealogies[nextID] = margRoot
            margIDs.append(nextID)
            nextID += 1
    
    return margIDs

#### DEBUGGING ####

def traverse(node, seenNodes, func, depth):
    """Tree/graph traversal."""
    if node in seenNodes:
        return

    seenNodes.append(node)
    func(node, depth)

    for child in node.children:
        traverse(child, seenNodes, func, depth+1)

def nodeInfo(node, depth):
    stderr.write(" "*depth + str(node.height) + "\n")
    stderr.write(" "*depth + "Parents: " + str(node.parents) + "\n")
    stderr.write(" "*depth + str(node.ancestralParents) + "\n")
    stderr.write(" "*depth + "Children: " + str(node.children) + "\n")
    stderr.write(" "*depth + str(node.ancestralChildren) + "\n\n")

###################


if __name__ == '__main__':

    parser = ArgumentParser(description="""Simulate ARG for given rho,
theta and track length (delta).""")

    parser.add_argument("n", help="Number of samples", type=int)
    parser.add_argument("rho", help="Recombination rate", type=float)
    parser.add_argument("theta", help="Coalescence rate", type=float)
    parser.add_argument("delta", help="Mean track length", type=float)

    parser.add_argument("--debug", help="Turn on debugging output.", action='store_true')

    parser.add_argument("-o","--output", metavar="file", type=FileType('w'), default=stdout,
                        help="""Send extended Newick representation of
ARG to given file. (Defaults to stdout.)""")

    parser.add_argument("-m","--marginal", metavar="file", type=FileType('w'),
                        help="Output marginal genealogies to given file.")

    parser.add_argument("-l","--loci", metavar="file", type=FileType('w'),
                        help="Output genealogy transition loci to given file.""")

    args = parser.parse_args(argv[1:])

    # Simulate ARG
    root = simulate(args.n, args.rho, args.theta, args.delta, debug=args.debug)

    # Output full ARG in extended newick
    args.output.write(root.getNewick() + "\n")

    # Print marginal genealogies and/or loci+IDs:
    if args.marginal != None or args.loci != None:
        margRoots, eventLoci, margIDs = getMarginalTrees(root)

        if args.marginal != None:
            for margRoot in margRoots:
                args.marginal.write(margRoot.getNewick(marginal=True) + "\n")

        if args.loci != None:
            args.loci.write("locus genealogyID\n")
            for i in range(len(margIDs)):
                args.loci.write("{} {}\n".format(eventLoci[i], margIDs[i]))

    #traverse(root, [], nodeInfo, 0)
