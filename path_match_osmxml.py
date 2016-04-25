"""
path_match_osmxml.py resolves a GTFS shapefile or trackpoint CSV file to an OpenStreetMap network
    series of links and outputs CSV data that expresses the shapefile points in terms of the underlying
    OSM topology.
    
    NOTE that this is a temporary feature addition until a more elegant multi-format scheme
    is working in the next version.
    
@author: Kenneth Perrine
@contact: kperrine@utexas.edu
@organization: Network Modeling Center, Center for Transportation Research,
    Cockrell School of Engineering, The University of Texas at Austin 
@version: 1.0

@copyright: (C) 2016, The University of Texas at Austin
@license: GPL v3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from __future__ import print_function
from nmc_mm_lib import gtfs, vista_network, path_engine, graph, compat
import sys, argparse, xml.etree.ElementTree as ET

def pathMatch(osmXML, shapePath, limitMap=None, trackpoints=False):
    # Default parameters, with explanations and cross-references to Perrine et al., 2015:
    pointSearchRadius = 1000    # "k": Radius (ft) to search from GTFS point to perpendicular VISTA links
    pointSearchPrimary = 350    # "k_p": Radius (ft) to search from GTFS point to new VISTA links    
    pointSearchSecondary = 200  # "k_s": Radius (ft) to search from VISTA perpendicular point to previous point
    limitLinearDist = 3800      # Path distance (ft) to allow new proposed paths from one point to another
    limitDirectDist = 3500      # Radius (ft) to allow new proposed paths from one point to another
    limitDirectDistRev = 500    # Radius (ft) to allow backtracking on an existing link (e.g. parking lot)
    distanceFactor = 1.0        # "f_d": Cost multiplier for Linear path distance
    driftFactor = 2.0           # "f_r": Cost multiplier for distance from GTFS point to its VISTA link
    nonPerpPenalty = 1.5        # "f_p": Penalty multiplier for GTFS points that aren't perpendicular to VISTA links
    limitClosestPoints = 12     # "q_p": Number of close-proximity points that are considered for each GTFS point 
    limitSimultaneousPaths = 8  # "q_e": Number of proposed paths to maintain during pathfinding stage
    
    maxHops = 12                # Maximum number of VISTA links to pursue in a path-finding operation
    
    # Read in the topology from the OSM XML file:
    print("INFO: Read OpenStreetMap topology from '%s'..." % osmXML, file=sys.stderr)
    vistaGraph = fillGraph(osmXML)
    
    # Read in the shapefile information:
    if not trackpoints:
        print("INFO: Read GTFS shapefile...", file=sys.stderr)
        gtfsShapes = gtfs.fillShapes(shapePath, vistaGraph.gps)
    else:
        # TODO: Complete
        pass
    
    # Initialize the path-finder:
    pathFinder = path_engine.PathEngine(pointSearchRadius, pointSearchPrimary, pointSearchSecondary, limitLinearDist,
                            limitDirectDist, limitDirectDistRev, distanceFactor, driftFactor, nonPerpPenalty, limitClosestPoints,
                            limitSimultaneousPaths)
    pathFinder.maxHops = maxHops
    
    # Begin iteration through each shape:
    shapeIDs = compat.listkeys(gtfsShapes)
    "@type shapeIDs: list<int>"
    shapeIDs.sort()
    gtfsNodesResults = {}
    "@type gtfsNodesResults: dict<int, list<path_engine.PathEnd>>"
    
    if limitMap is not None:
        for shapeID in limitMap:
            if shapeID not in shapeIDs:
                print("WARNING: Limit shape ID %d is not found in the shape file." % shapeID, file = sys.stderr)
    
    for shapeID in shapeIDs:
        "@type shapeID: int"
        
        if limitMap is not None and shapeID not in limitMap:
            continue
        
        print("INFO: -- Shape ID %d --" % shapeID, file = sys.stderr)
        
        # Find the path for the given shape:
        gtfsNodes = pathFinder.constructPath(gtfsShapes[shapeID], vistaGraph)
    
        # File this away as a result for later output:
        gtfsNodesResults[shapeID] = gtfsNodes
    return gtfsNodesResults

class _NodeRef:
    """
    Internal class for tying reference counts to nodes.
    """
    def __init__(self, node):
        self.node = node
        self.refCount = 0

def fillGraph(osmXML):
    """
    fillGraph fills up the Graph structure from the given OSM XML filename
    @return A Graph representing the underlying topology
    @rtype graph.GraphLib
    """
    
    print("INFO: Reading geographic input from OSM XML file '%s'..." % osmXML, file=sys.stderr)
    xmlTree = ET.parse(osmXML)
    xmlRoot = xmlTree.getroot()
    
    # Step 1: Figure out the geographic center of the network and create the Graph:
    nodeCollection = {}
    "@type nodeCollection: dict<str, _NodeRef>"
    rows = 0
    latSum = 0
    lonSum = 0
    for xmlNode in xmlRoot.iterfind("node"):
        nodeID = str(xmlNode.attrib["id"])
        lat = float(xmlNode.attrib["lat"])
        latSum += lat
        lon = float(xmlNode.attrib["lon"])
        lonSum += lon
        rows += 1
        
        node = graph.GraphNode(nodeID, lat, lon)
        nodeCollection[nodeID] = _NodeRef(node)

    if rows == 0:
        print("ERROR: The OSM XML file '%s' doesn't have any node data." % osmXML, file=sys.stderr)
        return None
    graphLib = graph.GraphLib(latSum / rows, lonSum / rows)

    # Step 2: Collect nodes to OSM ways:
    # First, prepare reference counts for nodes to later figure out whether they're intersections.
    for xmlWay in xmlRoot.iterfind("way"):
        highwayFlag = False
        for tag in xmlWay.iterfind("tag"):
            key = tag.attrib("k")
            if key == "highway":
                highwayFlag = True
        if highwayFlag:
            errFlag = False
            nodeCount = 0
            for osmNode in xmlWay.iterfind("nd"):
                if osmNode.attrib("ref") not in nodeCollection:
                    print("WARNING: In OSM Way %s, Node %s is not found." % (xmlWay.attrib("id"),
                        osmNode.attrib("ref")), file=sys.stderr)
                    errFlag = True
                    break
                nodeCount += 1
            if errFlag:
                continue
            if len(nodeCount) < 2:
                print("WARNING: The OSM Way %s has too few elements." % xmlWay.attrib("id"), file=sys.stderr)
                continue
            for osmNode in xmlWay.iterfind("nd"):
                # Increment the reference count for this node:
                nodeCollection[osmNode.attrib("ref")].refCount += 1

    # Next, create links between all Nodes that are endpoints or referenced more than once (which implies
    # intersections or joints between two OSM ways): 
    for xmlWay in xmlRoot.iterfind("way"):
        highwayFlag = False
        onewayFlag = False
        streetName = ""
        for tag in xmlWay.iterfind("tag"):
            key = tag.attrib("k")
            if key == "highway":
                highwayFlag = True
            elif key == "oneway":
                if tag.attrib("v") == "yes" or tag.attrib("v") == "1":
                    onewayFlag = True
            if key == "name":
                streetName = tag.attrib("v")
        if highwayFlag:
            nodeIDList = []
            errFlag = False
            for osmNode in xmlWay.iterfind("nd"):
                nodeIDList.append(osmNode.attrib("ref"))
                if nodeIDList[-1] not in nodeCollection:
                    errFlag = True
                    break
            if errFlag:
                continue
            if len(nodeIDList) < 2:
                continue
            for nodeIndex in range(0, len(nodeIDList)):
                nodeRef = nodeCollection[nodeIDList[nodeIndex]]
                addFlag = False
                if nodeIndex == 0 or nodeIndex == len(nodeIDList) - 1 or \
                        nodeRef.refCount > 1:
                    addFlag = True
                if addFlag and nodeID not in graphLib.nodeMap:
                    graphLib.addNode(nodeRef.node)
            prevNodeID = nodeIDList[0]
            nodeRef = nodeCollection[nodeIDList[0]]
            vertices = [graph.GraphLinkVertex(nodeRef.node.gpsLat, nodeRef.node.gpsLng)]
            for nodeIndex in range(1, len(nodeIDList)):
                nodeRef = nodeCollection[nodeIDList[nodeIndex]]
                vertices.append(graph.GraphLinkVertex(nodeRef.node.gpsLat, nodeRef.node.gpsLng))
                if nodeIndex == len(nodeIDList) - 1 or nodeRef.refCount > 1:
                    linkID = nodeIDList[prevNodeID] + "-" + nodeIDList[nodeIndex]
                    link = graph.GraphLink(linkID, nodeCollection[prevNodeID], nodeCollection[nodeIDList[nodeIndex]])
                    link.streetName = streetName
                    link.addVertices(vertices)
                    prevNodeID = nodeIDList[nodeIndex]
                    graphLib.addLink(link)
                    if nodeIndex < len(nodeIDList) - 1:
                        vertices = [graph.GraphLinkVertex(nodeRef.node.gpsLat, nodeRef.node.gpsLng)]
            if not onewayFlag:
                prevNodeID = nodeIDList[-1]
                nodeRef = nodeCollection[nodeIDList[-1]]
                vertices = [graph.GraphLinkVertex(nodeRef.node.gpsLat, nodeRef.node.gpsLng)]
                for nodeIndex in range(len(nodeIDList) - 2, 0, -1):
                    nodeRef = nodeCollection[nodeIDList[nodeIndex]]
                    vertices.append(graph.GraphLinkVertex(nodeRef.node.gpsLat, nodeRef.node.gpsLng))
                    if nodeIndex == 0 or nodeCollection[nodeIDList[nodeIndex]].refCount > 1:
                        linkID = nodeIDList[prevNodeID] + "-" + nodeIDList[nodeIndex]
                        link = graph.GraphLink(linkID, nodeCollection[prevNodeID], nodeCollection[nodeIDList[nodeIndex]])
                        link.streetName = streetName
                        link.addVertices(vertices)
                        prevNodeID = nodeIDList[nodeIndex]
                        graphLib.addLink(link)
                        if nodeIndex > 0:
                            vertices = [graph.GraphLinkVertex(nodeRef.node.gpsLat, nodeRef.node.gpsLng)]
            
    # Optimize the graph for lookups:
    graphLib.generateQuadSet()
        
    # There we are.
    return graphLib

def main(argv):
    # Initialize from command-line parameters:
    parser = argparse.ArgumentParser(description="resolves a GTFS shapefile or trackpoint CSV file to an OpenStreetMap network " +
        "series of links and outputs CSV data that expresses the shapefile points in terms of the underlying OSM topology.")
    parser.add_argument("osmXML", help="Filename for the OSM XML file")
    parser.add_argument("shapePath", help="Path to GTFS files or trackpoint file")
    parser.add_argument("-L", "--interlinks", action="store_true", help="Outputs intermediate links and positions")
    parser.add_argument("-T", "--trackpoints", action="store_true", help="Read trackpoint CSV file rather than a GTFS shapefile")
    args = parser.parse_args()

    gtfsNodesResults = pathMatch(args.osmXML, args.shapePath, trackpoints=args.trackpoints)
    
    # Extract useful information:
    print("INFO: -- Final --", file = sys.stderr)
    print("INFO: Print output...", file = sys.stderr)
    path_engine.dumpStandardHeader()

    shapeIDs = compat.listkeys(gtfsNodesResults)
    "@type shapeIDs: list<int>"
    shapeIDs.sort()
    for shapeID in shapeIDs:
        "@type shapeID: int"
        path_engine.dumpStandardInfo(gtfsNodesResults[shapeID])
        
# Boostrap:
if __name__ == '__main__':
    main(sys.argv)
