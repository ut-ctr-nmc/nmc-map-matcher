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
from nmc_mm_lib import gtfs, path_engine, graph, compat
from io import StringIO
import sys, argparse, xml.sax

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
    print("INFO: Read OpenStreetMap topology...", file=sys.stderr)
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

class Pass1Handler(xml.sax.ContentHandler):
    def __init__(self):
        super().__init__()
        self.nodeIDList = []
        self.highwayFlag = False
        self.wayID = 0
        
        # The outputs:
        self.wayIDSet = set()
        self.nodeIDRefs = {}
        
    def startElement(self, tag, attributes):
        if tag == "way":
            self.highwayFlag = False
            self.nodeIDList.clear()
            self.wayID = int(attributes["id"]) 
        elif self.wayID and tag == "nd":
            self.nodeIDList.append(int(attributes["ref"]))
        elif self.wayID and tag == "tag":
            if attributes["k"] == "highway" or \
                    (attributes["k"] == "oneway" and (attributes["v"] == "yes" or attributes["v"] == "1")):
                self.highwayFlag = True
                
    def endElement(self, tag):
        if self.wayID and tag == "way":
            if self.highwayFlag:
                for nodeID in self.nodeIDList:
                    if nodeID not in self.nodeIDRefs:
                        self.nodeIDRefs[nodeID] = 0
                    self.nodeIDRefs[nodeID] += 1
                self.wayIDSet.add(self.wayID)
            self.wayID = 0

def parseOSMSpeedVal(inStr):
    if inStr.endswith(" mph"):
        return int(inStr[:-4]) # Parse numeric part of mph
    else:
        try:
            return int(int(inStr) * 0.621371) # Convert km to mph
        except ValueError:
            return inStr
           
class Pass2Handler(xml.sax.ContentHandler):
    def __init__(self, nodeIDRefs, wayIDSet):
        super().__init__()
        
        # From Pass 1:
        self.nodeIDRefs = nodeIDRefs
        self.wayIDSet = wayIDSet
        
        # Making GraphLib and Nodes:
        self.latSum = 0
        self.lngSum = 0
        self.rows = 0
        self.nodeCollection = {}
        self.orderingFlag = False
        self.linkID = 1
        
        # Current Way:
        self.wayID = 0
        self.oneWayFlag = False
        self.streetName = ""
        self.speedLimit = None
        self.speedLimitRev = None
        self.nodeRefList = []
        
        # The output:        
        self.graphLib = None
     
    def startElement(self, tag, attributes):
        if tag == "node":
            nodeID = int(attributes["id"])
            if nodeID in self.nodeIDRefs:
                if self.graphLib and not self.orderingFlag:
                    print("WARNING: The OSM XML file is assumed to be structured such that nodes come before ways. The center-point may not be the mean of the entire set.")
                    self.orderingFlag = True
                lat = float(attributes["lat"])
                lng = float(attributes["lon"])
                self.latSum += lat
                self.lngSum += lng
                self.rows += 1
                node = graph.GraphNode(nodeID, lat, lng)
                nodeRef = _NodeRef(node)                
                nodeRef.refCount = self.nodeIDRefs[nodeID]
                self.nodeCollection[nodeID] = nodeRef
        elif tag == "way":
            # We are assuming that all of the ways come after the nodes.
            self.wayID = int(attributes["id"])
            if self.wayID in self.wayIDSet:
                if not self.graphLib:
                    self.graphLib = graph.GraphLib(self.latSum / self.rows, self.lngSum / self.rows)
                self.oneWayFlag = False
                self.streetName = ""
                self.speedLimit = None
                self.nodeRefList.clear()
            else:
                self.wayID = 0
        elif self.wayID and tag == "nd":
            nodeID = int(attributes["ref"])
            if nodeID in self.nodeCollection:
                self.nodeRefList.append(self.nodeCollection[nodeID])
            else:
                print("WARNING: Node %d is attempted to be referenced from Way %d, but it hadn't been defined." % (nodeID, self.wayID), file=sys.stderr)
        elif self.wayID and tag == "tag":
            if attributes["k"] == "oneway" and (attributes["v"] == "yes" or attributes["v"] == "1"):
                self.oneWayFlag = True
            elif attributes["k"] == "name":
                self.streetName = attributes["v"]
            elif attributes["k"] == "maxspeed":
                self.speedLimit = parseOSMSpeedVal(attributes["v"])
            elif attributes["k"] == "maxspeed:forward":
                self.speedLimit = parseOSMSpeedVal(attributes["v"])
            elif attributes["k"] == "maxspeed:backward":
                self.speedLimitRev = parseOSMSpeedVal(attributes["v"])            
                
    def endElement(self, tag):
        if tag == "way" and self.wayID:
            if len(self.nodeRefList) >= 2:
                # Process the beginning:
                prevNodeRef = self.nodeRefList[0]
                vertices = [graph.GraphLinkVertex(prevNodeRef.node.gpsLat, prevNodeRef.node.gpsLng)]
                vertices[-1].id = prevNodeRef.node.id
                if prevNodeRef.node.id not in self.graphLib.nodeMap:
                    self.graphLib.addNode(prevNodeRef.node)
                    
                # Process the middle, creating individual links among places where nodes branch off to multiple paths:
                for nodeIndex in range(1, len(self.nodeRefList)):
                    nodeRef = self.nodeRefList[nodeIndex]
                    vertices.append(graph.GraphLinkVertex(nodeRef.node.gpsLat, nodeRef.node.gpsLng))
                    vertices[-1].id = nodeRef.node.id
                    # Check if we need to end a link:
                    if nodeIndex == len(self.nodeRefList) - 1 or nodeRef.refCount > 1:
                        if nodeRef.node.id not in self.graphLib.nodeMap:
                            self.graphLib.addNode(nodeRef.node)
                        linkID = self.linkID
                        self.linkID += 1
                        # Or, if we want to label according to nodes: "wayID:startNode-endNode"
                        linkID = str(self.wayID) + ":" + str(prevNodeRef.node.id) + "-" + str(nodeRef.node.id)
                        link = graph.GraphLink(linkID, prevNodeRef.node, nodeRef.node, self.graphLib)
                        link.metadata["streetName"] = self.streetName
                        link.metadata["speedLimit"] = self.speedLimit
                        link.addVertices(vertices)
                        self.graphLib.addLink(link)
                        
                        # We also need to add in the complementary reverse link:
                        if not self.oneWayFlag:
                            # Create a copy of the vertex list, as it gets modified in addVertices().
                            revVertices = []
                            for vertex in reversed(vertices):
                                revVertices.append(graph.GraphLinkVertex(vertex.lat, vertex.lng))
                                revVertices[-1].id = vertex.id
                            linkID = str(self.wayID) + ":" + str(nodeRef.node.id) + "-" + str(prevNodeRef.node.id)
                            link = graph.GraphLink(linkID, nodeRef.node, prevNodeRef.node, self.graphLib)
                            link.metadata["streetName"] = self.streetName
                            link.metadata["speedLimit"] = self.speedLimitRev if self.speedLimitRev is not None else self.speedLimit
                            link.addVertices(revVertices)
                            self.graphLib.addLink(link)                            
                        
                        if nodeIndex < len(self.nodeRefList) - 1:
                            # Reset the vertex list for the next segment:
                            prevNodeRef = nodeRef
                            vertices = [graph.GraphLinkVertex(nodeRef.node.gpsLat, nodeRef.node.gpsLng)]
                            vertices[-1].id = nodeRef.node.id
            else:
                print("WARNING: Two or more nodes must be defined for Way %d to be used." % self.wayID, file=sys.stderr)
            self.wayID = 0

def fillGraph(osmXML, isString=False):
    """
    fillGraph fills up the Graph structure from the given OSM XML filename
    @return A Graph representing the underlying topology
    @rtype graph.GraphLib
    """
    
    if not isString:
        print("INFO: Reading geographic input from OSM XML file '%s'..." % osmXML, file=sys.stderr)
    parser = xml.sax.make_parser()
    parser.setFeature(xml.sax.handler.feature_namespaces, 0)
    print("INFO: - Pass 1 of 2...", file=sys.stderr)
    pass1Handler = Pass1Handler()
    osmXMLPass = osmXML if not isString else StringIO(osmXML)
    xml.sax.parse(osmXMLPass, pass1Handler)
    
    print("INFO: - Pass 2 of 2...", file=sys.stderr)
    pass2Handler = Pass2Handler(pass1Handler.nodeIDRefs, pass1Handler.wayIDSet)
    osmXMLPass = osmXML if not isString else StringIO(osmXML)
    xml.sax.parse(osmXMLPass, pass2Handler)
    
    # Optimize the graph for lookups:
    pass2Handler.graphLib.generateQuadSet()
        
    # There we are.
    return pass2Handler.graphLib

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
    print("INFO: -- Final --", file=sys.stderr)
    print("INFO: Print output...", file=sys.stderr)
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
