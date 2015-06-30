"""
osm_xml.py is an underlying node-link data representation of a
    roadway network coded within OpenStreetMap XML files.
@author: Kenneth Perrine
@contact: kperrine@utexas.edu
@organization: Network Modeling Center, Center for Transportation Research,
    Cockrell School of Engineering, The University of Texas at Austin 
@version: 1.1

@copyright: (C) 2015, The University of Texas at Austin
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
from nmc_mm_lib import geo_data_provider, graph
import sys, xml.etree.ElementTree as ET

class _NodeRef:
    """
    Internal class for tying reference counts to nodes.
    """
    def __init__(self, node):
        self.node = node
        self.refCount = 0
            
class OSMxml(geo_data_provider.GeoDataProvider):
    """
    @ivar osmFilename: str representing the OpenStreetMap XML filename
    """    
    def __init__(self):
        self.osmFilename = ""
        
        super(geo_data_provider.GeoDataProvider, self).__init__()
    
    def provideCmdLineOpts(self, options):
        """
        Stuffs new command-line options into the given argparse object.
        @type options: argparse
        """
        options.add_argument("-fo", nargs=1, type=str, metavar="OSMFILENAME", dest="osmFilename",
            help="If using OpenStreetMap XML for underlying topology input, specifies the XML file.")

    def checkCmdLineOpts(self, args):
        """
        Checks the given command-line options and determines if there are sufficient options
        for successfully loading in data for this class.
        @type args: argparse.Namespace
        @return True if sufficient command-line options are provided for this class
        """
        # Use "or" here to make better error messages later if one is missing.
        return args.osmFilename is not None
    
    def readCmdLineOpts(self, args):
        """
        Reads in all relevant command-line options and emits an error if there are insufficient
        options provided.
        @type args: argparse.Namespace
        @return False if there are insufficient options provided, or True if successful.
        @rtype boolean
        """
        self.osmFilename = args.osmFilename
        return True
        
    def readData(self):
        """
        Using information from the command-line options, reads in data from the geographic data
        source that's handled by this class.
        @return A graph representing the geographic data
        @rtype graph.GraphLib
        """
        print("INFO: Reading geographic input from OpenStreetMap XML file '%s'..." % self.osmFilename, file=sys.stderr)
        xmlTree = ET.parse(self.osmFilename)
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
            print("ERROR: The OSM XML file '%s' doesn't have any node data." % self.osmFilename, file=sys.stderr)
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
                    addFlag = False
                    if nodeIndex == 0 or nodeIndex == len(nodeIDList) - 1 or \
                            nodeCollection[nodeIDList[nodeIndex]].refCount > 1:
                        addFlag = True
                    if addFlag and not graphLib.hasNodeID(nodeID):
                        graphLib.addNode(nodeCollection[node])
                prevNodeID = nodeIDList[0]
                for nodeIndex in range(1, len(nodeIDList)):
                    if nodeIndex == len(nodeIDList) - 1 or \
                            nodeCollection[nodeIDList[nodeIndex]].refCount > 1:
                        linkID = nodeIDList[prevNodeID] + "-" + nodeIDList[nodeIndex]
                        link = graph.GraphLink(linkID, nodeCollection[prevNodeID], nodeCollection[nodeIDList[nodeIndex]])
                        link.streetName = streetName
                        # TODO: Support link curvature.
                        prevNodeID = nodeIDList[nodeIndex]
                        graphLib.addLink(link)
                if not onewayFlag:
                    prevNodeID = nodeIDList[-1]
                    for nodeIndex in range(len(nodeIDList) - 2, 0, -1):
                        if nodeIndex == 0 or nodeCollection[nodeIDList[nodeIndex]].refCount > 1:
                            linkID = nodeIDList[prevNodeID] + "-" + nodeIDList[nodeIndex]
                            link = graph.GraphLink(linkID, nodeCollection[prevNodeID], nodeCollection[nodeIDList[nodeIndex]])
                            link.streetName = streetName
                            # TODO: Support link curvature.
                            prevNodeID = nodeIDList[nodeIndex]
                            graphLib.addLink(link)
            
        # There we are.
        return graphLib
