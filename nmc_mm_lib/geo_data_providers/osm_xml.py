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
            nodeCollection[nodeID] = node

        if rows == 0:
            print("ERROR: The OSM XML file '%s' doesn't have any node data." % self.osmFilename, file=sys.stderr)
            return None
        graphLib = graph.GraphLib(latSum / rows, lonSum / rows)

        # Step 2: Create links out of OSM ways:            
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
                        print("WARNING: In OSM Way %s (for '%s'), Node %s is not found." % (xmlWay.attrib("id"),
                            streetName, nodeIDList[-1]), file=sys.stderr)
                        errFlag = True
                        break
                if errFlag:
                    continue
                if len(nodeIDList) < 2:
                    print("WARNING: The OSM Way %s (for '%s') has too few elements." % (xmlWay.attrib("id"),
                        streetName), file=sys.stderr)
                    continue
                for nodeID in nodeIDList:
                    if not graphLib.hasNodeID(nodeID):
                        graphLib.addNode(nodeCollection[node])
                for nodeIndex in range(1, len(nodeIDList)):
                    linkID = nodeIDList[nodeIndex - 1] + "-" + nodeIDList[nodeIndex]
                    link = graph.GraphLink(linkID, nodeCollection[nodeIDList[nodeIndex - 1]], nodeCollection[nodeIDList[nodeIndex]])
                    
                        
                        
    for (pugi::xml_node nd = way.child("nd"); nd; nd = nd.next_sibling("nd"))
        ei.nodes.push_back(nd.attribute("ref").as_uint(0));

            
    for (pugi::xml_node tag = way.child("tag"); tag; tag = tag.next_sibling("tag"))
        ei.attr.insert(std::pair<std::string, std::string>(tag.attribute("k").value(), tag.attribute("v").value()));

            if "highway" in xmlWay.attrib:
                
            
            
            


        
        
        with open(self.nodesFilename, 'r') as nodesFile:
            # Sanity check:
            fileLine = nodesFile.readline()
            if not fileLine.startswith(NODES_HEADER1) and not fileLine.startswith(NODES_HEADER2):
                print("ERROR: The nodes file '%s' doesn't have the expected header:" % self.nodesFilename, file=sys.stderr)
                print('       "%s" or "%s"' % (NODES_HEADER2, NODES_HEADER1), file=sys.stderr)
                return None
            hasType = fileLine.startswith(NODES_HEADER1)
            
            # Step 1: Figure out the geographic center of the network and create the Graph:
            rows = 0
            x = 0
            y = 0
            for fileLine in nodesFile:
                if len(fileLine) == 0:
                    continue
                lineElems = fileLine.split(',')
                if hasType:
                    if int(lineElems[1] != 1):
                        continue
                    x += float(lineElems[2])
                    y += float(lineElems[3])
                else:
                    x += float(lineElems[1])
                    y += float(lineElems[2])
                rows += 1
            
            if rows == 0:
                print("ERROR: The nodes file '%s' doesn't have any data." % self.nodesFilename, file=sys.stderr)
                return None
            graphLib = nmc_mm_lib.graph.GraphLib(y / rows, x / rows)
            
            # Step 2: Rewind, and fill out the nodes:
            # TODO: If there is a monster nodes file, then it may be worth while to save the values
            # from the first pass into an array and avoid the second read-through.
            nodesFile.seek(0)
            nodesFile.readLine()
            for fileLine in nodesFile:
                if len(fileLine) == 0:
                    continue
                lineElems = fileLine.split(',')
                nodeID = int(lineElems[0])
                if hasType:
                    if int(lineElems[1] != 1):
                        continue
                    x = float(lineElems[2])
                    y = float(lineElems[3])
                else:
                    x = float(lineElems[1])
                    y = float(lineElems[2])
                node = nmc_mm_lib.graph.GraphNode(nodeID, y, x)
                graphLib.addNode(node)

        # Step 3: Fill out the links:
        with open(self.linksFilename, 'r') as linksFile:
            # Sanity check:
            fileLine = linksFile.readline()
            if not fileLine.startswith(LINKS_HEADER1) and not fileLine.startswith(LINKS_HEADER2):
                print("ERROR: The links file '%s' doesn't have the expected header:" % self.linksFilename, file=sys.stderr)
                print('       "%s" or "%s"' % (LINKS_HEADER2, LINKS_HEADER1), file=sys.stderr)
                return None
            hasType = fileLine.startswith(LINKS_HEADER1)

            for fileLine in linksFile:
                if len(fileLine) == 0:
                    continue
                lineElems = fileLine.split(',')
                linkID = int(lineElems[0])
                if hasType:
                    if int(lineElems[1] != 1):
                        continue
                    source = int(lineElems[2])
                    dest = int(lineElems[3])
                else:
                    source = int(lineElems[1])
                    dest = int(lineElems[2])
                link = nmc_mm_lib.graph.GraphLink(linkID, graphLib.nodeMap[source], graphLib.nodeMap[dest])
                graphLib.addLink(link)
            
        # There we are.
        return graphLib
