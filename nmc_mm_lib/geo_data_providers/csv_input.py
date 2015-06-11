"""
csv_input.py is an underlying node-link data representation of a
    roadway network coded within specifically-defined CSV files.
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
import nmc_mm_lib.graph, nmc_mm_lib.geo_data_provider, sys

NODES_HEADER1 = "id,type,x,y"
NODES_HEADER2 = "id,x,y"
LINKS_HEADER1 = "id,type,source,destination"
LINKS_HEADER2 = "id,source,destination"

class CSVInput(nmc_mm_lib.geo_data_provider.GeoDataProvider):
    """
    @ivar nodesFilename: str representing the nodes CSV filename
    @ivar linksFilename: str representing the links CSV filename
    """
    
    def __init__(self):
        "@type self.connection: psycopg2.connection"
        self.nodesFilename = ""
        self.linksFilename = ""
        
        super(nmc_mm_lib.geo_data_provider.GeoDataProvider, self).__init__()
    
    def provideCmdLineOpts(self, options):
        """
        Stuffs new command-line options into the given argparse object.
        @type options: argparse
        """
        options.add_argument("-fn", nargs=1, type=str, metavar="NODESFILENAME", dest="nodesFilename",
            help="If using CSV files for input, specifies the nodes CSV file")
        options.add_argument("-fl", nargs=1, type=str, metavar="LINKSFILENAME", dest="linksFilename",
            help="If using CSV files for input, specifies the links CSV file")

    def checkCmdLineOpts(self, args):
        """
        Checks the given command-line options and determines if there are sufficient options
        for successfully loading in data for this class.
        @type args: argparse.Namespace
        @return True if sufficient command-line options are provided for this class
        """
        return args.nodesFilename is not None or args.linksFilename is not None
    
    def readCmdLineOpts(self, args):
        """
        Reads in all relevant command-line options and emits an error if there are insufficient
        options provided.
        @type args: argparse.Namespace
        @return False if there are insufficient options provided, or True if successful.
        @rtype boolean
        """
        self.nodesFilename = args.nodesFilename
        if self.nodesFilename is None:
            print("ERROR: When using CSV files for input, a nodes filename -fn must be specified.", file=sys.stderr)
            return False
        self.linksFilename = args.linksFilename
        if self.linksFilename is None:
            print("ERROR: When using CSV files for input, a links filename -fl must be specified.", file=sys.stderr) 
            return False
        return True
    
    def readData(self):
        """
        Using information from the command-line options, reads in data from the geographic data
        source that's handled by this class.
        @return A graph representing the geographic data
        @rtype graph.GraphLib
        """
        print("INFO: Reading geographic input from CSV files '%s' and '%s'..." % (self.nodesFilename, self.linksFilename), file=sys.stderr)
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
