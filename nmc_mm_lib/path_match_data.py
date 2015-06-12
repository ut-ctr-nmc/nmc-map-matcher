"""
path_match_data.py handles the inputting and outputting of path match data.
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
import gtfs, path_engine, geo_data_provider, linear, sys

HEADER = "shapeID,shapeSeq,shapeType,linkID,linkDist,totalDist,numLinksTrav,linksTrav"

def provideCmdLineOpts(options):
    """
    Stuffs new command-line options into the given argparse object.
    @type options: argparse
    """
    options.add_argument("-p", nargs=1, type=str, metavar="PATHMATCHFILE", dest="pathMatchFilename",
        help="Path match CSV file that is the output of path_match.py and others", required=True)

def restorePathMatch(args):
    """
    Performs loading operations for the underlying geometry, shape file, and path match file.
    geo_data_provider.loadGeoDataProviders() and provideCmdLineOpts() must have already been called.
    @type args: argparse.Namespace
    @rtype (graph.GraphLib, dict<int, list<ShapesEntry>>, dict<int, list<path_engine.PathEnd>>, set<?>) 
    """
    # Get the underlying geometry:
    print("INFO: Reading underlying topology...", file=sys.stderr)
    graph = geo_data_provider.readData()
        
    # Read in the shapefile information:
    print("INFO: Read GTFS shapefile...", file = sys.stderr)
    shapes = gtfs.fillShapes(shapePath, graph.gps)

    # Read the path-match file:
    print("INFO: Read the path-match file '%s'..." % args.pathMatchFilename, file=sys.stderr)
    with open(args.pathMatchFilename, 'r') as inFile:
        nodes = _readPathMatch(graph, shapes, inFile)
        "@type nodes: dict<int, list<path_engine.PathEnd>>"

    # Filter out the unused shapes:
    unusedShapeIDs = set()
    for shapeID in shapes.keys():
        if shapeID not in nodes:
            del shapes[shapeID]
            unusedShapeIDs.add(shapeID)

    return (graph, shapes, nodes, unusedShapeIDs)

def _readPathMatch(graph, shapes, inFile, shapeIDMaker = lambda x: int(x)):
    """
    readPathMatch reconstructs the tree entries that PathEngine had created from the
    path match file given by path inFile.
    @type graph: graph.GraphLib
    @type shapes: dict<int, list<gtfs.ShapesEntry>>
    @type inFile: file
    @type shapeIDMaker: function
    @return A dictionary of shapeID to a list of PathEnds
    @rtype dict<int, list<path_engine.PathEnd>>
    """
    ret = {}
    "@type ret: dict<int, list<path_engine.PathEnd>>"

    # Sanity check:
    fileLine = inFile.readline()
    if not fileLine.startswith(HEADER):
        print("ERROR: The path match file doesn't have the expected header.", file = sys.stderr)
        return None
        
    # Storage place for sequence numbers:
    shapeSeqs = {}
    "@type shapeSeqs: dict<int, int>"
    
    hintSeqs = {}
    "@type hintSeqs: dict<int, int>"
        
    # Go through the lines of the file:
    for fileLine in inFile:
        if len(fileLine) > 0:
            lineElems = fileLine.split(',')
            shapeID = shapeIDMaker(lineElems[0])
            shapeSeq = int(lineElems[1])
            hintFlag = int(lineElems[2]) != 0
            linkID = int(lineElems[3])
            linkDist = float(lineElems[4])
            distTotal = float(lineElems[5])
            linksTravCount = int(lineElems[6])
            if linksTravCount >= 0:
                linksTrav = linksTravCount * [None]
            else:
                # If linksTravCount is -1, then that signifies that we are restarting.  Deal with it later.
                linksTrav = []
            "@type linksTrav: list<graph.GraphNode>"
                        
            # Get the variable-length link list that happens at the end:
            contFlag = False
            for index in range(0, len(linksTrav)):
                linksTravID = int(lineElems[index + 7])
                if linksTravID not in graph.linkMap:
                    print("WARNING: The path match file refers to a nonexistent link ID %d." % linksTravID, file = sys.stderr)
                    contFlag = True
                    break
                linksTrav[index] = graph.linkMap[linksTravID]
            if contFlag:
                # This is run if the break above is run.
                continue

            # Resolve the shapeID list:
            if shapeID not in shapes:
                print("WARNING: The path match file refers to a nonexistent shape ID %s." % str(shapeID), file = sys.stderr)
                continue
            shapeElems = shapes[shapeID]
            "@type shapeElems: list<gtfs.ShapesEntry>"
            
            # Set up the shape index cache to reduce linear searching later on:
            if not hintFlag:
                if shapeID not in shapeSeqs:
                    shapeSeqs[shapeID] = -1
            else:
                if shapeID not in hintSeqs:
                    hintSeqs[shapeID] = -1
                
            # Resolve the link object:
            if linkID not in graph.linkMap:
                print("WARNING: The path match file refers to a nonexistent link ID %d." % linkID, file = sys.stderr)
                continue
            link = graph.linkMap[linkID]
            "@type link: graph.GraphLink"
            
            # Resolve the shape entry:
            shapeEntry = None
            "@type shapeEntry: gtfs.ShapesEntry"
            if not hintFlag:
                for index in range(shapeSeqs[shapeID] + 1, len(shapeElems)):
                    if shapeElems[index].shapeSeq == shapeSeq:
                        shapeEntry = shapeElems[index]
                        shapeSeqs[shapeID] = index
                        break
                if shapeEntry is None:
                    print("WARNING: No GTFS shape entry for shape ID: %s, seq: %d; check for out of order."
                          % (str(shapeID), shapeSeq), file = sys.stderr)
                    continue
            else:
                # Reconstruct the hint:
                # TO DO: Without the hint file, we won't get the exact same GPS coordinates as were used for the hint.
                if shapeSeq > hintSeqs[shapeID]:
                    hintSeqs[shapeID] = shapeSeq
                    pointX = link.origNode.coordX + (link.destNode.coordX - link.origNode.coordX) * linkDist / link.distance
                    pointY = link.origNode.coordY + (link.destNode.coordY - link.origNode.coordY) * linkDist / link.distance
                    (lat, lng) = graph.gps.feet2gps(pointX, pointY)
                    shapeEntry = gtfs.ShapesEntry(shapeID, shapeSeq, lat, lng, True)
                    (shapeEntry.pointX, shapeEntry.pointY) = (pointX, pointY)
                else:
                    print("WARNING: Hint entry may be out of order: Shape: %s, seq: %d.  Hint will not be used."
                          % (str(shapeID), shapeSeq), file = sys.stderr)
                    continue
            
            # Recalculate parameters needed for the tree node:
            (pointX, pointY) = graph.gps.gps2feet(shapeEntry.lat, shapeEntry.lng)
            (distRef, distLinear, perpFlag) = linear.pointDist(pointX, pointY, link.origNode.coordX, link.origNode.coordY,
                                                               link.destNode.coordX, link.destNode.coordY)
            pointOnLink = graph.PointOnLink(link, distLinear, not perpFlag, distRef)
            newEntry = path_engine.PathEnd(shapeEntry, pointOnLink)
            newEntry.totalCost = distTotal # TotalCost won't be available.
            newEntry.totalDist = distTotal
            if linksTravCount == -1:
                # We are restarting the path and don't have complete information up to this point.
                newEntry.restart = True
            
            # Reconstruct the node list:
            newEntry.routeInfo = linksTrav
            
            if shapeID not in ret:
                ret[shapeID] = []
            else:
                # Restore the previous tree entry linkage:
                newEntry.prevTreeNode = ret[shapeID][len(ret[shapeID]) - 1]
            ret[shapeID].append(newEntry)

    # Return the tree nodes:
    return ret

