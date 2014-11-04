"""
dump_gps.py outputs GPS information for GTFS shapefiles and VISTA points.
@author: Kenneth Perrine
@contact: kperrine@utexas.edu
@organization: Network Modeling Center, Center for Transportation Research,
    Cockrell School of Engineering, The University of Texas at Austin 
@version: 1.0

@copyright: (C) 2014, The University of Texas at Austin
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
import transit_gtfs, sys

def dumpGPS(gtfsNodes, vistaGraph, outFile = sys.stdout):
    """
    Takes a GTFS node set and outputs a CSV format of GPS points and other information.
    @type gtfsNodes: list<path_engine.PathEnd>
    @type vistaGraph: graph.GraphLib  
    """
    print("shapeID,shapeSeq,linkID,linkDist,gtfsLat,gtfsLng,vistaLat,vistaLng,vistaNodeLat,vistaNodeLng", file = outFile)

    shapeIDs = gtfsNodes.keys()
    shapeIDs.sort()
    for shapeID in shapeIDs:
        gtfsNodeList = gtfsNodes[shapeID]
        "@type gtfsNodeList: list<path_engine.PathEnd>"
        for gtfsNode in gtfsNodeList:
            "@type gtfsNode: path_engine.PathEnd"
            (vistaLat, vistaLng) = vistaGraph.GPS.feet2gps(gtfsNode.pointOnLink.pointX, gtfsNode.pointOnLink.pointY) 
            
            outStr = "%d,%d,%d,%g,%g,%g,%g,%g,%g,%g" % (gtfsNode.shapeEntry.shapeID, gtfsNode.shapeEntry.shapeSeq,
                            gtfsNode.pointOnLink.link.id, gtfsNode.pointOnLink.dist, gtfsNode.shapeEntry.lat, gtfsNode.shapeEntry.lng,
                            vistaLat, vistaLng, gtfsNode.pointOnLink.link.origNode.gpsLat, gtfsNode.pointOnLink.link.origNode.gpsLng)
            print(outStr, file = outFile)

def syntax():
    """
    Print usage information
    """
    print("dump_gps outputs GPS information for GTFS shapefiles and VISTA points.")
    print("Usage:")
    print("  python dump_gps.py dbServer network user password shapePath pathMatchFile")
    sys.exit(0)

def main(argv):
    # Initialize from command-line parameters:
    if len(argv) < 7:
        syntax()
    dbServer = argv[1]
    networkName = argv[2]
    userName = argv[3]
    password = argv[4]
    shapePath = argv[5]
    pathMatchFilename = argv[6]
        
    # Restore the stuff that was built with path_match:
    (vistaGraph, gtfsShapes, gtfsNodes, unusedShapeIDs) = transit_gtfs.restorePathMatch(dbServer, networkName, userName, password,
                                                                        shapePath, pathMatchFilename)
    print("INFO: Output CSV...", file = sys.stderr)
    dumpGPS(gtfsNodes, vistaGraph)
    print("INFO: Done.", file = sys.stderr)
    
# Boostrap:
if __name__ == '__main__':
    main(sys.argv)
