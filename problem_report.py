"""
problem_report.py outputs GPS information for GTFS shapefiles reports potential
    problems with VISTA path matching
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

PERP_DIST = 300.0
"@var PERP_DIST: Maximum allowed distance in ft from VISTA line for perpendicular match" 
NONPERP_DIST = 150.0
"@var NONPERP_DIST: Maximum allowed distance in ft from VISTA line for non-perpendicular match" 

def problemReport(gtfsNodes, vistaGraph, byTripFlag=False, outFile=sys.stdout):
    """
    Takes a GTFS node set and outputs a CSV format of GPS points where there are indications of problems.
    @type gtfsNodes: dict<?, path_engine.PathEnd>
    @type vistaGraph: graph.GraphLib  
    """
    strStart = "shapeID," if not byTripFlag else "tripID,"
    print(strStart + "shapeSeq,linkID,linkDist,problemCode,gtfsCoords,vistaCoords", file = outFile)

    shapeIDs = gtfsNodes.keys()
    shapeIDs.sort()
    for shapeID in shapeIDs:
        gtfsNodeList = gtfsNodes[shapeID]
        "@type gtfsNodeList: list<path_engine.PathEnd>"
        for gtfsNode in gtfsNodeList:
            "@type gtfsNode: path_engine.PathEnd"
            (vistaLat, vistaLng) = vistaGraph.gps.feet2gps(gtfsNode.pointOnLink.pointX, gtfsNode.pointOnLink.pointY) 
            
            # Determine whether we have a problem to report:
            problemCode = 0
            if gtfsNode.restart:
                problemCode = 1
            elif not gtfsNode.pointOnLink.nonPerpPenalty and gtfsNode.pointOnLink.refDist > PERP_DIST:
                problemCode = 2
            elif gtfsNode.pointOnLink.nonPerpPenalty and gtfsNode.pointOnLink.refDist > NONPERP_DIST:
                problemCode = 3
            
            outStr = "%s,%d,%d,%g,%d,%s,%s" % (str(gtfsNode.shapeEntry.shapeID), gtfsNode.shapeEntry.shapeSeq,
                            gtfsNode.pointOnLink.link.id if gtfsNode.pointOnLink.link is not None else -1,
                            gtfsNode.pointOnLink.dist, problemCode, str(gtfsNode.shapeEntry.lat)
                            + " " + str(gtfsNode.shapeEntry.lng), str(vistaLat) + " " + str(vistaLng))
            print(outStr, file = outFile)

def syntax():
    """
    Print usage information
    """
    print("problem_report outputs GPS information for GTFS shapefiles reports potential problems with VISTA path matching.")
    print("Problem codes: 0: OK; 1: Path restarted; 2: For perpendicular; 3: For nonperpendicular")
    print("Usage:")
    print("  python problem_report.py dbServer network user password shapePath pathMatchFile")
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
    (vistaGraph, gtfsShapes, gtfsNodes, unusedShapeIDs) = transit_gtfs.restorePathMatch(dbServer, networkName,
        userName, password, shapePath, pathMatchFilename)
    print("INFO: Output CSV...", file = sys.stderr)
    problemReport(gtfsNodes, vistaGraph)
    print("INFO: Done.", file = sys.stderr)
    
# Boostrap:
if __name__ == '__main__':
    main(sys.argv)
