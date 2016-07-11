"""
problem_report_osmxml.py reports potential problems with OpenStreetMap path matching
    
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
import path_match_osmxml, problem_report
from nmc_mm_lib import gtfs, path_engine, compat
import sys, argparse

PERP_DIST = 300.0
"@var PERP_DIST: Maximum allowed distance in ft from OSM line for perpendicular match" 
NONPERP_DIST = 150.0
"@var NONPERP_DIST: Maximum allowed distance in ft from OSM line for non-perpendicular match" 

def restorePathMatch(osmXML, shapePath, pathMatchFilename, trackpoints=False):
    """
    This is like transit_gtfs.restorePathMatch(), but geared for OSM XML underlying topology.
    As before, this is a temporary fix until the input structure is completed in the next version.
    """
    
    # Read in the OSM XML:
    print ("INFO: Importing OSM XML from '%s'..." % osmXML, file=sys.stderr)
    underlyingTopo = path_match_osmxml.fillGraph(osmXML)
    
    
    
    
    # DEBUG
    import json
    outp = {}
    nodes = {}
    underlyingTopo._flatten()
    for node in underlyingTopo.nodeMap.values():
        entry = {}
        entry["ol"] = compat.listkeys(node.outgoingLinkMap)
        nodes[node.id] = entry
    outp["nodes"] = nodes
    links = {}
    for link in underlyingTopo.linkMap.values():
        entry = {}
        entry["start"] = link.origNode.id
        entry["end"] = link.destNode.id
        entry["vc"] = len(link.vertices)
        entry["sn"] = link.streetName
        vertices = []
        for vertex in link.vertices:
            vertices.append(vertex.id)
        entry["v"] = vertices
        links[link.id] = entry
    outp["links"] = links
    print(json.dumps(outp, sort_keys=True, indent=4, separators=(',', ': ')))
    raise Exception


    
    
    if not trackpoints:
        # Read in the shapefile information:
        print("INFO: Read GTFS shapefile...", file=sys.stderr)
        gtfsShapes = gtfs.fillShapes(shapePath, underlyingTopo.gps)
    else:
        # TODO: Complete
        pass

    # Read the path-match file:
    print("INFO: Read the path-match file '%s'..." % pathMatchFilename, file=sys.stderr)
    with open(pathMatchFilename, 'r') as inFile:
        gtfsNodes = path_engine.readStandardDump(underlyingTopo, gtfsShapes, inFile)
        "@type gtfsNodes: dict<int, list<path_engine.PathEnd>>"

    # Filter out the unused shapes:
    unusedShapeIDs = set()
    for shapeID in compat.listkeys(gtfsShapes):
        if shapeID not in gtfsNodes:
            del gtfsShapes[shapeID]
            unusedShapeIDs.add(shapeID)

    return underlyingTopo, gtfsShapes, gtfsNodes, unusedShapeIDs

def main(argv):
    # Initialize from command-line parameters:
    parser = argparse.ArgumentParser(description="problem_report outputs GPS information for GTFS shapefiles, reporting " +
        "potential problems with VISTA path matching. Problem codes: 0: OK; 1: Path restarted; 2: For perpendicular; " +
        "3: For nonperpendicular; 4: Intermediate link")
    parser.add_argument("osmXML", help="Filename for the OSM XML file")
    parser.add_argument("shapePath", help="Path to GTFS files or trackpoint file")
    parser.add_argument("pathMatchFile", help="Path match file that was previously output by path_match.py or others")
    parser.add_argument("-L", "--interlinks", action="store_true", help="Outputs intermediate links and positions")
    parser.add_argument("-T", "--trackpoints", action="store_true", help="Read trackpoint CSV file rather than a GTFS shapefile")
    args = parser.parse_args()
        
    # Restore the stuff that was built with path_match:
    underlyingTopo, gtfsShapes, gtfsNodes, unusedShapeIDs = restorePathMatch(args.osmXML, args.shapePath, args.pathMatchFile, trackpoints=args.trackpoints)
    print("INFO: Output CSV...", file=sys.stderr)
    problem_report.problemReport(gtfsShapes, gtfsNodes, underlyingTopo, showLinks=args.interlinks)
    print("INFO: Done.", file = sys.stderr)
    
# Boostrap:
if __name__ == '__main__':
    main(sys.argv)
