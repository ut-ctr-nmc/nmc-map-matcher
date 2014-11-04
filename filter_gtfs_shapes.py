"""
filter_gtfs_shapes.py filters the shapefile contents to remove duplicate
    shapes and unwanted routes. It outputs to the console a new, reduced shapefile.
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
from nmc_mm_lib import gtfs, graph
import path_refine, difflib, sys

SEQUENCE_CUTOFF = 0.6

def syntax():
    """
    Print usage information
    """
    print("filter_gtfs_shapes.py filters the shapefile contents to remove duplicate shapes")
    print("and unwanted routes. It outputs to the console a new, reduced shapefile.")
    print("Usage:")
    print("  python filter_gtfs_shapes.py shapePath  [-x excludeRouteFile]")
    sys.exit(0)

def filterSimilarity(gtfsShapes):
    """
    Compares among all entries and figures out which ones are enough of duplicates. Keeps the
    longer lists.
    """
    shapeIDs = gtfsShapes.keys()
    
    excludedIDs = set()
    "@type excluded: set<int>"
    
    for origIndex in range(len(shapeIDs)):
        "@type origIndex: int"
        if shapeIDs[origIndex] not in excludedIDs:
            for targetIndex in range(len(shapeIDs)):
                "@type targetIndex: int"
                if (origIndex != targetIndex) and (shapeIDs[targetIndex] not in excludedIDs) \
                        and (len(gtfsShapes[shapeIDs[origIndex]]) >= len(gtfsShapes[shapeIDs[targetIndex]])):
                    s = difflib.SequenceMatcher()
                    s.set_seqs([(shapeEntry.lat, shapeEntry.lng) for shapeEntry in gtfsShapes[shapeIDs[origIndex]]],
                               [(shapeEntry.lat, shapeEntry.lng) for shapeEntry in gtfsShapes[shapeIDs[targetIndex]]])
                    if s.ratio() > SEQUENCE_CUTOFF:
                        excludedIDs.add(shapeIDs[targetIndex])
                        print("INFO: Shape ID %d is kept, where Shape ID %d is a duplicate (%.3g)" \
                              % (shapeIDs[origIndex], shapeIDs[targetIndex], s.ratio()), file = sys.stderr)
    ret = dict(gtfsShapes)
    for shapeID in excludedIDs:
        del ret[shapeID]
    return ret

def main(argv):
    # Initialize from command-line parameters:
    if (len(argv) < 1) or (argv[1].lower == "-h") or (argv[1].lower == "--help"):
        syntax()
    shapePath = argv[1]
    routeRestrictFilename = None
    if len(argv) > 1:
        i = 2
        while i < len(argv):
            if argv[i] == "-x" and i < len(argv) - 1:
                routeRestrictFilename = argv[i + 1]
                i += 1
            i += 1

    # Create a fake GPS coordinate:
    graph = graph.GraphLib(0, 0)
    
    # Read in the shapefile information:
    print("INFO: Read GTFS shapefile...", file = sys.stderr)
    gtfsShapes = gtfs.fillShapes(shapePath, graph.GPS)

    # Filter shapes according to exclusion file:
    if routeRestrictFilename is not None:
        gtfsShapes = path_refine.filterRoutes(gtfsShapes, shapePath, gtfsShapes, routeRestrictFilename, True)

    # Similarity search:
    gtfsShapes = filterSimilarity(gtfsShapes)
    
    # Extract useful information:
    print("INFO: Print output...", file = sys.stderr)
    shapeIDs = gtfsShapes.keys()
    "@type shapeIDs: list<int>"
    shapeIDs.sort()
    
    print("shape_id,shape_pt_lat,shape_pt_lon,shape_pt_sequence,shape_dist_traveled")
    for shapeID in shapeIDs:
        "@type shapeID: int"
        for shapeEntry in gtfsShapes[shapeID]:
            "@type shapeEntry: gtfs.ShapesEntry"
            print("%d,%f,%f,%d," % (shapeEntry.shapeID, shapeEntry.lat, shapeEntry.lng, shapeEntry.shapeSeq))
        
# Boostrap:
if __name__ == '__main__':
    main(sys.argv)
