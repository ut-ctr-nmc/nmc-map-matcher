"""
path_match.py resolves a GTFS shapefile to a VISTA network series of links and
    outputs a CSV format of data
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
from nmc_mm_lib import gtfs, vista_network, path_engine, compat
import sys

def syntax():
    """
    Print usage information
    """
    print("path_match.py resolves a GTFS shapefile to a VISTA network series of links and")
    print("outputs a CSV format of data.")
    print("Usage:")
    print("  python path_match.py dbServer network user password shapePath")
    sys.exit(0)

def pathMatch(dbServer, networkName, userName, password, shapePath, limitMap=None, pickleOutTopo=None):
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
    
    # Get the database connected:
    print("INFO: Connect to database...", file = sys.stderr)
    database = vista_network.connect(dbServer, userName, password, networkName)
    
    # Read in the topology from the VISTA database:
    print("INFO: Read topology from database...", file = sys.stderr)
    vistaGraph = vista_network.fillGraph(database)
    
    if pickleOutTopo:
        print("INFO: Writing out underlying topology to %s" % pickleOutTopo, file=sys.stderr)
        #sys.setrecursionlimit(3800)
        with open(pickleOutTopo, "wb") as pickleFile:
            vistaGraph.serialize(pickleFile)
        pickleFile.close()
        del pickleFile
    
    # Read in the shapefile information:
    print("INFO: Read GTFS shapefile...", file = sys.stderr)
    gtfsShapes = gtfs.fillShapes(shapePath, vistaGraph.gps)
    
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

def main(argv):
    # Initialize from command-line parameters:
    if len(argv) < 6:
        syntax()
    dbServer = argv[1]
    networkName = argv[2]
    userName = argv[3]
    password = argv[4]
    shapePath = argv[5]
    pickleOutTopo = None
    
    i = 6
    while i < len(argv):
        if argv[i] == "--ptout" and i < len(argv) - 1:
            pickleOutTopo = argv[i + 1]
            i += 1
        i += 1
    
    gtfsNodesResults = pathMatch(dbServer, networkName, userName, password, shapePath, pickleOutTopo=pickleOutTopo)
    
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
