"""
path_refine.py takes a path match file as input and reevaluates areas at restarts.
Outputs another path match CSV
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
from nmc_mm_lib import gtfs, path_engine, compat
import sys, transit_gtfs

def syntax():
    """
    Print usage information
    """
    print("path_refine takes a path match file as input and reevaluates areas at restarts. Outputs another path match CSV")
    print("Usage:")
    print("  python path_refine.py dbServer network user password shapePath pathMatchFile [-r filterRouteFile]")
    sys.exit(0)

def filterRoutes(gtfsNodes, shapePath, gtfsShapes, routeRestrictFilename, inclusiveFlag = False):
    """
    filterRoutes removes all GTFS nodes whose routes are not listed in the given file.
    @type gtfsNodes: dict<int, list<path_engine.PathEnd>>
    @type shapePath: str
    @type gtfsShapes: gtfs.ShapesEntry
    @type routeRestrictFilename: str
    @rtype dict<int, list<path_engine.PathEnd>>
    """
    # TODO: This type of filtering capability would probably best go in another centralized
    # location where multiple methods can all act upon the filter list and save processing time.

    if (routeRestrictFilename is None) or (len(routeRestrictFilename) == 0):
        return gtfsNodes

    print("INFO: Read GTFS routesfile...", file = sys.stderr)
    gtfsRoutes = gtfs.fillRoutes(shapePath)
    "@type gtfsRoutes: dict<int, gtfs.RoutesEntry>"
    
    print("INFO: Read GTFS tripsfile...", file = sys.stderr)
    (gtfsTrips, unusedTripIDs) = gtfs.fillTrips(shapePath, gtfsShapes, gtfsRoutes)
    "@type gtfsTrips: dict<int, gtfs.TripsEntry>"
    
    routeIDSet = set()
    "@type routeIDSet: set<int>"
    with open(routeRestrictFilename, 'r') as inFile:
        # Go through the lines of the file:
        for fileLine in inFile:
            if len(fileLine) > 0:                
                routeIDSet.add(int(fileLine))
    
    shapeIDSet = set()
    "@type shapeIDSet: set<int>"

    # Resolve shapeIDs from the routeIDs by getting all of the tripIDs that use the routes.
    if inclusiveFlag:
        matchingTrips = [gtfsTrip for gtfsTrip in gtfsTrips.values() if gtfsTrip.route.routeID not in routeIDSet]
    else:
        matchingTrips = [gtfsTrip for gtfsTrip in gtfsTrips.values() if gtfsTrip.route.routeID in routeIDSet]
    for gtfsTrip in matchingTrips:
        "@type gtfsTrip: gtfs.TripsEntry"
        shapeIDSet.add(gtfsTrip.shapeEntries[0].shapeID)
        
    # Report shapes:
    shapeStr = ""
    shapeIDReport = list(shapeIDSet)
    shapeIDReport.sort()
    ctr = 0
    for shapeID in shapeIDReport:
        shapeStr += str(shapeID)
        ctr += 1
        if ctr < len(shapeIDReport):
            shapeStr += ", "
    print("INFO: Route restrictions select Shape IDs %s." % shapeStr, file = sys.stderr)

    # Create a new GTFS list with the shapeIDs that had been selected:
    ret = {}
    "@type ret: dict<int, list<path_engine.PathEnd>>"
    shapeIDs = list(shapeIDSet)
    shapeIDs.sort()
    for shapeID in shapeIDs:
        if shapeID not in gtfsNodes:
            print("WARNING: Restriction shape ID %s does not exist in the GTFS set." % str(shapeID), file = sys.stderr)
        else:
            ret[shapeID] = gtfsNodes[shapeID]
    return ret


def pathsRefine(gtfsNodes, vistaGraph):
    # Default parameters, with explanations and cross-references to Perrine et al., 2015:
    termRefactorRadius = 3000   # Radius (ft) to invalidate found points at either end of a restart.
    pointSearchRadius = 1600    # "k": Radius (ft) to search from GTFS point to perpendicular VISTA links
    pointSearchPrimary = 1600   # "k_p": Radius (ft) to search from GTFS point to new VISTA links    
    pointSearchSecondary = 200  # "k_s": Radius (ft) to search from VISTA perpendicular point to previous point
    limitLinearDist = 6200      # Path distance (ft) to allow new proposed paths from one point to another
    limitDirectDist = 6200      # Radius (ft) to allow new proposed paths from one point to another
    limitDirectDistRev = 500    # Radius (ft) to allow backtracking on an existing link (e.g. parking lot)
    distanceFactor = 1.0        # "f_d": Cost multiplier for Linear path distance
    driftFactor = 1.5           # "f_r": Cost multiplier for distance from GTFS point to its VISTA link
    nonPerpPenalty = 1.5        # "f_p": Penalty multiplier for GTFS points that aren't perpendicular to VISTA links
    limitClosestPoints = 25     # "q_p": Number of close-proximity points that are considered for each GTFS point 
    limitSimultaneousPaths = 25 # "q_e": Number of proposed paths to maintain during pathfinding stage

    maxHops = 8                 # Maximum number of VISTA links to pursue in a path-finding operation
        
    # Initialize the path-finder:
    pathFinder = path_engine.PathEngine(pointSearchRadius, pointSearchPrimary, pointSearchSecondary, limitLinearDist,
                            limitDirectDist, limitDirectDistRev, distanceFactor, driftFactor, nonPerpPenalty, limitClosestPoints,
                            limitSimultaneousPaths)
    pathFinder.setRefineParams(termRefactorRadius)
    pathFinder.maxHops = maxHops
    
    # Begin iteration through each shape:
    shapeIDs = compat.listkeys(gtfsNodes)
    "@type shapeIDs: list<int>"
    shapeIDs.sort()
    gtfsNodesResults = {}
    "@type gtfsNodesResults: dict<int, list<path_engine.PathEnd>>"
    
    for shapeID in shapeIDs:
        "@type shapeID: int"
        
        print("INFO: -- Shape ID %s --" % str(shapeID), file = sys.stderr)
        
        # Find the path for the given shape:
        gtfsNodesRevised = pathFinder.refinePath(gtfsNodes[shapeID], vistaGraph) 
    
        # File this away as a result for later output:
        gtfsNodesResults[shapeID] = gtfsNodesRevised
    return gtfsNodesResults

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
    routeRestrictFilename = None
    if len(argv) > 6:
        i = 7
        while i < len(argv):
            if argv[i] == "-r" and i < len(argv) - 1:
                routeRestrictFilename = argv[i + 1]
                i += 1
            i += 1
    
    # Restore the stuff that was built with path_match:
    (vistaGraph, gtfsShapes, gtfsNodes, unusedShapeIDs) = transit_gtfs.restorePathMatch(dbServer, networkName,
        userName, password, shapePath, pathMatchFilename)
    # TODO: We don't do anything with unusedShapeIDs right now.
    
    # Filter down the routes that we're interested in:
    if routeRestrictFilename is not None:
        gtfsNodes = filterRoutes(gtfsNodes, shapePath, gtfsShapes, routeRestrictFilename)

    print("INFO: Refining paths.", file = sys.stderr)
    gtfsNodesResults = pathsRefine(gtfsNodes, vistaGraph)
    "@type gtfsNodesResults: dict<int, list<path_engine.PathEnd>>"
    
    print("INFO: -- Final --", file = sys.stderr)
    print("INFO: Print output...", file = sys.stderr)
    path_engine.dumpStandardHeader()

    shapeIDs = compat.listkeys(gtfsNodesResults)
    "@type shapeIDs: list<int>"
    shapeIDs.sort()
    for shapeID in shapeIDs:
        "@type shapeID: int"
        path_engine.dumpStandardInfo(gtfsNodesResults[shapeID])
        
    print("INFO: Done.", file = sys.stderr)

# Boostrap:
if __name__ == '__main__':
    main(sys.argv)
