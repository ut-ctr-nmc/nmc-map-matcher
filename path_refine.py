"""
path_refine.py takes a path match file as input and reevaluates areas at hints
    and restarts. Outputs another path match CSV
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
from nmc_mm_lib import gtfs, path_engine
import sys, operator, transit_gtfs

def syntax():
    """
    Print usage information
    """
    print("path_refine takes a path match file as input and reevaluates areas at hints and")
    print("restarts. Outputs another path match CSV")
    print("Usage:")
    print("  python path_refine.py dbServer network user password shapePath pathMatchFile [-h hintFile] [-r filterRouteFile]")
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

def fillHints(filename, shapePath, gtfsShapes, GPS, unusedShapeIDs):
    """
    fillHints retrieves hints that guide where paths should be drawn. 
    @type filePath: str
    @type shapePath: str
    @type gtfsShapes: gtfs.ShapesEntry
    @type GPS: GPS.GPS
    @return A map of stop_id to a hint_entry.
    @rtype dict<int, path_engine.ShapesEntry>
    """
    ret = {}
    "@type ret: dict<int, list<path_engine.ShapesEntry>>"
    
    if (filename is None) or (len(filename) == 0):
        return ret
    
    # Retrieve trips file information that will get us the routes-to-shape mapping that we want. 
    print("INFO: Read GTFS routesfile...", file = sys.stderr)
    gtfsRoutes = gtfs.fillRoutes(shapePath)
    "@type gtfsRoutes: dict<int, gtfs.RoutesEntry>"
    
    print("INFO: Read GTFS tripsfile...", file = sys.stderr)
    (gtfsTrips, unusedTripIDs) = gtfs.fillTrips(shapePath, gtfsShapes, gtfsRoutes, unusedShapeIDs)
    "@type gtfsTrips: dict<int, gtfs.TripsEntry>"
        
    with open(filename, 'r') as inFile:
        # Sanity check:
        fileLine = inFile.readline()
        if not fileLine.startswith("route_id,hint_seq,lat,lon"):
            print("ERROR: The hints file '%s' doesn't have the expected header." % filename, file = sys.stderr)
            return None
        
        # Go through the lines of the file:
        for fileLine in inFile:
            if len(fileLine) > 0:
                lineElems = fileLine.split(',')
                
                routeID = int(lineElems[0])
                hintSeq = int(lineElems[1])
                gpsLat = float(lineElems[2])
                gpsLng = float(lineElems[3])
                
                # Resolve shapeIDs from the routeIDs by getting all of the tripIDs that use the route.
                matchingTrips = [gtfsTrip for gtfsTrip in gtfsTrips.values() if gtfsTrip.route.routeID == routeID]
                shapeIDSet = set()
                for gtfsTrip in matchingTrips:
                    "@type gtfsTrip: gtfs.TripsEntry"
                    shapeID = gtfsTrip.shapeEntries[0].shapeID
                    if shapeID not in shapeIDSet: # Only write in a hint once per shape.
                        newEntry = gtfs.ShapesEntry(shapeID, hintSeq, gpsLat, gpsLng, True)
                        (newEntry.pointX, newEntry.pointY) = GPS.gps2feet(gpsLat, gpsLng)
                
                        if shapeID not in ret:
                            ret[shapeID] = list()
                        ret[shapeID].append(newEntry)
                        shapeIDSet.add(shapeID)

    # Sort the hints so that they will be intercepted in the correct order:                    
    for shapeID in ret:
        ret[shapeID].sort(key = operator.attrgetter('shapeSeq'))
    
    # Return the hints file contents:
    return ret

def pathsRefine(gtfsNodes, hintEntries, vistaGraph):
    # Default parameters, with explanations and cross-references to Perrine et al., 2015:
    hintRefactorRadius = 1000   # Radius (ft) to invalidate surrounding found points.
    termRefactorRadius = 3000   # Radius (ft) to invalidate found points at either end of a restart.
    pointSearchRadius = 1600    # "k": Radius (ft) to search from GTFS point to perpendicular VISTA links
    pointSearchPrimary = 1600   # "k_p": Radius (ft) to search from GTFS point to new VISTA links    
    pointSearchSecondary = 200  # "k_s": Radius (ft) to search from VISTA perpendicular point to previous point
    limitLinearDist = 6200      # Path distance (ft) to allow new proposed paths from one point to another
    limitDirectDist = 6200      # Radius (ft) to allow new proposed paths from one point to another
    limitDirectDistRev = 500    # Radius (ft) to allow backtracking on an existing link (e.g. parking lot)
    distanceFactor = 1.0        # "f_d": Cost multiplier for Linear path distance
    driftFactor = 1.5           # "f_r": Cost multiplier for distance from GTFS point to its VISTA link
    nonPerpPenalty = 1.5        # "f_n": Penalty multiplier for GTFS points that aren't perpendicular to VISTA links
    limitClosestPoints = 25     # "q_p": Number of close-proximity points that are considered for each GTFS point 
    limitSimultaneousPaths = 25 # "q_e": Number of proposed paths to maintain during pathfinding stage

    maxHops = 8                 # Maximum number of VISTA links to pursue in a path-finding operation
    limitHintClosest = 4        # Number of hint closest points and closest previous track points
        
    # Initialize the path-finder:
    pathFinder = path_engine.PathEngine(pointSearchRadius, pointSearchPrimary, pointSearchSecondary, limitLinearDist,
                            limitDirectDist, limitDirectDistRev, distanceFactor, driftFactor, nonPerpPenalty, limitClosestPoints,
                            limitSimultaneousPaths)
    pathFinder.setRefineParams(hintRefactorRadius, termRefactorRadius)
    pathFinder.maxHops = maxHops
    pathFinder.limitHintClosest = limitHintClosest
    
    # Begin iteration through each shape:
    shapeIDs = gtfsNodes.keys()
    "@type shapeIDs: list<int>"
    shapeIDs.sort()
    gtfsNodesResults = {}
    "@type gtfsNodesResults: dict<int, list<path_engine.PathEnd>>"
    
    for shapeID in shapeIDs:
        "@type shapeID: int"
        
        print("INFO: -- Shape ID %s --" % str(shapeID), file = sys.stderr)
        
        # Find the path for the given shape:
        gtfsNodesRevised = pathFinder.refinePath(gtfsNodes[shapeID], vistaGraph, 
            hintEntries[shapeID] if shapeID in hintEntries else list()) 
    
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
    hintFilename = None
    routeRestrictFilename = None
    if len(argv) > 6:
        i = 7
        while i < len(argv):
            if argv[i] == "-h" and i < len(argv) - 1:
                hintFilename = argv[i + 1]
                i += 1
            elif argv[i] == "-r" and i < len(argv) - 1:
                routeRestrictFilename = argv[i + 1]
                i += 1
            i += 1
    
    # Restore the stuff that was built with path_match:
    (vistaGraph, gtfsShapes, gtfsNodes, unusedShapeIDs) = transit_gtfs.restorePathMatch(dbServer, networkName,
        userName, password, shapePath, pathMatchFilename)
    # TODO: We don't do anything with unusedShapeIDs right now.
    
    # Restore the hint file if it is specified:
    if hintFilename is not None:
        print("INFO: Read hint file...", file = sys.stderr)
    else:
        print("INFO: No hint file was specified.", file = sys.stderr)
    hintEntries = fillHints(hintFilename, shapePath, gtfsShapes, vistaGraph.GPS, unusedShapeIDs)
    "@type hintEntries: dict<int, path_engine.ShapesEntry>"

    # Filter down the routes that we're interested in:
    if routeRestrictFilename is not None:
        gtfsNodes = filterRoutes(gtfsNodes, shapePath, gtfsShapes, routeRestrictFilename)

    print("INFO: Refining paths.", file = sys.stderr)
    gtfsNodesResults = pathsRefine(gtfsNodes, hintEntries, vistaGraph)
    "@type gtfsNodesResults: dict<int, list<path_engine.PathEnd>>"
    
    print("INFO: -- Final --", file = sys.stderr)
    print("INFO: Print output...", file = sys.stderr)
    path_engine.dumpStandardHeader()

    shapeIDs = gtfsNodesResults.keys()
    "@type shapeIDs: list<int>"
    shapeIDs.sort()
    for shapeID in shapeIDs:
        "@type shapeID: int"
        path_engine.dumpStandardInfo(gtfsNodesResults[shapeID])
        
    print("INFO: Done.", file = sys.stderr)

# Boostrap:
if __name__ == '__main__':
    main(sys.argv)
