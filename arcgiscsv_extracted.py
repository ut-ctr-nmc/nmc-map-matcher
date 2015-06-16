"""
arcgiscsv_extracted.py matches an ArcGIS CSV format GPS tracks to a VISTA network series
    of links and outputs a CSV format of data in the standard output format.
    Note the presence of hard-coded global parameters in pathMatch().
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
from datetime import datetime, timedelta
from nmc_mm_lib import gtfs, vista_network, path_engine
import operator, sys, csv

# A module that deals with reading CSV files extracted from an ArcGIS CSV export format.
DIVISOR = 1000L * long(1 << 48)

def syntax():
    """
    Print usage information
    """
    print("arcgiscsv_extracted.py matches an ArcGIS text format GPS track to a VISTA network series")
    print("of links and outputs a CSV format of data in the standard output format.")
    print("Usage:")
    print("  python gdb_extracted.py dbServer network user password arcgiscsvFile")
    sys.exit(0)

def fillFromFile(filename, GPS):
    """
    fillFromFile retrieves the ArcGIS information from the CSV dump and builds up a series of shapes
    varying by userID from it
    @type filename: str
    @type GPS: GPS.GPS
    @return A map of userID to a list of shape entries
    @rtype dict<str, list<gtfs.ShapesEntry>>
    """
    ret = {}
    "@type ret: dict<str, list<gtfs.ShapesEntry>>"
    with open(filename, 'r') as inFile:
        csvReader = csv.DictReader(inFile)

        firstRun = True
        for fileLine in csvReader:
            if firstRun:
                # Sanity check:
                if not ("RouteId" in fileLine and "Datafile" in fileLine and "Latitude" in fileLine
                        and "Longitude" in fileLine and "OBJECTID" in fileLine and "UtcDateTim" in fileLine):
                    print("ERROR: The ArcGIS CSV file %s doesn't have the expected header." % filename, file = sys.stderr)
                    return None
                firstRun = False
        
            routeID = int(fileLine["RouteId"].replace(",", ""))
            if routeID > 0:
                identifier = fileLine["Datafile"] + "(" + str(routeID) + ")"
                newEntry = gtfs.ShapesEntry(identifier, int(fileLine["OBJECTID"].replace(",", "")), float(fileLine["Latitude"]),
                                    float(fileLine["Longitude"]), False)
                (newEntry.pointX, newEntry.pointY) = GPS.gps2feet(newEntry.lat, newEntry.lng)
                newEntry.time = datetime.strptime(fileLine["UtcDateTim"], '%m/%d/%Y %I:%M:%S %p')
                
# Temporary set to CST: 
                timeDelta = timedelta(hours=-6)
                newEntry.time += timeDelta
                
                # Keep only the time:
                newEntry.time = datetime.strptime("%02d:%02d:%02d" % (newEntry.time.hour, newEntry.time.minute,
                                                                      newEntry.time.second), '%H:%M:%S')
                
                if newEntry.shapeID not in ret:
                    ret[newEntry.shapeID] = []
                ret[newEntry.shapeID].append(newEntry)

    # Ensure that the lists are sorted:
    for shapesEntries in ret.values():
        "@type shapesEntries: list<ShapesEntry>"
        shapesEntries.sort(key = operator.attrgetter('shapeSeq'))
                    
    # Return the shapes file contents:
    return ret

def pathMatch(dbServer, networkName, userName, password, filename, limitMap = None):
    # Default parameters, with explanations and cross-references to Perrine et al., 2015:
    pointSearchRadius = 1000    # "k": Radius (ft) to search from GTFS point to perpendicular VISTA links
    pointSearchPrimary = 350    # "k_p": Radius (ft) to search from GTFS point to new VISTA links    
    pointSearchSecondary = 200  # "k_s": Radius (ft) to search from VISTA perpendicular point to previous point
    limitLinearDist = 500      # Path distance (ft) to allow new proposed paths from one point to another
    limitDirectDist = 400      # Radius (ft) to allow new proposed paths from one point to another
    limitDirectDistRev = 300    # Radius (ft) to allow backtracking on an existing link (e.g. parking lot)
    distanceFactor = 1.0        # "f_d": Cost multiplier for Linear path distance
    driftFactor = 1.5           # "f_r": Cost multiplier for distance from GTFS point to its VISTA link
    nonPerpPenalty = 1.5        # "f_p": Penalty multiplier for GTFS points that aren't perpendicular to VISTA links
    limitClosestPoints = 10     # "q_p": Number of close-proximity points that are considered for each GTFS point 
    limitSimultaneousPaths = 8  # "q_e": Number of proposed paths to maintain during pathfinding stage
    
    maxHops = 8                # Maximum number of VISTA links to pursue in a path-finding operation
    
    # Get the database connected:
    print("INFO: Connect to database...", file = sys.stderr)
    database = vista_network.connect(dbServer, userName, password, networkName)
    
    # Read in the topology from the VISTA database:
    print("INFO: Read topology from database...", file = sys.stderr)
    vistaGraph = vista_network.fillGraph(database)
    
    # Read in the GPS track information:
    print("INFO: Read ArcGIS CSV GPS track...", file = sys.stderr)
    gpsTracks = fillFromFile(filename, vistaGraph.gps)
    
    # Initialize the path-finder:
    pathFinder = path_engine.PathEngine(pointSearchRadius, pointSearchPrimary, pointSearchSecondary, limitLinearDist,
                            limitDirectDist, limitDirectDistRev, distanceFactor, driftFactor, nonPerpPenalty, limitClosestPoints,
                            limitSimultaneousPaths)
    pathFinder.maxHops = maxHops
    
    # Begin iteration through each shape:
    datafileIDs = gpsTracks.keys()
    "@type datafileIDs: list<str>"
    datafileIDs.sort()
    nodesResults = {}
    "@type nodesResults: dict<str, list<path_engine.PathEnd>>"
    
    if limitMap is not None:
        for datafileID in limitMap:
            if datafileID not in datafileIDs:
                print("WARNING: Limit datafile ID %d is not found in the shape file." % datafileID, file = sys.stderr)
    
    for datafileID in datafileIDs:
        "@type datafileID: int"
        
        if limitMap is not None and datafileID not in limitMap:
            continue
        
        print("INFO: -- Datafile %s --" % datafileID, file = sys.stderr)
        
        # Find the path for the given shape:
        gtfsNodes = pathFinder.constructPath(gpsTracks[datafileID], vistaGraph)
    
        # File this away as a result for later output:
        nodesResults[datafileID] = gtfsNodes
    return nodesResults

def main(argv):
    # Initialize from command-line parameters:
    if len(argv) < 6:
        syntax()
    dbServer = argv[1]
    networkName = argv[2]
    userName = argv[3]
    password = argv[4]
    filename = argv[5]
    
    gtfsNodesResults = pathMatch(dbServer, networkName, userName, password, filename)
    
    # Extract useful information:
    print("INFO: -- Final --", file = sys.stderr)
    print("INFO: Print output...", file = sys.stderr)
    path_engine.dumpStandardHeader()

    datafileIDs = gtfsNodesResults.keys()
    "@type datafileIDs: list<str>"
    datafileIDs.sort()
    for datafileID in datafileIDs:
        "@type datafileID: str"
        path_engine.dumpStandardInfo(gtfsNodesResults[datafileID])
        
# Boostrap:
if __name__ == '__main__':
    main(sys.argv)
