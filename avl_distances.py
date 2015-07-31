"""
avl_distances.py performs a distance analysis over a given shape and
outputs the distance traveled from the start of the shape.
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
from nmc_mm_lib import gtfs, path_engine, graph
import problem_report, transit_gtfs, sys, csv, math
from datetime import datetime, timedelta

problemReport = False
"@var problemReport is set to true when the -p parameter is specified."

def syntax(exitCode):
    """
    Print usage information
    """
    print("avl_distances performs a distance analysis over a given shape and")
    print("outputs the distance traveled from the start of the shape")
    print()
    print("Usage:")
    print("  python avl_distances.py dbServer network user password shapePath")
    print("    pathMatchFile -a avlCSVFile [-r routeID] [-h headsign] [-p] [-s]")
    print()
    print("where:")
    print("  -a is the AVL CSV file.")
    print("  -r is the route ID of interest within the AVL CSV file.")
    print("  -h is the headsign of interest within the AVL CSV file. (Route ID and")
    print("     headsign form a unique ID).")
    print("  -p outputs a problem report on the stop matches")
    print("  -s outputs distances for stops that are in the GTFS set")
    sys.exit(exitCode)

def readAVLCSV(avlCSVFile, gtfsTrips, gps, routeID=None, routeHeadsign=None):
    """
    Reads AVL entries from the given CSV file and uses existing bus stop architecture to "fake" AVL entries
    as bus stops.
    @type avlCSVFile: str
    @param gtfsTrips: All valid GTFS trips
    @type gtfsTrips: dict<int, gtfs.TripsEntry>
    @param gps: gps.GPS
    @type routeID: int
    @type routeHeadsign: str
    @return A dictionary of trips to lists of StopTimesEntry
    @rtype dict<gtfs.TripsEntry, list<gtfs.StopTimesEntry>>
    """
    ret = {}
    "@type ret: dict<gtfs.TripsEntry, list<gtfs.StopTimesEntry>>"
    with open(avlCSVFile, 'r') as inFile:
        csvReader = csv.DictReader(inFile)
        prevTime = None
        prevRouteID = None
        prevTripID = None
        prevRouteHeadsign = None
        duplicateMsgFlag = False
        duplicateTimes = set()
        "@type duplicateTimes: set<str>"
        previousTripIDs = set()
        "@type previousTripIDs: set<str>"
        ctr = 0

        firstRun = True
        for fileLine in csvReader:
            if firstRun:
                # Sanity check:
                if not (all(x in fileLine for x in ["vehicle_id", "dist_traveled", "speed", "lon", "route_id", "trip_headsign",
                                                    "timestamp", "lat", "trip_id"])):
                    print("ERROR: The AVL CSV file %s doesn't have the expected header." % avlCSVFile, file=sys.stderr)
                    return None
                firstRun = False
        
            if (routeID is None or int(fileLine["route_id"]) == routeID) and (routeHeadsign is None or fileLine["trip_headsign"] == routeHeadsign):
                if prevRouteID is not None and prevRouteID != int(fileLine["route_id"]) or \
                        prevRouteHeadsign is not None and prevRouteHeadsign != fileLine["trip_headsign"]:
                    if not duplicateMsgFlag:
                        print("WARNING: Only one unique route ID and trip headsign from the AVL CSV file are allowed to be processed at once. "
                              "There was ambiguity at route ID %s, trip headsign %s." % (fileLine["route_id"], fileLine["trip_headsign"]), file=sys.stderr)
                        duplicateMsgFlag = True
                    continue
                prevRouteID = int(fileLine["route_id"])
                prevRouteHeadsign = fileLine["trip_headsign"]
                
                if prevTripID is None or int(fileLine["trip_id"]) != prevTripID:
                    if int(fileLine["trip_id"]) in previousTripIDs:
                        print("WARNING: In the AVL CSV input, Trip ID %s cannot be continued after going to another Trip ID." % fileLine["trip_id"], file=sys.stderr)
                        # TODO: Allow these to be out of order.
                        continue
                    previousTripIDs.add(int(fileLine["trip_id"]))
                    prevTripID = int(fileLine["trip_id"])
                    ctr = 0
                    prevTime = None
                
                # The parsing of the date is set up for the format: YYYY-MM-DDTHH:MM:SS-/+HH:MM
                timeParts = fileLine["timestamp"].split("T")
                datePart = datetime.strptime(timeParts[0], "%Y-%m-%d")
                
                # TODO: We shall ignore the time zone part right now and assume that other time references are in the current time zone.
                for searchChar in ["+", "-"]:
                    tzPos = timeParts[1].find(searchChar)
                    if tzPos > 0:
                        timeParts[1] = timeParts[1][0:tzPos]
                        break                
                timePart = datetime.strptime(timeParts[1], "%H:%M:%S")
                ourTime = datePart + timedelta(hours=timePart.hour, minutes=timePart.minute, seconds=timePart.second)
                
                if prevTime is not None and ourTime < prevTime:
                    if int(fileLine["trip_id"]) not in duplicateTimes:
                        print("WARNING: A non-increasing timestamp was discovered in the AVL CSV file %s, Trip %s; ignoring." % (avlCSVFile,
                            fileLine["trip_id"]), file=sys.stderr)
                        duplicateTimes.add(int(fileLine["trip_id"]))
                    continue
                prevTime = ourTime
                
                if int(fileLine["trip_id"]) not in gtfsTrips:
                    if int(fileLine["trip_id"]) not in duplicateTimes:
                        print("WARNING: Trip ID %s from the AVL CSV file is not found in the GTFS set." % fileLine["trip_id"], file=sys.stderr)
                        duplicateTimes.add(int(fileLine["trip_id"]))
                    continue
                gtfsTrip = gtfsTrips[int(fileLine["trip_id"])]
                
                # Here we fabricate fake stops for each AVL point:
                stop = gtfs.StopsEntry(ctr, fileLine["speed"], float(fileLine["lat"]), float(fileLine["lon"]))
                stop.pointX, stop.pointY = gps.gps2feet(stop.gpsLat, stop.gpsLng)
                stopTime = gtfs.StopTimesEntry(gtfsTrip, stop, ctr)
                stopTime.arrivalTime = ourTime
                stopTime.departureTime = ourTime
                ctr += 1
                
                # Store the stop results:
                if gtfsTrips[int(fileLine["trip_id"])] not in ret:
                    ret[gtfsTrips[int(fileLine["trip_id"])]] = []
                ret[gtfsTrips[int(fileLine["trip_id"])]].append(stopTime)
                    
    # Return the fake stop times:
    return ret

def dumpAVLDistances(gtfsTrips, gtfsStopTimes, gtfsNodes, vistaNetwork, stopSearchRadius, problemReport, stopsFlag=False,
                     outFile=sys.stdout):
    """
    dumpAVLDistances writes out AVL distances and speeds at each position along AVL paths.
    @type gtfsTrips: dict<int, gtfs.TripsEntry>
    @param gtfsStopTimes: A dictionary of Trip values to lists of StopTimesEntry
    @type gtfsStopTimes: dict<gtfs.TripsEntry, list<gtfs.StopTimesEntry>>
    @type gtfsNodes: dict<int, list<path_engine.PathEnd>>
    @type vistaNetwork: graph.GraphLib
    @type stopSearchRadius: float
    @type problemReport: bool
    @param stopsFlag: Set this true to format the output as though gtfsStopTimes holds the bus stop information
            rather than the "fake stop" AVL data.
    @type stopsFlag: bool
    @type outFile: file
    @return A mapping of tripID to AVL points-on-links
    @rtype dict<int, graph.PointOnLink>
    """
    # Set up the output:
    ret = {}
    "@type ret: dict<int, list<graph.PointOnLink>>"
        
    # Initialize the path engine for use later:
    pathEngine = path_engine.PathEngine(stopSearchRadius, stopSearchRadius, stopSearchRadius, sys.float_info.max, sys.float_info.max,
                                        stopSearchRadius, 1, 1.5, 1.5, sys.maxint, sys.maxint)
    pathEngine.limitClosestPoints = 12
    pathEngine.limitSimultaneousPaths = 12
    pathEngine.maxHops = 50
    pathEngine.logFile = None # Suppress the log outputs for the path engine; enough stuff will come from other sources.

    problemReportNodes = {}
    "@type problemReportNodes: dict<str, path_engine.PathEnd>"

    # Output header:
    if not problemReport:    
        if not stopsFlag:
            print("tripID,distance,timestamp,speed", file=outFile)
        else:
            print("tripID,stopID,stopSeq,distance,arrival,departure,name", file=outFile)
            
    tripIDs = [trip.tripID for trip in gtfsStopTimes.keys()]
    tripIDs.sort()
    for tripID in tripIDs:
        if gtfsTrips[tripID].shapeEntries[0].shapeID not in gtfsNodes:
            # This happens if the incoming files contain a subset of all available topology.
            print("WARNING: Skipping route for trip %d because no points are available." % tripID, file=sys.stderr)
            continue
        
        stopTimes = gtfsStopTimes[gtfsTrips[tripID]]
        "@type stopTimes: list<gtfs.StopTimesEntry>"
        
        # Step 1: Find the longest distance of contiguous valid links within the shape for each trip:
        # Step 2: Ignore routes that are entirely outside our valid time interval.
        ourGTFSNodes, _ = transit_gtfs.treeContiguous(gtfsNodes[gtfsTrips[tripID].shapeEntries[0].shapeID], vistaNetwork, 
            stopTimes)
        if ourGTFSNodes is None:
            continue        
            
        # Step 3: Match up stops to that contiguous list:
        # At this point, we're doing something with this.
        print("INFO: -- Matching stops for trip %d --" % tripID, file=sys.stderr)
        vistaSubset, _ = transit_gtfs.buildSubset(ourGTFSNodes, vistaNetwork)
        linkList = vistaSubset.linkMap.values()
        startPointOnLink = None
        endPointOnLink = None
        if len(linkList) > 0:
            # Express the exact start and end locations in terms of links that were created for the subset.
            startPointOnLink = graph.PointOnLink(linkList[0], ourGTFSNodes[0].pointOnLink.dist, ourGTFSNodes[0].pointOnLink.nonPerpPenalty,
                ourGTFSNodes[0].pointOnLink.refDist)
            startPointOnLink.pointX, startPointOnLink.pointY = ourGTFSNodes[0].pointOnLink.pointX, ourGTFSNodes[0].pointOnLink.pointY 
            endPointOnLink = graph.PointOnLink(linkList[-1], ourGTFSNodes[-1].pointOnLink.dist, ourGTFSNodes[-1].pointOnLink.nonPerpPenalty,
                ourGTFSNodes[-1].pointOnLink.refDist)
            endPointOnLink.pointX, endPointOnLink.pointY = ourGTFSNodes[-1].pointOnLink.pointX, ourGTFSNodes[-1].pointOnLink.pointY 

        # Then, prepare the stops as GTFS shapes entries:
        print("INFO: Mapping stops to VISTA network...", file=sys.stderr)
        gtfsShapes, _ = transit_gtfs.prepareMapStops(ourGTFSNodes, stopTimes, False)
        #gtfsShapes, gtfsStopsLookup = transit_gtfs.prepareMapStops(ourGTFSNodes, stopTimes)

        # Find a path through our prepared node map subset:
        resultTree = pathEngine.constructPath(gtfsShapes, vistaSubset, startPointOnLink, endPointOnLink)
        #resultTree = pathEngine.constructPath(gtfsShapes, vistaSubset)
        "@type resultTree: list<path_engine.PathEnd>"
        
        # Strip off the dummy ends:
        del resultTree[-1]
        del resultTree[0]
        if len(resultTree) > 0:
            resultTree[0].prevTreeNode = None
        
        # So now we should have one tree entry per matched stop.

        # Deal with Problem Report:
        # TODO: The Problem Report will include all nodes on each path regardless of valid time interval;
        # However; we will not have gotten here if the trip was entirely outside of it. 
        if problemReport:
            problemReportNodes[tripID if not stopsFlag else 0] = transit_gtfs.assembleProblemReport(resultTree, vistaNetwork) 
        
        # Dump out the output:
        distance = 0
        ctr = 0
        lastPos = ourGTFSNodes[0].pointOnLink.dist
        lastLink = ourGTFSNodes[0].pointOnLink.link
        lastPoint = ourGTFSNodes[0].pointOnLink
        
        index = 0
        for pathEnd in resultTree:
            "@type pathEnd: path_engine.PathEnd"
            if not pathEnd.restart:
                procLinks = pathEnd.routeInfo[:]
                procLinks.append(pathEnd.pointOnLink.link)
                for link in procLinks:
                    if link.id != lastLink.id: # Watch out! Link objects for vistaSubset are different instanciations than those in PathEnds.
                        distance += lastLink.distance - lastPos
                        lastPos = 0
                    else:
                        distance += pathEnd.pointOnLink.dist - lastPos
                        lastPos = pathEnd.pointOnLink.dist
                    lastLink = link
            else:
                print("WARNING: Trip ID %d, stop index %d had restarted; using direct distance." % (tripID, index), file=sys.stderr)
                distance += math.sqrt((pathEnd.pointOnLink.pointX - lastPoint.pointX) ** 2 + (pathEnd.pointOnLink.pointY - lastPoint.pointY) ** 2)
                lastLink = pathEnd.pointOnLink.link
            lastPoint = pathEnd.pointOnLink
            index += 1

            if not problemReport:            
                if not stopsFlag:
                    # Recall that because of transit_gtfs.prepareMapStops() the shapeID is the trips ID.
                    # And, recall that because of readAVLCSV(), stopName is the speed value. 
                    print ("%d,%f,%s,%s" % (pathEnd.shapeEntry.shapeID, distance, stopTimes[ctr].arrivalTime.strftime("%Y-%m-%dT%H:%M:%S"),
                        stopTimes[ctr].stop.stopName), file=outFile);
                else:
                    print ("%d,%d,%d,%f,%s,%s,%s" % (pathEnd.shapeEntry.shapeID, stopTimes[ctr].stop.stopID, stopTimes[ctr].stopSeq, distance,
                        stopTimes[ctr].arrivalTime.strftime("%H:%M:%S"), stopTimes[ctr].departureTime.strftime("%H:%M:%S"), 
                        stopTimes[ctr].stop.stopName), file=outFile);
                    # TODO: Note that the GTFS stopTimes input uses hours greater than 23 to express next early morning service.
                    # In gtfs.fillStopTimes(), this has been adapted to datetime by incrementing the day and doing a mod 24 on the hours.
                    # Already, stop times are stored relative to the day the epoch 1/1/1900. Maybe a new function call will do the trick.
            ctr += 1
        ret[tripID] = resultTree
    
    # Deal with Problem Report:
    if problemReport:
        print("INFO: Output problem report CSV...", file = sys.stderr)
        problemReportNodesOut = {}
        for idVal in problemReportNodes:
            seqs = problemReportNodes[idVal].keys()
            seqs.sort()
            ourTgtList = []
            for seq in seqs:
                ourTgtList.append(problemReportNodes[idVal][seq])
            problemReportNodesOut[idVal] = ourTgtList                
        problem_report.problemReport(problemReportNodesOut, vistaNetwork, True)
    
    return ret 

def main(argv):
    global problemReport
    
    # Initialize from command-line parameters:
    if len(argv) < 7:
        syntax(1)
    dbServer = argv[1]
    networkName = argv[2]
    userName = argv[3]
    password = argv[4]
    shapePath = argv[5]
    pathMatchFilename = argv[6]
    routeID = None
    routeHeadsign = None
    avlCSVFile = None
    stopsFlag = False
    problemReport = False

    if len(argv) > 6:
        i = 7
        while i < len(argv):
            if argv[i] == "-a" and i < len(argv) - 1:
                avlCSVFile = argv[i + 1]
                i += 1
            elif argv[i] == "-r" and i < len(argv) - 1:
                routeID = int(argv[i + 1])
                i += 1
            elif argv[i] == "-h" and i < len(argv) - 1:
                routeHeadsign = argv[i + 1]
                i += 1
            elif argv[i] == "-s":
                stopsFlag = True
            elif argv[i] == "-p":
                problemReport = True
            i += 1
    
    if avlCSVFile is None:
        print("ERROR: No AVL CSV file is specified. You must use the -a parameter.", file=sys.stderr)
        syntax(1)
    
    # Default parameters:
    stopSearchRadius = 1200
    
    # Restore the stuff that was built with path_match:
    vistaGraph, gtfsShapes, gtfsNodes, unusedShapeIDs = transit_gtfs.restorePathMatch(dbServer, networkName, userName,
        password, shapePath, pathMatchFilename)
    
    # Read in the routes information:
    print("INFO: Read GTFS routesfile...", file=sys.stderr)
    gtfsRoutes = gtfs.fillRoutes(shapePath)
    "@type gtfsRoutes: dict<int, RoutesEntry>"
    
    # Read in the trips information:
    print("INFO: Read GTFS tripsfile...", file=sys.stderr)
    gtfsTrips, unusedTripIDs = gtfs.fillTrips(shapePath, gtfsShapes, gtfsRoutes, unusedShapeIDs)
    "@type gtfsTrips: dict<int, TripsEntry>"
    "@type unusedTripIDs: set<int>"
        
    # Read in the AVL CSV file and create fake stops off of each point so we can use the existing code.
    print("INFO: Read AVL file %s..." % avlCSVFile, file=sys.stderr)
    gtfsStopTimes = readAVLCSV(avlCSVFile, gtfsTrips, vistaGraph.gps, routeID, routeHeadsign)
    if len(gtfsStopTimes.keys()) == 0:
        print("ERROR: No AVL points were found. Check to make sure that your search criteria in -r and/or -h are present in your AVL file.", file=sys.stderr)
        sys.exit(1)
    if stopsFlag:
        # Evidently there's the possibility of more trips in the GTFS trips file than expressed in the AVL for
        # the same shapes. Add entries to unusedTripIDs.
        tripsDelSet = set(gtfsTrips.keys())
        for avlTrip in gtfsStopTimes.keys():
            tripsDelSet.remove(avlTrip.tripID)
        unusedTripIDs = unusedTripIDs.union(tripsDelSet)
        del gtfsStopTimes # We don't need this from the AVL data anymore; it can be huge.
        
        # Read in the stops information:
        print("INFO: Read GTFS stopsfile...", file=sys.stderr)
        gtfsStops = gtfs.fillStops(shapePath, vistaGraph.gps)
        "@type gtfsStops: dict<int, StopsEntry>"
        
        # Read stop times information:
        print("INFO: Read GTFS stop times...", file=sys.stderr)
        gtfsStopTimes = gtfs.fillStopTimes(shapePath, gtfsTrips, gtfsStops, unusedTripIDs)
        "@type gtfsStopTimes: dict<gtfs.TripsEntry, list<StopTimesEntry>>"
        
    # Output the route distance information:
    if not problemReport:
        print("INFO: Outputting AVL distance information...", file=sys.stderr)
    dumpAVLDistances(gtfsTrips, gtfsStopTimes, gtfsNodes, vistaGraph, stopSearchRadius, problemReport, stopsFlag)

    print("INFO: Done.", file=sys.stderr)

# Boostrap:
if __name__ == '__main__':
    main(sys.argv)
