"""
transit_gtfs.py outputs a series of CSV files in the current path that
    are used for VISTA analysis of transit paths. In short, map-matched
    transit paths are the basis for map-matched bus stops. These stops
    are mapped to the underlying VISTA network, but because of a
    limitation in VISTA, only one of these stops may be mapped to any
    one link. (Extra stops for the time being are dropped).
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
from nmc_mm_lib import gtfs, vista_network, path_engine, graph
import problem_report, sys, time
from datetime import datetime, timedelta

DWELLTIME_DEFAULT = 0
"@var DWELLTIME_DEFAULT: The dwell time to report in the bus_route_link.csv file output."

STOP_SEARCH_RADIUS = 800
"@var STOP_SEARCH_RADIUS: 'k': Radius (ft) to search from GTFS point to perpendicular VISTA links"

DISTANCE_FACTOR = 1.0
"@var DISTANCE_FACTOR: 'f_d': Cost multiplier for Linear path distance in stop matching"

DRIFT_FACTOR = 2.0
"@var DRIFT_FACTOR: 'f_r': Cost multiplier for distance from GTFS point to its VISTA link in stop matching"

NON_PERP_PENALTY = 1.5
"@var NON_PERP_PENALTY: 'f_p': Penalty multiplier for GTFS points that aren't perpendicular to VISTA links"

problemReport = False
"@var problemReport is set to true when the -p parameter is specified."

def syntax(exitCode):
    """
    Print usage information
    """
    print("transit_gtfs outputs a series of CSV files in the current path that")
    print("are used for VISTA analysis of transit paths.")
    print()
    print("Usage:")
    print("  python transit_gtfs.py dbServer network user password shapePath")
    print("    pathMatchFile -t refDateTime [-e endTime] {[-c serviceID]")
    print("    [-c serviceID] ...} [-u] [-w] [-p]")
    print()
    print("where:")
    print("  -t is the zero-reference time that all arrival time outputs are related to.")
    print("     (Note that the day is ignored.) Use the format HH:MM:SS.")
    print("  -e is the duration in seconds (86400 by default). -t and -e filter stops.")
    print("  -c restricts results to specific service IDs (default: none)")
    print("  -u excludes links upstream of the first valid stop")
    print("  -w, -wb, -we: widen both, widen begin, widen end: include entire routes that")
    print("     would otherwise be cut off by -t (begin) and/or -e (end). This will")
    print("     suggest new starting time/duration and record all times relative to that.")
    print("  -x, -xb, -xe: exclude both, exclude begin, exclude end: excludes entire")
    print("     entire routes that intersect with -t (begin) and/or -e (end).")
    print("  -p outputs a problem report on the stop matches")
    sys.exit(exitCode)

def restorePathMatch(dbServer, networkName, userName, password, shapePath, pathMatchFilename):
    # Get the database connected:
    print("INFO: Connect to database...", file=sys.stderr)
    database = vista_network.connect(dbServer, userName, password, networkName)
    
    # Read in the topology from the VISTA database:
    print("INFO: Read topology from database...", file=sys.stderr)
    vistaGraph = vista_network.fillGraph(database)
    
    # Read in the shapefile information:
    print("INFO: Read GTFS shapefile...", file=sys.stderr)
    gtfsShapes = gtfs.fillShapes(shapePath, vistaGraph.gps)

    # Read the path-match file:
    print("INFO: Read the path-match file '%s'..." % pathMatchFilename, file=sys.stderr)
    with open(pathMatchFilename, 'r') as inFile:
        gtfsNodes = path_engine.readStandardDump(vistaGraph, gtfsShapes, inFile)
        "@type gtfsNodes: dict<int, list<path_engine.PathEnd>>"

    # Filter out the unused shapes:
    unusedShapeIDs = set()
    for shapeID in gtfsShapes.keys():
        if shapeID not in gtfsNodes:
            del gtfsShapes[shapeID]
            unusedShapeIDs.add(shapeID)

    return (vistaGraph, gtfsShapes, gtfsNodes, unusedShapeIDs)

def _outHeader(tableName, userName, networkName, outFile):
    print("User,%s" % userName, file = outFile)
    print("Network,%s" % networkName, file = outFile)
    print("Table,public.bus_route", file = outFile)
    print(time.strftime("%a %b %d %H:%M:%S %Y"), file = outFile)
    print(file = outFile)

def dumpBusRoutes(gtfsTrips, userName, networkName, outFile = sys.stdout):
    """
    dumpBusRoutes dumps out a public.bus_route.csv file contents.
    @type gtfsTrips: dict<int, gtfs.TripsEntry>
    @type userName: str
    @type networkName: str
    @type outFile: file
    """
    _outHeader("public.bus_route", userName, networkName, outFile)
    print("\"id\",\"name\",", file = outFile)
    
    # Remember, we are treating each route as a trip.
    tripIDs = gtfsTrips.keys()
    tripIDs.sort()
    for tripID in tripIDs:
        append = ""
        if len(gtfsTrips[tripID].route.name) > 0:
            append = ": " + gtfsTrips[tripID].route.name
        if len(gtfsTrips[tripID].tripHeadsign) > 0:
            append += " " + gtfsTrips[tripID].tripHeadsign
        print("\"%d\",\"%s\"" % (tripID, gtfsTrips[tripID].route.shortName + append),
                file = outFile)

def treeContiguous(treeNodes, vistaNetwork, gtfsStopTimes=None, startTime=None, endTime=None):
    """
    treeContiguous finds the largest consecutive string of continuous points that fall within the
    given time range (if given). Used by dumpBusRouteLinks() and possibly others.
    @type treeNodes: list<path_engine.PathEnd>
    @type vistaNetwork: graph.GraphLib
    @type gtfsStopTimes: list<gtfs.StopTimesEntry>
    @type startTime: datetime
    @type endTime: datetime
    @return A list of tree entries that are extracted from the given shape, or None if failure.
    @rtype (list<path_engine.PathEnd>, int)
    """
    startIndex = -1
    curIndex = 0
    linkCount = 0
    totalLinks = 0
    
    longestStart = -1
    longestEnd = len(treeNodes)
    longestDist = sys.float_info.min
    longestLinkCount = 0
    
    while curIndex <= len(treeNodes):
        if curIndex == len(treeNodes) or curIndex == 0 or treeNodes[curIndex].restart:
            totalLinks += 1
            linkCount += 1
            if curIndex > startIndex and startIndex >= 0:
                # We have a contiguous interval.  See if it wins:
                if treeNodes[curIndex - 1].totalDist - treeNodes[startIndex].totalDist > longestDist:
                    longestStart = startIndex
                    longestEnd = curIndex
                    longestDist = treeNodes[curIndex - 1].totalDist - treeNodes[startIndex].totalDist
                    longestLinkCount = linkCount
                    linkCount = 0
                
            # This happens if it is time to start a new interval:
            startIndex = curIndex
        else:
            totalLinks += len(treeNodes[curIndex].routeInfo)
            linkCount += len(treeNodes[curIndex].routeInfo)
        curIndex += 1

    if longestStart >= 0:
        # We have a valid path.  See if it had been trimmed down and report it.
        if longestStart > 0 or longestEnd < len(treeNodes):
            print("WARNING: For shape ID %s from seq. %d through %d, %.2g%% of %d links will be used used because of restarts in the path match file.." \
                  % (str(treeNodes[longestStart].shapeEntry.shapeID), treeNodes[longestStart].shapeEntry.shapeSeq,
                     treeNodes[longestEnd - 1].shapeEntry.shapeSeq, 100 * float(longestLinkCount) / float(totalLinks),
                     totalLinks), file=sys.stderr)
        
        # Isolate the relevant VISTA tree nodes: (Assume from above that this is a non-zero length array)
        return treeNodes[longestStart:longestEnd], longestStart
    
    else:
        print("WARNING: No links for shape %s." % str(treeNodes[longestStart].shapeEntry.shapeID), file=sys.stderr)
        return None, 0

def buildSubset(treeNodes, vistaNetwork):
    """
    buildSubset builds a tree from a shape. Used by dumpBusRouteLinks() and possibly others.
    @type treeNodes: list<path_engine.PathEnd>
    @type vistaNetwork: graph.GraphLib
    @return The network subset and list of link IDs
    @rtype (graph.GraphLib, list<graph.GraphLink>)
    """
    # We are going to recreate a small VISTA network from ourGTFSNodes and then match up the stops to that.
    # First, prepare the small VISTA network:
    subset = graph.GraphLib(vistaNetwork.gps.latCtr, vistaNetwork.gps.lngCtr)
    
    # Build a list of links:
    outLinkList = []
    "@type outLinkList: list<graph.GraphLink>"
    
    # Plop in the start node:
    subsetNodePrior = graph.GraphNode(treeNodes[0].pointOnLink.link.origNode.id,
        treeNodes[0].pointOnLink.link.origNode.gpsLat, treeNodes[0].pointOnLink.link.origNode.gpsLng)
    "@type subsetNodePrior: graph.GraphNode"
    subsetNodePrior.coordX, subsetNodePrior.coordY = treeNodes[0].pointOnLink.link.origNode.coordX, treeNodes[0].pointOnLink.link.origNode.coordY
    subset.addNode(subsetNodePrior)
    prevLinkID = treeNodes[0].pointOnLink.link.id
    
    # Link together nodes as we traverse through them:
    for ourGTFSNode in treeNodes:
        "@type ourGTFSNode: path_engine.PathEnd"
        # There should only be one destination link per VISTA node because this comes form our tree.
        # If there is no link or we're repeating the first one, then there were no new links assigned.
        if len(ourGTFSNode.routeInfo) < 1 or (len(outLinkList) == 1 \
                and ourGTFSNode.routeInfo[0].id == treeNodes[0].pointOnLink.link.id):
            continue
        for link in ourGTFSNode.routeInfo:
            "@type link: graph.GraphLink"
        
            if link.id not in vistaNetwork.linkMap:
                print("WARNING: In finding bus route links, link ID %d is not found in the VISTA network." % link.id, file=sys.stderr)
                continue
            origVistaLink = vistaNetwork.linkMap[link.id]
            "@type origVistaLink: graph.GraphLink"
            
            # Create a new node, even if the node had been visited before. We are creating a single-path and need a separate instance:
            subsetNode = graph.GraphNode(origVistaLink.origNode.id, origVistaLink.origNode.gpsLat, origVistaLink.origNode.gpsLng)
            subsetNode.coordX, subsetNode.coordY = origVistaLink.origNode.coordX, origVistaLink.origNode.coordY
            subset.addNode(subsetNode)
                
            # We shall label our links as indices into the stage we're at in ourGTFSNodes links.  This will allow for access later.
            newLink = graph.GraphLink(prevLinkID, subsetNodePrior, subsetNode)
            subset.addLink(newLink)
            outLinkList.append(newLink)
            subsetNodePrior = subsetNode
            prevLinkID = link.id
            
    # And then finish off the graph with the last link:
    subsetNode = graph.GraphNode(ourGTFSNode.pointOnLink.link.destNode.id, ourGTFSNode.pointOnLink.link.destNode.gpsLat, ourGTFSNode.pointOnLink.link.destNode.gpsLng)
    subsetNode.coordX, subsetNode.coordY = ourGTFSNode.pointOnLink.link.destNode.coordX, ourGTFSNode.pointOnLink.link.destNode.coordY
    subset.addNode(subsetNode)
    newLink = graph.GraphLink(prevLinkID, subsetNodePrior, subsetNode)
    subset.addLink(newLink)
    outLinkList.append(newLink)
    
    return subset, outLinkList

def prepareMapStops(treeNodes, stopTimes, dummyFlag=True):
    """
    prepareMapStops maps stops information to an underlying path. Used by dumpBusRouteLinks() and possibly others.
    @type treeNodes: list<path_engine.PathEnd>
    @type stopTimes: list<gtfs.StopTimesEntry>
    @param dummyFlag: Set to true to add dummy entries to the start and end. This will allow proximity searching at these places.
    @type dummyFlag: bool
    @return Prepared stop information
    @rtype (list<gtfs.ShapesEntry>, dict<int, gtfs.StopTimesEntry>)
    """
    gtfsShapes = []
    gtfsStopsLookup = {}
    "@type gtfsStopsLookup: dict<int, gtfs.StopTimesEntry>"
    
    if dummyFlag:
        # Append an initial dummy shape to force routing through the path start:
        newShapesEntry = gtfs.ShapesEntry(stopTimes[0].trip.tripID, -1, treeNodes[0].pointOnLink.link.origNode.gpsLat,
            treeNodes[0].pointOnLink.link.origNode.gpsLng)
        newShapesEntry.pointX, newShapesEntry.pointY = treeNodes[0].pointOnLink.link.origNode.coordX, treeNodes[0].pointOnLink.link.origNode.coordY
        newShapesEntry.typeID = 1 # To signify that this is like a stop.
        gtfsShapes.append(newShapesEntry)
    
    # Append all of the stops:
    for gtfsStopTime in stopTimes:
        "@type gtfsStopTime: gtfs.StopTimesEntry"
        newShapesEntry = gtfs.ShapesEntry(gtfsStopTime.trip.tripID, gtfsStopTime.stopSeq, gtfsStopTime.stop.gpsLat, gtfsStopTime.stop.gpsLng)
        newShapesEntry.pointX, newShapesEntry.pointY = gtfsStopTime.stop.pointX, gtfsStopTime.stop.pointY
        newShapesEntry.typeID = 1 # To signify that this is a stop.
        gtfsShapes.append(newShapesEntry)
        gtfsStopsLookup[gtfsStopTime.stopSeq] = gtfsStopTime

    if dummyFlag:
        # Append a trailing dummy shape to force routing through the path end:
        newShapesEntry = gtfs.ShapesEntry(stopTimes[0].trip.tripID, -1, treeNodes[-1].pointOnLink.link.destNode.gpsLat,
            treeNodes[-1].pointOnLink.link.destNode.gpsLng)
        newShapesEntry.pointX, newShapesEntry.pointY = treeNodes[-1].pointOnLink.link.destNode.coordX, treeNodes[-1].pointOnLink.link.destNode.coordY
        newShapesEntry.typeID = 1 # To signify that this is like a stop.
        gtfsShapes.append(newShapesEntry)

    return gtfsShapes, gtfsStopsLookup

def zipGTFSWithStops(underlyingNetwork, leftResultTree, rightResultTree):
    """
    zipGTFSWithStops combines together two series of graph.PointOnLink lists into one. The two series must overlap each other,
    the right one being at most as long as the left. Also expresses all links in PointOnLinks in terms of the left PathEnds.
    Please take care to keep the ability to identify which list the original PathEnds are by setting the respective gtfs.ShapesEntry's typeID value.
    @type underlyingNetwork: graph.GraphLib
    @type leftResultTree: list<path_engine.PathEnd>
    @type rightResultTree: list<path_engine.PathEnd>
    @rtype list<path_engine.PathEnd>
    """
    superResultTree = []
    "@type superResultTree: list<path_engine.PathEnd>"
    
    leftIndex = 0
    leftRouteInfoIndex = 0
    leftLinkID = -1
    leftDist = 0
    rightIndex = 0
    rightLinkID = -1
    recordRightNext = False
    currentRouteInfoList = []
    while leftIndex < len(leftResultTree) or rightIndex < len(rightResultTree):
        rightLinkID = rightResultTree[rightIndex].pointOnLink.link.id if rightIndex < len(rightResultTree) else -1
        while True:            
            # Obtain the next left link ID here:
            if leftIndex >= len(leftResultTree):
                leftLinkID = -1
            elif not leftResultTree[leftIndex].routeInfo:
                leftLinkID = leftResultTree[leftIndex].pointOnLink.link.id
            else:
                leftLinkID = leftResultTree[leftIndex].routeInfo[leftRouteInfoIndex].id
            if leftIndex >= len(leftResultTree) or leftRouteInfoIndex < len(leftResultTree[leftIndex].routeInfo) - 1:
                leftDist = 0
            else:
                leftDist = leftResultTree[leftIndex].pointOnLink.dist

            # Are we needing to record the right PathEnd?
            if recordRightNext:
                # Are we sure we want to record at this time? (Check to see if we still need to evaluate the next left PathEnd before recording the right PathEnd
                # because the next left PathEnd still comes before the right one).
                if not (leftIndex < len(leftResultTree) and leftRouteInfoIndex >= len(leftResultTree[leftIndex].routeInfo) - 1 and not currentRouteInfoList and superResultTree \
                        and superResultTree[-1].pointOnLink.link.id == rightResultTree[rightIndex].pointOnLink.link.id and \
                        leftDist < rightResultTree[rightIndex].pointOnLink.dist):
                    # For elements that are recorded for the right PathEnds, we want to express the links with link objects that are tied with the
                    # underlying network, not the manufactured subsets.
                    superPathEnd = path_engine.PathEnd(rightResultTree[rightIndex].shapeEntry, rightResultTree[rightIndex].pointOnLink.copy())
                    superPathEnd.pointOnLink.link = underlyingNetwork.linkMap[superPathEnd.pointOnLink.link.id]
                    superPathEnd.totalCost, superPathEnd.totalDist, superPathEnd.routeInfo = rightResultTree[rightIndex].totalCost, rightResultTree[rightIndex].totalDist, currentRouteInfoList
                    if superResultTree:
                        superPathEnd.prevTreeNode = superResultTree[-1] 
                    superResultTree.append(superPathEnd)
                    recordRightNext = False
                    currentRouteInfoList = []
                    break
                            
            # Keep a record of links we visit.
            if leftResultTree[leftIndex].routeInfo and not (superResultTree and leftLinkID == superResultTree[-1].pointOnLink.link.id):
                currentRouteInfoList.append(leftLinkID)
            
            # See if we need to record a right PathEnd.
            if leftLinkID == rightLinkID:
                recordRightNext = True
                if leftRouteInfoIndex >= len(leftResultTree[leftIndex].routeInfo) - 1:
                    # This is the last element, so we're needing to compare the distances.
                    if leftResultTree[leftIndex].pointOnLink.dist > rightResultTree[rightIndex].pointOnLink.dist:
                        # The left PathEnd comes after the right PathEnd, so record the right one first.
                        continue
                # Otherwise, record the left PathEnd now and record the right PathEnd when we compare again.
        
            # Are we needing to record the left PathEnd?
            if leftRouteInfoIndex >= len(leftResultTree[leftIndex].routeInfo) - 1:
                # This is the last element, so we need to consider recording it.
                # Does the point overlap the right point? Or, does it go backwards from the previous point?
                if not (recordRightNext and leftDist == rightResultTree[rightIndex].pointOnLink.dist \
                        or not currentRouteInfoList and superResultTree and superResultTree[-1].pointOnLink.link.id == leftResultTree[leftIndex].pointOnLink.link.id \
                        and superResultTree[-1].pointOnLink.dist > leftDist):
                    # No, we don't want to skip the left one. (Otherwise, it would have been redundant).
                    superPathEnd = leftResultTree[leftIndex].cleanCopy()
                    superPathEnd.totalCost, superPathEnd.totalDist, superPathEnd.routeInfo = leftResultTree[leftIndex].totalCost, leftResultTree[leftIndex].totalDist, currentRouteInfoList
                    if superResultTree:
                        superPathEnd.prevTreeNode = superResultTree[-1] 
                    superResultTree.append(superPathEnd)
                    currentRouteInfoList = []
                    
                # Move on to the next left one.
                leftRouteInfoIndex = 0
                leftIndex += 1
            else:
                leftRouteInfoIndex += 1
            
            if leftIndex >= len(leftResultTree) and not recordRightNext:
                break
        rightIndex += 1
    return superResultTree
    
def assembleProblemReport(resultTree, vistaNetwork):
    """
    assembleProblemReport puts together updated path_engine.PathEnd objects that are used for a problem report.
    Used by dumpBusRouteLinks() and possibly others.
    @type resultTree: list<path_engine.PathEnd>
    @type vistaNetwork: graph.GraphLib
    @rtype dict<int, path_engine.PathEnd>
    """
    revisedNodeList = {}
    prevNode = None
    "@type revisedNodeList: dict<int, path_engine.PathEnd>"
    for stopNode in resultTree:
        # TODO: Do we need to isolate this to shape typeID 1?
        
        # Reconstruct a tree node in terms of the original network.
        # TODO: Check to make sure that resultTree[0].shapeEntry.shapeID is correct.
        newShape = gtfs.ShapesEntry(resultTree[0].shapeEntry.shapeID,
            stopNode.shapeEntry.shapeSeq, stopNode.shapeEntry.lat, stopNode.shapeEntry.lng, False)
        origLink = vistaNetwork.linkMap[stopNode.pointOnLink.link.id] 
        newPointOnLink = graph.PointOnLink(origLink, stopNode.pointOnLink.dist,
            stopNode.pointOnLink.nonPerpPenalty, stopNode.pointOnLink.refDist)
        newNode = path_engine.PathEnd(newShape, newPointOnLink)
        newNode.restart = stopNode.restart
        newNode.totalCost = stopNode.totalCost
        newNode.totalDist = stopNode.totalDist
        newNode.routeInfo = []
        for link in stopNode.routeInfo:
            newNode.routeInfo.append(vistaNetwork.linkMap[link.id])
        newNode.prevTreeNode = prevNode
        prevNode = newNode
        revisedNodeList[stopNode.shapeEntry.shapeSeq] = newNode
    return revisedNodeList 

def dumpBusRouteLinks(gtfsTrips, gtfsStops, gtfsStopTimes, gtfsNodes, vistaNetwork, stopSearchRadius, excludeUpstream, userName,
        networkName, startTime, endTime, widenBegin, widenEnd, excludeBegin, excludeEnd, outFile=sys.stdout):
    """
    dumpBusRouteLinks dumps out a public.bus_route_link.csv file contents. This also will remove all stop times and trips
    that fall outside of the valid evaluation interval as dictated by the exclusion parameters.
    @type gtfsTrips: dict<int, gtfs.TripsEntry>
    @type gtfsStops: dict<int, StopsEntry>
    @type gtfsStopTimes: dict<TripsEntry, list<gtfs.StopTimesEntry>>
    @type gtfsNodes: dict<int, list<path_engine.PathEnd>>
    @type vistaNetwork: graph.GraphLib
    @type stopSearchRadius: float
    @type excludeUpstream: boolean
    @type userName: str
    @type networkName: str
    @type startTime: datetime
    @type endTime: datetime
    @type widenBegin: bool
    @type widenEnd: bool
    @type excludeBegin: bool
    @type excludeEnd: bool
    @type outFile: file
    @return A mapping of stopID to points-on-links plus the start and end times adjusted for
            warm-up and cool-down (if widenBegin or widenEnd is True)
    @rtype (dict<int, graph.PointOnLink>, datetime, datetime)
    """
    _outHeader("public.bus_route_link", userName, networkName, outFile)
    print('"route","sequence","link","stop","dwelltime",', file = outFile)
    
    # Set up the output:
    ret = {}
    "@type ret: dict<int, graph.PointOnLink>"
        
    warmupStartTime = startTime
    cooldownEndTime = endTime

    # Initialize the path engine for use later:
    pathEngine = path_engine.PathEngine(stopSearchRadius, stopSearchRadius, stopSearchRadius, sys.float_info.max, sys.float_info.max,
                                        stopSearchRadius, DISTANCE_FACTOR, DRIFT_FACTOR, NON_PERP_PENALTY, sys.maxint, sys.maxint)
    pathEngine.limitClosestPoints = 8
    pathEngine.limitSimultaneousPaths = 6
    pathEngine.maxHops = sys.maxint
    pathEngine.logFile = None # Suppress the log outputs for the path engine; enough stuff will come from other sources.

    class _TripsBundle:
        """
        _TripsBundle is a small structure for tying visible identifiers to lists of trips that are unique to each
        shape/stop sequence.
        @ivar trips: list<gtfs.TripsEntry>
        @ivar label: str
        @ivar longestStart: int
        @ivar resultTree: list<path_engine.PathEnd>"
        @ivar subset: graph.GraphLib
        @ivar stopsLookup: dict<int, gtfs.StopTimesEntry>
        @ivar forceLinks: list<set<graph.GraphLink>
        """
        def __init__(self):
            self.trips = []
            self.label = ""
            self.longestStart = 0
            self.resultTree = None
            self.subset = None
            self.stopsLookup = None
            self.forceLinks = None
    
    print("INFO: Perform the initial bus stop matching...", file=sys.stderr)
    shapeStops = {}
    "@type shapeStops: dict<int, dict<tuple<int>, _TripsBundle>>"
    
    tripID = trip = stopsList = shapeID = stopsTuples = stopsTuple = tripsBundle = emptyShapes = flag = None
    tripIDs = gtfsTrips.keys()
    tripIDs.sort()
    for tripID in tripIDs:
        trip = gtfsTrips[tripID]
        shapeID = trip.shapeEntries[0].shapeID
        "@type trip: gtfs.TripsEntry"
        if shapeID not in gtfsNodes:
            # This happens if the incoming files contain a subset of all available topology.
            print("WARNING: Skipping route for trip %d because no points are available." % tripID, file=sys.stderr)
            continue
        
        # Identify all elements of analysis: shape ID and stop ID combination. For each, we later want to analyze
        # all trips that use that combination as a unit.
        stopsList = [shapeID]
        "@type stopsList: list<int>"
        for stopTimesEntry in gtfsStopTimes[trip]:
            "@type stopTimesEntry: list<gtfs.StopTimesEntry>"
            stopsList.append(stopTimesEntry.stop.stopID)
        stopsTuple = tuple(stopsList)
        "@type stopsTuple: tuple<int>"
        # So now we have a key that is composed of the shape ID and all stop IDs in a sequence to identify all unique
        # shape/stop ID combinations.
                    
        # Ignore trips that are entirely outside our valid time interval.
        flag = False
        # TODO: Trips with no stops will be ignored. Is that what we want to do?
        if not gtfsStopTimes:
            # This will happen if we don't have stops defined. In this case, we want to go ahead and process the bus_route_link
            # outputs because we don't know if the trip falls in or out of the valid time range.
            flag = True
        else:
            for stopEntry in gtfsStopTimes[trip]:
                if (startTime is None or stopEntry.arrivalTime >= startTime) and (endTime is None or stopEntry.arrivalTime <= endTime):
                    flag = True
                    break
        if not flag:
            # This will be done silently because (depending upon the valid interval) there could be
            # hundreds of these in a GTFS set.
            continue
        
        if shapeID not in shapeStops:
            shapeStops[shapeID] = {}
        if stopsTuple not in shapeStops[shapeID]:
            shapeStops[shapeID][stopsTuple] = _TripsBundle()
        shapeStops[shapeID][stopsTuple].trips.append(trip)    
    del tripID, trip, stopsList, shapeID, stopsTuples, stopsTuple, tripsBundle, emptyShapes, flag
    
    shapeID = stopsTuples = stopsTuple = tripsBundle = varCounter = resultTree = None
    allTripsBundles = {}
    "@type allTripsBundles: dict<tuple<int>, _TripsBundle>"    
    for shapeID, stopsTuples in shapeStops.iteritems():
        "@type shapeID: int"
        "@type stopsTuples: dict<tuple<int>, _TripsBundle>"
        
        varCounter = 1
        for stopsTuple, tripsBundle in stopsTuples.iteritems():
            "@type stopsTuple: tuple<int>"
            "@type tripsBundle: _TripsBundle"
        
            tripsBundle.label = str(shapeID) if len(stopsTuples) <= 1 else str(shapeID) + " (variation %d)" % varCounter
                
            # Step 1: Find the longest distance of contiguous valid links within each unique shape. And,
            # TODO: Note that only the longest contiguous series is found; if there are two or more sections, then they'll be ignored.
            print("INFO: -- Matching stops for Shape %s: %d trip(s) --" % (tripsBundle.label, len(tripsBundle.trips)), file=sys.stderr)
            ourGTFSNodes, tripsBundle.longestStart = treeContiguous(gtfsNodes[shapeID], vistaNetwork)
            "@type ourGTFSNodes: list<path_engine.PathEnd>"
            
            # Step 3: Build a new subset network with new links and nodes that represents the single-path
            # specified by the GTFS shape (for bus route):
            tripsBundle.subset, _ = buildSubset(ourGTFSNodes, vistaNetwork)
                    
            # Step 4: Match up stops to that contiguous list:
            print("INFO: Mapping stops to VISTA network...", file=sys.stderr)
            gtfsShapes, tripsBundle.stopsLookup = prepareMapStops(ourGTFSNodes, gtfsStopTimes[tripsBundle.trips[0]])
            # Here, trip[0] is a representative trip that is identical to all other trips in this bundle. The call above is only
            # getting the stop sequence, not the times of the stops.
    
            # Find a path through our prepared node map subset:
            resultTree = pathEngine.constructPath(gtfsShapes, tripsBundle.subset)
            "@type resultTree: list<path_engine.PathEnd>"
    
            # So now resultTree is one tree entry per matched stop plus dummy ends.
            # Strip off the dummy ends:
            del resultTree[-1]
            del resultTree[0]
            if len(resultTree) > 0:
                resultTree[0].prevTreeNode = None
                resultTree[0].routeInfo = []
                
            # Check if our matched path has any problems:
            if sum([pathEnd.restart for pathEnd in resultTree]) > 0:
                print("WARNING: Skipping analysis of Shape %s because the bus stop matching resulted in disjointed sections." % tripsBundle.label, file=sys.stderr)
                continue
    
            # Put together GTFS points and stop points into the same list and express all links in terms of the underlying network:
            tripsBundle.resultTree = zipGTFSWithStops(vistaNetwork, ourGTFSNodes, resultTree)
    
            # While we're at it, initialize the force links list for later:
            tripsBundle.forceLinks = [None] * len(tripsBundle.resultTree)
            
            # Store this for the next step, where we attempt to use the same bus stop location across all routes that use
            # that bus stop.
            allTripsBundles[stopsTuple] = tripsBundle
   
            varCounter += 1
    del shapeID, stopsTuples, stopsTuple, tripsBundle, varCounter, resultTree
    
    # At this point, all shapes represented in allTripsBundles are ones we want to analyze, and those extra ones in shapeStops
    # have disjoints. We'll probably want to output both at the end.
    print("INFO: -- End matching stops --", file=sys.stderr)
            
    # Now figure out where stop locations differ among multiple routes that share the same stop. This will involve
    # filling out a StopRecord for each stop and then resolving the discrepancies.
    
    print("INFO: Resolving discrepancies in bus stop locations across all routes...", file=sys.stderr)    
    class StopRecord:
        """
        StopRecord is a container for storing all of the information about stops and links so that 
        the link that is used across all routes may be set to be the same.
        @ivar linkCounts: Stores a reference count for each link referring to this stop. LinkID -> reference count.
        @type linkCounts: dict<int, int>
        @ivar referents: Identifies for each trip the index in the respective treeEntry list addresses the stop.
                        stop tuple -> index into tree list
        @type referents: dict<tuple<int>, set<int>>
        @ivar refCount: The number of referents there are among all trips.
        @type refCount: int
        """
        def __init__(self):
            self.linkCounts = {}
            self.referents = {}
            self.refCount = 0
    
    stopRecords = {}
    "@type stopRecords: dict<int, StopRecord>"

    shapeID = stopsTuples = stopsTuple = tripsBundle = linkID = stopRecord = treeEntry = treeEntryIndex = None
    for shapeID, stopsTuples in shapeStops.iteritems():
        "@type shapeID: int"
        "@type stopsTuples: dict<tuple<int>, _TripsBundle>"
        
        for stopsTuple, tripsBundle in stopsTuples.iteritems():
            "@type stopsTuple: tuple<int>"
            "@type tripsBundle: _TripsBundle"
            
            # Is this one we skipped earlier?
            if stopsTuple not in allTripsBundles:
                continue
            
            for treeEntryIndex, treeEntry in enumerate(tripsBundle.resultTree):
                "@type treeEntry: path_engine.PathEnd"
                if treeEntry.shapeEntry.typeID == 1:
                    # Only process those entries that are associated with stop points.
                    stopID = tripsBundle.stopsLookup[treeEntry.shapeEntry.shapeSeq].stop.stopID
                    if stopID not in stopRecords:
                        stopRecords[stopID] = StopRecord()
                    stopRecord = stopRecords[stopID]
                    linkID = treeEntry.pointOnLink.link.id
                    # Use the ID here so we have a common reference among all trips on which link is what.
                    
                    # Count that this link is matched to this stop.
                    if linkID not in stopRecord.linkCounts:
                        stopRecord.linkCounts[linkID] = 0
                    stopRecord.linkCounts[linkID] += len(tripsBundle.trips)
                    
                    # Identify where in the resultTree matched stops list this matched stop occurs.
                    if stopsTuple not in stopRecord.referents:
                        stopRecord.referents[stopsTuple] = set()
                    stopRecord.referents[stopsTuple].add(treeEntryIndex)
                    
                    # Increment the counter for the total number of references.
                    stopRecord.refCount += len(tripsBundle.trips)
    del shapeID, stopsTuples, stopsTuple, tripsBundle, linkID, stopRecord, treeEntry, treeEntryIndex
    
    # Go through each stop and check for discrepancies.
    
    # We need a new PathEngine that has parameters similar to that of the path_refine module because we'll be looking against the entire 
    # underlying network.
    # TODO: Put these into a centralized location.
    termRefactorRadius = 1500   # Radius (ft) to invalidate found points at either end of a restart.
    pointSearchRadius = 300    # "k": Radius (ft) to search from GTFS point to perpendicular VISTA links
    pointSearchPrimary = 300   # "k_p": Radius (ft) to search from GTFS point to new VISTA links    
    pointSearchSecondary = 150  # "k_s": Radius (ft) to search from VISTA perpendicular point to previous point
    limitLinearDist = 1500      # Path distance (ft) to allow new proposed paths from one point to another
    limitDirectDist = 1500      # Radius (ft) to allow new proposed paths from one point to another
    limitDirectDistRev = 500    # Radius (ft) to allow backtracking on an existing link (e.g. parking lot)
    distanceFactor = 1.0        # "f_d": Cost multiplier for Linear path distance
    driftFactor = 2.0           # "f_r": Cost multiplier for distance from GTFS point to its VISTA link
    nonPerpPenalty = 1.5        # "f_p": Penalty multiplier for GTFS points that aren't perpendicular to VISTA links
    limitClosestPoints = 12     # "q_p": Number of close-proximity points that are considered for each GTFS point 
    limitSimultaneousPaths = 8  # "q_e": Number of proposed paths to maintain during pathfinding stage

    maxHops = 8                 # Maximum number of VISTA links to pursue in a path-finding operation
        
    # Initialize the path-finder:
    pathEngine = path_engine.PathEngine(pointSearchRadius, pointSearchPrimary, pointSearchSecondary, limitLinearDist,
                            limitDirectDist, limitDirectDistRev, distanceFactor, driftFactor, nonPerpPenalty, limitClosestPoints,
                            limitSimultaneousPaths)
    pathEngine.setRefineParams(termRefactorRadius)
    pathEngine.maxHops = maxHops    
    pathEngine.logFile = None # Suppress the log outputs for the path engine; enough stuff will come from other sources.
    del termRefactorRadius, pointSearchRadius, pointSearchPrimary, pointSearchSecondary, limitLinearDist, limitDirectDist, limitDirectDistRev, \
        distanceFactor, driftFactor, nonPerpPenalty, limitClosestPoints, limitSimultaneousPaths, maxHops
    
    pathEngine.setRefineParams(STOP_SEARCH_RADIUS)
    for stopID, stopRecord in stopRecords.iteritems():
        print("INFO: -- Stop %d --" % stopID, file=sys.stderr)
        if len(stopRecord.linkCounts) <= 1:
            if len(stopRecord.linkCounts) == 1:
                print("INFO: This has 1 matched linkID %d. OK." % (stopRecord.linkCounts.keys()[0]), file=sys.stderr)
            else:
                print("INFO: This has no links associated with it. Skipping.", file=sys.stderr)
            continue # All routes use the same link for the stop. No discrepancies for this stop.

        print("INFO: There are currently %d matched link(s) among %d trip(s). Disambiguating..." % (len(stopRecord.linkCounts), sum(stopRecord.linkCounts.values())),
              file=sys.stderr)
                
        # For each trip that uses this stop, discover closest points.
        stopsEntry = gtfsStops[stopID]
        "@type stopsEntry: StopsEntry"
        
        # Use a StopRecord to keep track of which stop and trips use which matched links
        stopResolveRecord = StopRecord()

        # This will hold the scores for all of the reference distances. 
        scores = {}
        "@type scores: dict<int, float>"                        
        
        closestLinks = pointOnLink = linkID = None
        for stopsTuple, treeEntryIndices in stopRecord.referents.iteritems():
            closestLinks = vistaNetwork.findPointsOnLinks(stopsEntry.pointX, stopsEntry.pointY, pathEngine.pointSearchRadius,
                pathEngine.pointSearchPrimary, pathEngine.pointSearchSecondary, None, pathEngine.limitClosestPoints)            
            "@type closestLinks: list<graph.PointOnLink>"
            
            if not closestLinks:
                print("WARNING: No closest points found for Shape %s index %d" % (allTripsBundles[stopsTuple].label, treeEntryIndices[0]), file=sys.stderr)
                stopResolveRecord.referents[stopsTuple] = set() # Put empty list here to not break later.
            else:
                for pointOnLink in closestLinks:
                    linkID = pointOnLink.link.id
                    # Use the ID here so we have a common reference among all subset networks.
                    
                    # Count that this link is matched to this stop.
                    if linkID not in stopResolveRecord.linkCounts:
                        stopResolveRecord.linkCounts[linkID] = 0
                        scores[linkID] = pointOnLink.refDist # This should be the same among all subsets. 
                    stopResolveRecord.linkCounts[linkID] += len(treeEntryIndices)
                    
                    # Identify where in the resultTree matched stops list this matched stop occurs.
                    if stopsTuple not in stopResolveRecord.referents:
                        stopResolveRecord.referents[stopsTuple] = set()
                    stopResolveRecord.referents[stopsTuple] |= treeEntryIndices
        
            # Set the counter for the total number of references.
            stopResolveRecord.refCount += len(allTripsBundles[stopsTuple].trips)
                
        del stopsTuple, treeEntryIndices, closestLinks, pointOnLink, linkID
                
        # Store the maximum score so that we can invert them to get the sorting right: 
        scoreMax = max(scores.values())
                
        # Build a list for identifying the best links by the refDist perpendicular distance metric:
        sortList = []
        for linkID, linkCount in stopResolveRecord.linkCounts.iteritems():
            sortList.append((linkCount, scoreMax - scores[linkID], linkID))
        
        # sortList will now have the links with the most common usage and smallest original scores at the bottom: 
        sortList = sorted(sortList)
        # TODO: Do we need to truncate the list to speed up processing?

        # This will hold the scores for all of the path tests: 
        allScores = {} # stops tuple -> linkID -> float
        "@type allScores: dict<tuple<int>, dict<int, float>>"                        
        
        linkID = stopsTuple = treeEntryIndices = tripsBundle = forceLinks = prevRestart = lfFlag = None
        while sortList:
            # Try out this link by invalidating the stop point and doing a refine cycle in each trip:
            linkID = sortList[-1][2]
            
            for stopsTuple, treeEntryIndices in stopResolveRecord.referents.iteritems():
                tripsBundle = allTripsBundles[stopsTuple]
                forceLinks = tripsBundle.forceLinks[:]
                if stopsTuple not in allScores:
                    allScores[stopsTuple] = {}

                prevRestart = {}
                for treeEntryIndex in treeEntryIndices:
                    # Invalidate the respective path match point for this stop and trip. The refine() call
                    # will then reevaluate those points along the path.
                    prevRestart[treeEntryIndex] = tripsBundle.resultTree[treeEntryIndex].restart 
                    tripsBundle.resultTree[treeEntryIndex].restart = True
                            
                    # Now, place all updated candidate links into the forceLinks list. Normally, there will
                    # just be one candidate, but if a route has any loops in it and this stop is hit twice, then
                    # we need a set of link objects that shares the underlying link.
                    forceLinks[treeEntryIndex] = vistaNetwork.linkMap[linkID]
                    
                pathEngine.setForceLinks(forceLinks)

                # Refine the path and see what the score is:
                lfFlag = True
                sys.stderr.write(".")
                resultTreeRefined = pathEngine.refinePath(tripsBundle.resultTree, vistaNetwork) 

                # Did the refine just flat-out fail? (e.g. do we still have a restart in there?)
                if sum([pathEnd.restart for pathEnd in resultTreeRefined]) == 0:
                    # No, continue recording (otherwise ignore this test)
                    allScores[stopsTuple][linkID] = resultTreeRefined[-1].totalCost
                    
                # Reset the invalidation.
                for treeEntryIndex in treeEntryIndices:
                    tripsBundle.resultTree[treeEntryIndex].restart = prevRestart[treeEntryIndex]
 
                del resultTreeRefined
            del sortList[-1]
        if lfFlag:
            print(file=sys.stderr)
        del linkID, stopsTuple, treeEntryIndices, tripsBundle, forceLinks, prevRestart, lfFlag
            
        # Now normalize the link scores for each stops tuple and multiply by the reference count:
        stopsTuple = scores = score = scoreMax = linkID = None
        scoreSums = {} # linkID -> float
        "@type scoreSums: dict<int, float>"
        for stopsTuple, scores in allScores.iteritems():
            if scores:
                scoreMax = max(scores.values())
                for linkID in scores.iterkeys():
                    scores[linkID] = scoreMax - scores[linkID] # Flip it so that the high score is the best.
                    scores[linkID] /= scoreMax if scoreMax > 0.0 else 1.0
                    if linkID not in scoreSums:
                        scoreSums[linkID] = 0.0
                    scoreSums[linkID] += scores[linkID] * len(stopRecord.referents[stopsTuple]) # One for each trip.
        
        # Use the scores for each link to rank the links:
        sortList = []
        for linkID, score in scoreSums.iteritems():
            sortList.append((score, linkID))
            
        sortList = sorted(sortList)
        # Now the best-scored link is at the end of the list.
        del stopsTuple, scores, score, scoreMax, linkID

        # Assign links to stops in order of descending popularity until all references are exhausted:
        stopsTuple = treeEntryIndices = tripsBundle = trip = treeEntryIndex = resultTreeRefined = prevRestart = None
        usedTripIDs = set();
        "@type usedTripIDs: set<int>"
        usedTrips = set();
        "@type usedTrips: set<tuple<<int>>"
        referenceCnt = 0

        firstLink = True        
        while sortList and referenceCnt < stopResolveRecord.refCount:
            # Find trips that use the linkID
            for stopsTuple, treeEntryIndices in stopResolveRecord.referents.iteritems():
                if stopsTuple not in usedTrips and sortList[-1][1] in allScores[stopsTuple]:
                    # This highest scoring linkID is in this trip. Invalidate the corresponding point in the tree list.
                    tripsBundle = allTripsBundles[stopsTuple] 
                    
                    for treeEntryIndex in treeEntryIndices:
                        # Skip it if we already use the correct linkID:
                        if tripsBundle.resultTree[treeEntryIndex].pointOnLink.link.id != sortList[-1][1]:
                            prevRestart = tripsBundle.resultTree[treeEntryIndex].restart;
                            tripsBundle.resultTree[treeEntryIndex].restart = True;
                            
                            # Modify the original forceLinks to keep all previous forced links in place. 
                            tripsBundle.forceLinks[treeEntryIndex] = vistaNetwork.linkMap[sortList[-1][1]]
                            
                            pathEngine.setForceLinks(tripsBundle.forceLinks)
                            
                            # Refine the path:
                            resultTreeRefined = pathEngine.refinePath(tripsBundle.resultTree, vistaNetwork)
                        else:
                            # Shortcut the refine process in cases where the link is the current one.
                            resultTreeRefined = tripsBundle.resultTree
    
                        # Did the refine fail? (Or is it unchanged?)
                        if tripsBundle.resultTree[treeEntryIndex].pointOnLink.link.id == sortList[-1][1] or sum([pathEnd.restart for pathEnd in resultTreeRefined]) == 0:
                            # No-- save the result.
                            tripsBundle.resultTree = resultTreeRefined
                            usedTrips.add(stopsTuple)
                            for trip in tripsBundle.trips:
                                usedTripIDs.add(trip.tripID)
                            referenceCnt += len(tripsBundle.trips)
                            
                            # Check to see if this is being saved to another linkID than the first one.
                            if not firstLink:
                                print("WARNING: The final match of Link %d for Shape %s is a different link."
                                    % (sortList[-1][1], tripsBundle.label), file=sys.stderr)
                        else:
                            print("WARNING: The match of Link %d for Shape %s failed."
                                % (sortList[-1][1], tripsBundle.label), file=sys.stderr)
                            tripsBundle.resultTree[treeEntryIndex].restart = prevRestart
                                
            if usedTripIDs:
                if firstLink:
                    print("INFO: %d trip(s) were matched to the best-scoring link %d." % (len(usedTripIDs), sortList[-1][1]), file=sys.stderr)
                firstLink = False
            
            del sortList[-1]
        if not sortList and referenceCnt < stopResolveRecord.refCount:
            # We have run out of candidates too early:
            tripIDList = ""
            tripIDCount = 0
            for stopsTuple in stopResolveRecord.referents.iterkeys(): 
                for trip in allTripsBundles[stopsTuple].trips:
                    if trip.tripID not in usedTripIDs:
                        if len(tripIDList) > 0:
                            tripIDList += ", "
                        tripIDList += str(trip.tripID)
                        tripIDCount += 1
            print("WARNING: %d trip(s) could not be matched; leaving unchanged: (tripID(s) %s)." % (tripIDCount, tripIDList), file=sys.stderr)
            del tripIDList, tripIDCount
        del stopsTuple, usedTripIDs, usedTrips, firstLink, trip, treeEntryIndices, tripsBundle, treeEntryIndex, resultTreeRefined, prevRestart, referenceCnt
            
    problemReportNodes = {}
    "@type problemReportNodes: dict<tuple<int>, path_engine.PathEnd>"
    
    print("INFO: -- End stop --", file=sys.stderr)
    print("INFO: Performing final output...", file=sys.stderr)
    for shapeID, stopsTuples in shapeStops.iteritems():
        for stopsTuple, tripsBundle in stopsTuples.iteritems():
            # Deal with Problem Report:
            # TODO: The Problem Report will include all nodes on each path regardless of valid time interval;
            # However; we will not have gotten here if the trip was entirely outside of it. 
            if problemReport:
                problemReportNodes[stopsTuple] = assembleProblemReport(tripsBundle.resultTree, vistaNetwork)
                    
            # Walk through our output link list and see when in time the resultTree entries occur. Keep those
            # that fall within our given time interval and entirely bail out on this trip if we are entirely
            # outside of the time range. We do this here because of the possibility that a route is shortened
            # because we are trying to match to, say, a subnetwork of a regional network. We had to have done
            # the steps above in order to know this.
            for trip in tripsBundle.trips:
                stopMatches = []
                "@type stopMatches: list<path_engine.PathEnd>"
                
                # Make a stop times lookup for this trip: shapeSeq -> StopTimesEntry
                stopTimesLookup = {}
                "@type stopTimes: dict<int, gtfs.StopTimesEntry>"
                for stopTimesEntry in gtfsStopTimes[trip]:
                    stopTimesLookup[stopTimesEntry.stopSeq] = stopTimesEntry
                    # This may not be necessary if stop sequences are always starting at 0 and sequential, but if they aren't,
                    # then this dictionary will always work. 
                
                rejectFlag = False
                for treeEntry in tripsBundle.resultTree:
                    "@type treeEntry: path_engine.PathEnd"
                    if treeEntry.shapeEntry.typeID == 1:
                        gtfsStopTime = stopTimesLookup[treeEntry.shapeEntry.shapeSeq]
                        if excludeBegin and gtfsStopTime.arrivalTime < startTime or excludeEnd and gtfsStopTime.arrivalTime > endTime:
                            # Throw away this entire route because it is excluded and part of it falls outside:
                            print("INFO: Trip %d excluded because activity happens outside of the valid time range.", trip.tripID, file=sys.stderr)
                            del stopMatches[:]
                            rejectFlag = True
                            break
                        elif (widenBegin or gtfsStopTime.arrivalTime >= startTime) and (widenEnd or gtfsStopTime.arrivalTime <= endTime):
                            stopMatches.append(treeEntry)
                        
                # Then, output the results if we had not been rejected:
                foundStopSet = set()
                if not rejectFlag:
                    if len(stopTimesLookup) > 0 and len(stopMatches) == 0:
                        # TODO: Because of a continue further above, this should never happen. 
                        print("INFO: No stops fall within the valid time range.", file=sys.stderr)
                    outSeqCtr = tripsBundle.longestStart
                    minTime = warmupStartTime
                    maxTime = cooldownEndTime
                    foundValidStop = False
                    stopMatchIndex = 0
                    for treeEntry in tripsBundle.resultTree:
                        # First, output the links leading up to this stop:
                        if len(treeEntry.routeInfo) - 1 > 0:
                            for routeInfoElem in treeEntry.routeInfo[0:-1]:
                                print('"%d","%d","%d",,,' % (trip.tripID, outSeqCtr, routeInfoElem.id), file=outFile)
                                outSeqCtr += 1
                            
                        if stopMatchIndex < len(stopMatches) and treeEntry == stopMatches[stopMatchIndex]:
                            foundStopSet.add(treeEntry.shapeEntry.shapeSeq) # Check off this stop sequence.
                            foundValidStop = True
                            stopID = stopTimesLookup[treeEntry.shapeEntry.shapeSeq].stop.stopID
                            print('"%d","%d","%d","%d","%d",' % (trip.tripID, outSeqCtr, treeEntry.pointOnLink.link.id,
                                stopID, DWELLTIME_DEFAULT), file=outFile)
                            if stopID in ret and ret[stopID].link.id != treeEntry.pointOnLink.link.id:
                                print("WARNING: stopID %d is attempted to be assigned to linkID %d, but it had already been assigned to linkID %d." \
                                    % (stopID, treeEntry.pointOnLink.link.id, ret[stopID].link.id), file=sys.stderr)
                                # TODO: This is a tricky problem. This means that among multiple bus routes, the same stop had been
                                # found to best fit two different links. I don't exactly know the best way to resolve this, other
                                # than (for NMC analyses) to create a "fake" stop that's tied with the new link. 
                            else:
                                ret[stopID] = treeEntry.pointOnLink
                                
                            # Check on the minimum/maximum time range:
                            gtfsStopTime = stopTimesLookup[treeEntry.shapeEntry.shapeSeq]
                            minTime = min(gtfsStopTime.arrivalTime, minTime)
                            maxTime = max(gtfsStopTime.arrivalTime, maxTime)
                            stopMatchIndex += 1
                        else:
                            # The linkID has nothing to do with any points in consideration.  Report it without a stop:
                            if foundValidStop or not excludeUpstream:
                                print('"%d","%d","%d",,,' % (trip.tripID, outSeqCtr, treeEntry.pointOnLink.link.id), file=outFile)
                        outSeqCtr += 1
                        # TODO: For start time estimation (as reported in the public.bus_frequency.csv output), it may be
                        # ideal to keep track of linear distance traveled before the first valid stop.
                        
                    # Widen out the valid interval if needed:
                    warmupStartTime = min(minTime, warmupStartTime)
                    cooldownEndTime = max(maxTime, cooldownEndTime)
        
                # Are there any stops left over?  If so, report them to say that they aren't in the output file.
                stopTimes = gtfsStopTimes[gtfsTrips[trip.tripID]]
                "@type stopTimes: list<gtfs.StopTimesEntry>"
                startGap = -1
                endGap = -1
                for gtfsStopTime in stopTimes:
                    "@type gtfsStopTime: gtfs.StopTimesEntry"
                    flag = False
                    if gtfsStopTime.stopSeq not in foundStopSet:
                        # This stop is unaccounted for:
                        if startGap < 0:
                            startGap = gtfsStopTime.stopSeq
                        endGap = gtfsStopTime.stopSeq
                        
                        # Old message is very annoying, especially if the underlying topology is a subset of shapefile
                        # geographic area and there's a ton of them. That's why there is the new range message as shown below.
                        # print("WARNING: Trip tripID %d, stopID %d stop seq. %d will not be in the bus_route_link file." % (tripID,
                        #    gtfsStopTime.stop.stopID, gtfsStopTime.stopSeq), file=sys.stderr)
                        
                        if problemReport:
                            revisedNodeList = problemReportNodes[gtfsTrips[tripID].shapeEntries[0].shapeID]  
                            if gtfsStopTime.stopSeq not in revisedNodeList:
                                # Make a dummy "error" node for reporting.
                                newShape = gtfs.ShapesEntry(gtfsTrips[tripID].shapeEntries[0].shapeID,
                                    gtfsStopTime.stopSeq, gtfsStopTime.stop.gpsLat,gtfsStopTime.stop.gpsLng, False)
                                newPointOnLink = graph.PointOnLink(None, 0)
                                newPointOnLink.pointX = gtfsStopTime.stop.pointX
                                newPointOnLink.pointY = gtfsStopTime.stop.pointY
                                newNode = path_engine.PathEnd(newShape, newPointOnLink)
                                newNode.restart = True
                                revisedNodeList[gtfsStopTime.stopSeq] = newNode
                    else:
                        flag = True
                    if (flag or gtfsStopTime.stopSeq == stopTimes[-1].stopSeq) and startGap >= 0:
                        subStr = "seqs %d-%d" % (startGap, endGap) if startGap != endGap else "seq %d" % startGap
                        print("WARNING: Because of time exclusion or subset network, tripID %d, stop %s will not be in the bus_route_link file." % (trip.tripID, subStr),
                            file=sys.stderr)
                        startGap = -1

    # Deal with Problem Report:
    if problemReport:
        print("INFO: Output problem report CSV...", file=sys.stderr)
        problemReportNodesOut = {}
        for stopsTuple, prNode in problemReportNodes.iteritems():
            seqs = prNode.keys()
            seqs.sort()
            ourTgtList = []
            for seq in seqs:
                ourTgtList.append(prNode[seq])
            problemReportNodesOut[stopsTuple] = ourTgtList                
        problem_report.problemReport(problemReportNodesOut, vistaNetwork)
    
    return ret, warmupStartTime, cooldownEndTime 

def dumpBusStops(gtfsStops, stopLinkMap, userName, networkName, outFile = sys.stdout):
    """
    dumpBusRouteStops dumps out a public.bus_route_link.csv file contents.
    @type gtfsStops: dict<int, StopsEntry>
    @type stopLinkMap: dict<int, graph.PointOnLink>
    @type userName: str
    @type networkName: str
    @type outFile: file
    """
    _outHeader("public.bus_stop", userName, networkName, outFile)
    print('"id","link","name","location",', file = outFile)
    
    # Iterate through the stopLinkMap:
    for stopID in stopLinkMap:
        "@type stopID: int"
        pointOnLink = stopLinkMap[stopID]
        "@type pointOnLink: graph.PointOnLink"
        print('"%d","%d","%s","%d"' % (stopID, pointOnLink.link.id, gtfsStops[stopID].stopName, int(pointOnLink.dist)), file = outFile) 

def main(argv):
    global problemReport
    excludeUpstream = False
    
    # Initialize from command-line parameters:
    if len(argv) < 7:
        syntax(1)
    dbServer = argv[1]
    networkName = argv[2]
    userName = argv[3]
    password = argv[4]
    shapePath = argv[5]
    pathMatchFilename = argv[6]
    endTimeInt = 86400
    refTime = None
    widenBegin = False
    widenEnd = False
    excludeBegin = False
    excludeEnd = False
    
    restrictService = set()
    "@type restrictService: set<string>"

    if len(argv) > 6:
        i = 7
        while i < len(argv):
            if argv[i] == "-t" and i < len(argv) - 1:
                refTime = datetime.strptime(argv[i + 1], '%H:%M:%S')
                i += 1
            elif argv[i] == "-e" and i < len(argv) - 1:
                endTimeInt = int(argv[i + 1])
                i += 1
            elif argv[i] == "-c" and i < len(argv) - 1:
                restrictService.add(argv[i + 1])
                i += 1
            elif argv[i] == "-u":
                excludeUpstream = True
            elif argv[i] == "-w":
                widenBegin = True
                widenEnd = True
            elif argv[i] == "-wb":
                widenBegin = True
            elif argv[i] == "-we":
                widenEnd = True
            elif argv[i] == "-x":
                excludeBegin = True
                excludeEnd = True
            elif argv[i] == "-xb":
                excludeBegin = True
            elif argv[i] == "-xe":
                excludeEnd = True
            elif argv[i] == "-p":
                problemReport = True
            i += 1
    
    if refTime is None:
        print("ERROR: No reference time is specified. You must use the -t parameter.", file=sys.stderr)
        syntax(1)
    endTime = refTime + timedelta(seconds = endTimeInt)
    
    if widenBegin and excludeBegin:
        print("ERROR: Widening (-w or -wb) and exclusion (-x or -xb) cannot be used together.")
        syntax(1)    
    if widenEnd and excludeEnd:
        print("ERROR: Widening (-w or -we) and exclusion (-x or -xe) cannot be used together.")
        syntax(1)
    
    # Restore the stuff that was built with path_match:
    (vistaGraph, gtfsShapes, gtfsNodes, unusedShapeIDs) = restorePathMatch(dbServer, networkName, userName,
        password, shapePath, pathMatchFilename)
    
    # Read in the routes information:
    print("INFO: Read GTFS routesfile...", file=sys.stderr)
    gtfsRoutes = gtfs.fillRoutes(shapePath)
    "@type gtfsRoutes: dict<int, RoutesEntry>"
    
    # Read in the stops information:
    print("INFO: Read GTFS stopsfile...", file=sys.stderr)
    gtfsStops = gtfs.fillStops(shapePath, vistaGraph.gps)
    "@type gtfsStops: dict<int, StopsEntry>"
    
    # Read in the trips information:
    print("INFO: Read GTFS tripsfile...", file=sys.stderr)
    (gtfsTrips, unusedTripIDs) = gtfs.fillTrips(shapePath, gtfsShapes, gtfsRoutes, unusedShapeIDs, restrictService)
    "@type gtfsTrips: dict<int, TripsEntry>"
    "@type unusedTripIDs: set<int>"
        
    # Read stop times information:
    print("INFO: Read GTFS stop times...", file=sys.stderr)
    gtfsStopTimes = gtfs.fillStopTimes(shapePath, gtfsTrips, gtfsStops, unusedTripIDs)
    "@type gtfsStopTimes: dict<TripsEntry, list<StopTimesEntry>>"
        
    # Output the routes_link file:
    print("INFO: Dumping public.bus_route_link.csv...", file=sys.stderr)
    with open("public.bus_route_link.csv", 'w') as outFile:
        (stopLinkMap, newStartTime, newEndTime) = dumpBusRouteLinks(gtfsTrips, gtfsStops, gtfsStopTimes, gtfsNodes, vistaGraph,
            STOP_SEARCH_RADIUS, excludeUpstream, userName, networkName, refTime, endTime, widenBegin, widenEnd,
            excludeBegin, excludeEnd, outFile)
        "@type stopLinkMap: dict<int, graph.PointOnLink>"
    
    # Filter only to bus stops and stop times that are used in the routes_link output:
    gtfsStopsFilterList = [gtfsStopID for gtfsStopID in gtfsStops if gtfsStopID not in stopLinkMap]
    for gtfsStopID in gtfsStopsFilterList:
        del gtfsStops[gtfsStopID]
    del gtfsStopsFilterList
    
    # Then, output the output the stop file:
    print("INFO: Dumping public.bus_stop.csv...", file=sys.stderr)
    with open("public.bus_stop.csv", 'w') as outFile:
        dumpBusStops(gtfsStops, stopLinkMap, userName, networkName, outFile)
        
    print("INFO: Dumping public.bus_frequency.csv...", file=sys.stderr)
    validTrips = {}
    "@type validTrips: dict<int, gtfs.TripsEntry>"
    with open("public.bus_frequency.csv", 'w') as outFile:
        _outHeader("public.bus_frequency", userName, networkName, outFile)
        print("\"route\",\"period\",\"frequency\",\"offsettime\",\"preemption\"", file = outFile)
        
        # Okay, here we iterate through stops until we get to the first defined one. That will
        # then affect the offsettime. (This is needed because of the idea that we want to start
        # a bus in the simulation wrt the topology that supports it, skipping those stops that
        # may fall outside the topology.)
        totalCycle = int((newEndTime - newStartTime).total_seconds()) 
        tripIDs = gtfsTrips.keys()
        tripIDs.sort()
        for tripID in tripIDs:
            stopsEntries = gtfsStopTimes[gtfsTrips[tripID]]
            for gtfsStopTime in stopsEntries:
                if gtfsStopTime.stop.stopID in gtfsStops:
                    # Here is a first valid entry! Use this offset value.
                    
                    # TODO: This could be inaccurate because the offset is that of the first
                    # valid stop time encountered in the underlying topology, not approximated
                    # to the first valid link encountered. While this isn't a big deal for an
                    # area with a high stop density, it could be a problem for limited-stop
                    # service where there happens to be a low density around where the bus
                    # first appears in the underlying topology.
                    stopTime = gtfsStopTime.arrivalTime
                    
                    # Adjust for cases where we need to add a day.
                    if stopTime < newStartTime: # Assume that we're working just within a day.
                        stopTime += timedelta(days = int((newStartTime - stopTime).total_seconds()) / 86400 + 1)
                    print("%d,1,%d,%d,0" % (tripID, totalCycle, int((stopTime - newStartTime).total_seconds())),
                        file = outFile)
                    validTrips[tripID] = gtfsTrips[tripID] # Record as valid.
                    break
                    
                # A byproduct of this scheme is that no bus_frequency entry will appear for
                # routes that don't have stops in the underlying topology.

    # Output the routes file:
    print("INFO: Dumping public.bus_route.csv...", file=sys.stderr)
    with open("public.bus_route.csv", 'w') as outFile:
        dumpBusRoutes(validTrips, userName, networkName, outFile)

    # Finally, define one period that spans the whole working time, which all of the individually
    # defined routes (again, one route per trip) will operate in.
    print("INFO: Dumping public.bus_period.csv...", file=sys.stderr)
    with open("public.bus_period.csv", 'w') as outFile:
        _outHeader("public.bus_period", userName, networkName, outFile)
        print("\"id\",\"starttime\",\"endtime\"", file = outFile)
        # The start time printed here is relative to the reference time.
        print("1,0,%d" % endTimeInt, file = outFile)
        
    if widenBegin or widenEnd:
        # Report the implicit adjustment in times because of warmup or cooldown:
        startTimeDiff = refTime - newStartTime
        endTimeDiff = newEndTime - endTime
        print("INFO: Widening requires start %d sec. earlier and duration %d sec. longer." % (startTimeDiff.total_seconds(),
            endTimeDiff.total_seconds() + startTimeDiff.total_seconds()), file=sys.stderr)
        totalTimeDiff = newEndTime - newStartTime
        print("INFO: New time reference is %s, duration %d sec." % (newStartTime.strftime("%H:%M:%S"), totalTimeDiff.total_seconds()),
            file=sys.stderr)

    print("INFO: Done.", file=sys.stderr)

# Boostrap:
if __name__ == '__main__':
    main(sys.argv)
