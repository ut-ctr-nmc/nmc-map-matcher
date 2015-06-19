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
"@var DWELLTIME_DEFAULT is the dwell time to report in the bus_route_link.csv file output."

problemReport = False
"@var problemReport is set to true when the -p parameter is specified."

def syntax():
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
    sys.exit(0)

def restorePathMatch(dbServer, networkName, userName, password, shapePath, pathMatchFilename):
    # Get the database connected:
    print("INFO: Connect to database...", file = sys.stderr)
    database = vista_network.connect(dbServer, userName, password, networkName)
    
    # Read in the topology from the VISTA database:
    print("INFO: Read topology from database...", file = sys.stderr)
    vistaGraph = vista_network.fillGraph(database)
    
    # Read in the shapefile information:
    print("INFO: Read GTFS shapefile...", file = sys.stderr)
    gtfsShapes = gtfs.fillShapes(shapePath, vistaGraph.gps)

    # Read the path-match file:
    print("INFO: Read the path-match file '%s'..." % pathMatchFilename, file = sys.stderr)
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

class StopMatch:
    """
    StopMatch represents a matched stop, used within dumpBusRouteLinks.
    
    @type bestTreeEntry: path_engine.PathEnd
    @type matchCtr: int
    @type linkID: int
    """
    def __init__(self, linkID):
        """
        Sets member variables to zero and initialized.
        @type linkID: int
        """
        self.bestTreeEntry = None
        self.matchCtr = 0
        self.linkID = linkID

def dumpBusRouteLinks(gtfsTrips, gtfsStopTimes, gtfsNodes, vistaNetwork, stopSearchRadius, excludeUpstream, userName,
        networkName, startTime, endTime, widenBegin, widenEnd, excludeBegin, excludeEnd, outFile = sys.stdout):
    """
    dumpBusRouteLinks dumps out a public.bus_route_link.csv file contents. This also will remove all stop times and trips
    that fall outside of the valid evaluation interval as dictated by the exclusion parameters.
    @type gtfsTrips: dict<int, gtfs.TripsEntry>
    @type gtfsStopTimes: dict<TripsEntry, list<StopTimesEntry>>
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
                                        stopSearchRadius, 1, 1, 1, sys.maxint, sys.maxint)
    pathEngine.limitClosestPoints = 8
    pathEngine.limitSimultaneousPaths = 6
    pathEngine.maxHops = 12
    pathEngine.logFile = None # Suppress the log outputs for the path engine; enough stuff will come from other sources.

    problemReportNodes = {}
    "@type problemReportNodes: dict<?, path_engine.PathEnd>"
    
    tripIDs = gtfsTrips.keys()
    tripIDs.sort()
    for tripID in tripIDs:
        if gtfsTrips[tripID].shapeEntries[0].shapeID not in gtfsNodes:
            # This happens if the incoming files contain a subset of all available topology.
            print("WARNING: Skipping route for trip %d because no points are available." % tripID, file = sys.stderr)
            continue
        
        treeNodes = gtfsNodes[gtfsTrips[tripID].shapeEntries[0].shapeID]
        "@type treeNodes: list<path_engine.PathEnd>"
        
        # Step 1: Find the longest distance of contiguous valid links within the shape for each trip:
        startIndex = -1
        curIndex = 0
        linkCount = 0
        totalLinks = 0
        
        longestStart = -1
        longestEnd = len(treeNodes)
        longestDist = sys.float_info.min
        longestLinkCount = 0
        
        while curIndex <= len(treeNodes):
            if (curIndex == len(treeNodes)) or (curIndex == 0) or treeNodes[curIndex].restart:
                totalLinks += 1
                linkCount += 1
                if (curIndex > startIndex) and (startIndex >= 0):
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
            if (longestStart > 0) or (longestEnd < len(treeNodes)):
                print("WARNING: For shape ID %s from seq. %d through %d, %.2g%% of %d links will be used." \
                      % (str(treeNodes[longestStart].shapeEntry.shapeID), treeNodes[longestStart].shapeEntry.shapeSeq,
                         treeNodes[longestEnd - 1].shapeEntry.shapeSeq, 100 * float(longestLinkCount) / float(totalLinks),
                         totalLinks), file = sys.stderr)
            
            # Step 2: Ignore routes that are entirely outside our valid time interval.
            flag = False
            if len(gtfsStopTimes[gtfsTrips[tripID]]) == 0:
                # This will happen if we don't have stops defined. In this case, we want to go ahead and process the bus_route_link
                # outputs because we don't know if the trip falls in or out of the valid time range.
                flag = True
            else:
                for stopEntry in gtfsStopTimes[gtfsTrips[tripID]]:
                    if stopEntry.arrivalTime >= startTime and stopEntry.arrivalTime <= endTime:
                        flag = True
                        break
            if not flag:
                # This will be done silently because (depending upon the valid interval) there could be
                # hundreds of these in a GTFS set.
                continue
            
            # Step 3: Match up stops to that contiguous list:
            # At this point, we're doing something with this.
            print("INFO: -- Matching stops for trip %d --" % tripID, file = sys.stderr)
        
            stopTimes = gtfsStopTimes[gtfsTrips[tripID]]
            "@type stopTimes: list<gtfs.StopTimesEntry>"
            
            # Isolate the relevant VISTA tree nodes: (Assume from above that this is a non-zero length array)
            ourGTFSNodes = treeNodes[longestStart:longestEnd]
            
            # We are going to recreate a small VISTA network from ourGTFSNodes and then match up the stops to that.
            # First, prepare the small VISTA network:
            vistaSubset = graph.GraphLib(vistaNetwork.gps.latCtr, vistaNetwork.gps.lngCtr)
            vistaNodePrior = None
            "@type vistaNodePrior: graph.GraphNode"
            
            # Build a list of links:
            outLinkIDList = []
            "@type outLinkList: list<int>"
            
            # Plop in the start node:
            vistaNodePrior = graph.GraphNode(ourGTFSNodes[0].pointOnLink.link.origNode.id,
                ourGTFSNodes[0].pointOnLink.link.origNode.gpsLat, ourGTFSNodes[0].pointOnLink.link.origNode.gpsLng)
            vistaSubset.addNode(vistaNodePrior)
            outLinkIDList.append(ourGTFSNodes[0].pointOnLink.link.id)
            
            # Link together nodes as we traverse through them:
            for ourGTFSNode in ourGTFSNodes:
                "@type ourGTFSNode: path_engine.PathEnd"
                # There should only be one destination link per VISTA node because this comes form our tree.
                # If there is no link or we're repeating the first one, then there were no new links assigned.
                if (len(ourGTFSNode.routeInfo) < 1) or ((len(outLinkIDList) == 1) \
                        and (ourGTFSNode.routeInfo[0].id == ourGTFSNodes[0].pointOnLink.link.id)):
                    continue
                for link in ourGTFSNode.routeInfo:
                    "@type link: graph.GraphLink"
                
                    if link.id not in vistaNetwork.linkMap:
                        print("WARNING: In finding bus route links, link ID %d is not found in the VISTA network." % link.id, file = sys.stderr)
                        continue
                    origVistaLink = vistaNetwork.linkMap[link.id]
                    "@type origVistaLink: graph.GraphLink"
                    
                    if origVistaLink.origNode.id not in vistaSubset.nodeMap:
                        # Create a new node:
                        vistaNode = graph.GraphNode(origVistaLink.origNode.id, origVistaLink.origNode.gpsLat, origVistaLink.origNode.gpsLng)
                        vistaSubset.addNode(vistaNode)
                    else:
                        # The path evidently crosses over itself.  Reuse an existing node.
                        vistaNode = vistaSubset.nodeMap[origVistaLink.origNode.id]
                        
                    # We shall label our links as indices into the stage we're at in ourGTFSNodes links.  This will allow for access later.
                    if outLinkIDList[-1] not in vistaSubset.linkMap:
                        vistaSubset.addLink(graph.GraphLink(outLinkIDList[-1], vistaNodePrior, vistaNode))
                    vistaNodePrior = vistaNode
                    outLinkIDList.append(link.id)
                    
            # And then finish off the graph with the last link:
            if ourGTFSNode.pointOnLink.link.destNode.id not in vistaSubset.nodeMap:
                vistaNode = graph.GraphNode(ourGTFSNode.pointOnLink.link.destNode.id, ourGTFSNode.pointOnLink.link.destNode.gpsLat, ourGTFSNode.pointOnLink.link.destNode.gpsLng)
                vistaSubset.addNode(vistaNode)
            if outLinkIDList[-1] not in vistaSubset.linkMap:
                vistaSubset.addLink(graph.GraphLink(outLinkIDList[-1], vistaNodePrior, vistaNode))
            
            # Then, prepare the stops as GTFS shapes entries:
            print("INFO: Mapping stops to VISTA network...", file = sys.stderr)
            gtfsShapes = []
            gtfsStopsLookup = {}
            "@type gtfsStopsLookup: dict<int, gtfs.StopTimesEntry>"
            
            # Append an initial dummy shape to force routing through the path start:
            gtfsShapes.append(gtfs.ShapesEntry(-1, -1, ourGTFSNodes[0].pointOnLink.link.origNode.gpsLat,
                                                ourGTFSNodes[0].pointOnLink.link.origNode.gpsLng))
            
            # Append all of the stops:
            for gtfsStopTime in stopTimes:
                "@type gtfsStopTime: gtfs.StopTimesEntry"
                gtfsShapes.append(gtfs.ShapesEntry(-1, gtfsStopTime.stopSeq, gtfsStopTime.stop.gpsLat, gtfsStopTime.stop.gpsLng))
                gtfsStopsLookup[gtfsStopTime.stopSeq] = gtfsStopTime

            # Append a trailing dummy shape to force routing through the path end:
            gtfsShapes.append(gtfs.ShapesEntry(-1, -1, ourGTFSNodes[-1].pointOnLink.link.destNode.gpsLat,
                                                ourGTFSNodes[-1].pointOnLink.link.destNode.gpsLng))
        
            # Find a path through our prepared node map subset:
            resultTree = pathEngine.constructPath(gtfsShapes, vistaSubset)
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
                revisedNodeList = {}
                prevNode = None
                "@type revisedNodeList = list<path_engine.PathEnd>"
                for stopNode in resultTree:
                    # Reconstruct a tree node in terms of the original network.
                    newShape = gtfs.ShapesEntry(gtfsTrips[tripID].shapeEntries[0].shapeID,
                        stopNode.shapeEntry.shapeSeq, stopNode.shapeEntry.lat, stopNode.shapeEntry.lng, False)
                    origLink = vistaNetwork.linkMap[stopNode.pointOnLink.link.id] 
                    newPointOnLink = graph.PointOnLink(origLink, stopNode.pointOnLink.dist,
                        stopNode.pointOnLink.nonPerpPenalty, stopNode.pointOnLink.refDist)
                    newNode = path_engine.PathEnd(newShape, newPointOnLink)
                    newNode.restart = False
                    newNode.totalCost = stopNode.totalCost
                    newNode.totalDist = stopNode.totalDist
                    newNode.routeInfo = []
                    for link in stopNode.routeInfo:
                        newNode.routeInfo.append(vistaNetwork.linkMap[link.id])
                    newNode.prevTreeNode = prevNode
                    prevNode = newNode
                    revisedNodeList[stopNode.shapeEntry.shapeSeq] = newNode
                problemReportNodes[gtfsTrips[tripID].shapeEntries[0].shapeID] = revisedNodeList 
        
            # Walk through our output link list and see where the resultTree entries occur:
            resultIndex = 0
            stopMatches = []
            "@type stopMatches: list<StopMatch>"
            rejectFlag = False
            for linkID in outLinkIDList:
                curResultIndex = resultIndex
                # This routine will advance resultIndex only if a stop is found for linkID, and will exit out when
                # no more stops are found for linkID.
                stopMatch = StopMatch(linkID)
                "@type stopMatch: StopMatch" 
                stopMatches.append(stopMatch)
                while curResultIndex < len(resultTree):
                    if resultTree[curResultIndex].pointOnLink.link.id == linkID:
                        # Only pay attention to this stop if it is within the valid time range:
                        gtfsStopTime = gtfsStopsLookup[resultTree[resultIndex].shapeEntry.shapeSeq]
                        if excludeBegin and gtfsStopTime.arrivalTime < startTime or excludeEnd and gtfsStopTime.arrivalTime > endTime:
                            # Throw away this entire route because it is excluded and part of it falls outside:
                            print("INFO: Excluded because of activity outside of the valid time range.", file = sys.stderr)
                            del stopMatches[:]
                            rejectFlag = True
                            break
                        elif (widenBegin or gtfsStopTime.arrivalTime >= startTime) and (widenEnd or gtfsStopTime.arrivalTime <= endTime):
                            if (stopMatch.bestTreeEntry is None) \
                                    or (resultTree[resultIndex].pointOnLink.refDist < stopMatch.bestTreeEntry.pointOnLink.refDist):
                                # Log the best match:
                                stopMatch.bestTreeEntry = resultTree[resultIndex]
                            stopMatch.matchCtr += 1
                        resultIndex = curResultIndex + 1
                    curResultIndex += 1
                    if (stopMatch and stopMatch.matchCtr == 0) \
                            or ((curResultIndex < len(resultTree)) and (resultTree[resultIndex].pointOnLink.link.id == linkID)):
                        continue
                    # We have gotten to the end of matched link(s). 
                    break
                if rejectFlag:
                    break
            
            # Then, output the results out if we are supposed to.
            foundStopSet = set()
            "@type foundStopSet: set<int>"
            if not rejectFlag:
                outSeqCtr = longestStart
                minTime = warmupStartTime
                maxTime = cooldownEndTime
                foundValidStop = False
                for stopMatch in stopMatches:    
                    if stopMatch.matchCtr > 1:
                        # Report duplicates:
                        print("WARNING: %d stops have been matched for TripID %d, LinkID %d. Keeping Stop %d, Stop Seq %d" \
                            % (stopMatch.matchCtr, tripID, stopMatch.linkID, gtfsStopsLookup[stopMatch.bestTreeEntry.shapeEntry.shapeSeq].stop.stopID,
                            stopMatch.bestTreeEntry.shapeEntry.shapeSeq), file = sys.stderr)
                        # TODO: This is a problem because VISTA only allows one stop per link. So, the stop that is closest to
                        # the link is the one that is the winner and the rest are ignored. We don't yet do anything intelligent with dwell
                        # times, etc.
                    if stopMatch.matchCtr > 0:
                        # Report the best match:
                        foundStopSet.add(stopMatch.bestTreeEntry.shapeEntry.shapeSeq) # Check off this stop sequence.
                        foundValidStop = True
                        print('"%d","%d","%d","%d","%d",' % (tripID, outSeqCtr, stopMatch.linkID,
                            gtfsStopsLookup[stopMatch.bestTreeEntry.shapeEntry.shapeSeq].stop.stopID, DWELLTIME_DEFAULT), file = outFile)
                        if gtfsStopsLookup[stopMatch.bestTreeEntry.shapeEntry.shapeSeq].stop.stopID in ret \
                                and ret[gtfsStopsLookup[stopMatch.bestTreeEntry.shapeEntry.shapeSeq].stop.stopID].link.id \
                                    != stopMatch.bestTreeEntry.pointOnLink.link.id:
                            print("WARNING: stopID %d is attempted to be assigned to linkID %d, but it had already been assigned to linkID %d." \
                                % (gtfsStopsLookup[stopMatch.bestTreeEntry.shapeEntry.shapeSeq].stop.stopID, stopMatch.bestTreeEntry.pointOnLink.link.id,
                                   ret[gtfsStopsLookup[stopMatch.bestTreeEntry.shapeEntry.shapeSeq].stop.stopID].link.id), file = sys.stderr)
                            # TODO: This is a tricky problem. This means that among multiple bus routes, the same stop had been
                            # found to best fit two different links. I don't exactly know the best way to resolve this, other
                            # than (for NMC analyses) to create a "fake" stop that's tied with the new link. 
                        else:
                            ret[gtfsStopsLookup[stopMatch.bestTreeEntry.shapeEntry.shapeSeq].stop.stopID] = stopMatch.bestTreeEntry.pointOnLink
                            
                        # Check on the minimum/maximum time range:
                        gtfsStopTime = gtfsStopsLookup[stopMatch.bestTreeEntry.shapeEntry.shapeSeq]
                        minTime = min(gtfsStopTime.arrivalTime, minTime)
                        maxTime = max(gtfsStopTime.arrivalTime, maxTime)
                    else:
                        # The linkID has nothing to do with any points in consideration.  Report it without a stop:
                        if foundValidStop or not excludeUpstream:
                            print('"%d","%d","%d",,,' % (tripID, outSeqCtr, stopMatch.linkID), file = outFile)
                    outSeqCtr += 1
                    # TODO: For start time estimation (as reported in the public.bus_frequency.csv output), it may be
                    # ideal to keep track of linear distance traveled before the first valid stop.
                    
                # Widen out the valid interval if needed:
                warmupStartTime = min(minTime, warmupStartTime)
                cooldownEndTime = max(maxTime, cooldownEndTime)

            # Are there any stops left over?  If so, report them to say that they aren't in the output file.
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
                    #    gtfsStopTime.stop.stopID, gtfsStopTime.stopSeq), file = sys.stderr)
                    
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
                    subStr = "Seqs. %d-%d" % (startGap, endGap) if startGap != endGap else "Seq. %d" % startGap
                    print("WARNING: Trip ID %d, Stop %s will not be in the bus_route_link file." % (tripID, subStr),
                        file = sys.stderr)
                    startGap = -1
        else:
            print("WARNING: No links for tripID %d." % tripID, file = sys.stderr)

    # Deal with Problem Report:
    if problemReport:
        print("INFO: Output problem report CSV...", file = sys.stderr)
        problemReportNodesOut = {}
        for shapeID in problemReportNodes:
            seqs = problemReportNodes[shapeID].keys()
            seqs.sort()
            ourTgtList = []
            for seq in seqs:
                ourTgtList.append(problemReportNodes[shapeID][seq])
            problemReportNodesOut[shapeID] = ourTgtList                
        problem_report.problemReport(problemReportNodesOut, vistaNetwork)
    
    return (ret, warmupStartTime, cooldownEndTime) 

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
        syntax()
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
        print("ERROR: No reference time is specified. You must use the -t parameter.", file = sys.stderr)
        syntax(1)
    endTime = refTime + timedelta(seconds = endTimeInt)
    
    if widenBegin and excludeBegin:
        print("ERROR: Widening (-w or -wb) and exclusion (-x or -xb) cannot be used together.")
        syntax(1)    
    if widenEnd and excludeEnd:
        print("ERROR: Widening (-w or -we) and exclusion (-x or -xe) cannot be used together.")
        syntax(1)
    
    # Default parameters:
    stopSearchRadius = 800
    
    # Restore the stuff that was built with path_match:
    (vistaGraph, gtfsShapes, gtfsNodes, unusedShapeIDs) = restorePathMatch(dbServer, networkName, userName,
        password, shapePath, pathMatchFilename)
    
    # Read in the routes information:
    print("INFO: Read GTFS routesfile...", file = sys.stderr)
    gtfsRoutes = gtfs.fillRoutes(shapePath)
    "@type gtfsRoutes: dict<int, RoutesEntry>"
    
    # Read in the stops information:
    print("INFO: Read GTFS stopsfile...", file = sys.stderr)
    gtfsStops = gtfs.fillStops(shapePath, vistaGraph.gps)
    "@type gtfsStops: dict<int, StopsEntry>"
    
    # Read in the trips information:
    print("INFO: Read GTFS tripsfile...", file = sys.stderr)
    (gtfsTrips, unusedTripIDs) = gtfs.fillTrips(shapePath, gtfsShapes, gtfsRoutes, unusedShapeIDs, restrictService)
    "@type gtfsTrips: dict<int, TripsEntry>"
    "@type unusedTripIDs: set<int>"
        
    # Read stop times information:
    print("INFO: Read GTFS stop times...", file = sys.stderr)
    gtfsStopTimes = gtfs.fillStopTimes(shapePath, gtfsTrips, gtfsStops, unusedTripIDs)
    "@type gtfsStopTimes: dict<TripsEntry, list<StopTimesEntry>>"
        
    # Output the routes_link file:
    print("INFO: Dumping public.bus_route_link.csv...", file = sys.stderr)
    with open("public.bus_route_link.csv", 'w') as outFile:
        (stopLinkMap, newStartTime, newEndTime) = dumpBusRouteLinks(gtfsTrips, gtfsStopTimes, gtfsNodes, vistaGraph,
            stopSearchRadius, excludeUpstream, userName, networkName, refTime, endTime, widenBegin, widenEnd,
            excludeBegin, excludeEnd, outFile)
        "@type stopLinkMap: dict<int, graph.PointOnLink>"
    
    # Filter only to bus stops and stop times that are used in the routes_link output:
    gtfsStopsFilterList = [gtfsStopID for gtfsStopID in gtfsStops if gtfsStopID not in stopLinkMap]
    for gtfsStopID in gtfsStopsFilterList:
        del gtfsStops[gtfsStopID]
    del gtfsStopsFilterList
    
    # Then, output the output the stop file:
    print("INFO: Dumping public.bus_stop.csv...", file = sys.stderr)
    with open("public.bus_stop.csv", 'w') as outFile:
        dumpBusStops(gtfsStops, stopLinkMap, userName, networkName, outFile)
        
    print("INFO: Dumping public.bus_frequency.csv...", file = sys.stderr)
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
    print("INFO: Dumping public.bus_route.csv...", file = sys.stderr)
    with open("public.bus_route.csv", 'w') as outFile:
        dumpBusRoutes(validTrips, userName, networkName, outFile)

    # Finally, define one period that spans the whole working time, which all of the individually
    # defined routes (again, one route per trip) will operate in.
    print("INFO: Dumping public.bus_period.csv...", file = sys.stderr)
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
            endTimeDiff.total_seconds() + startTimeDiff.total_seconds()), file = sys.stderr)
        totalTimeDiff = newEndTime - newStartTime
        print("INFO: New time reference is %s, duration %d sec." % (newStartTime.strftime("%H:%M:%S"), totalTimeDiff.total_seconds()),
            file = sys.stderr)

    print("INFO: Done.", file = sys.stderr)

# Boostrap:
if __name__ == '__main__':
    main(sys.argv)
