"""
path_engine.py contains logic for matching up GTFS paths to VISTA paths and the
    writing and reading of them.
@author: Kenneth Perrine
@contact: kperrine@utexas.edu
@organization: Network Modeling Center, Center for Transportation Research,
    Cockrell School of Engineering, The University of Texas at Austin 
@version: 1.0

@copyright: (C) 2026, The University of Texas at Austin
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

from collections.abc import Hashable, Iterable, Sequence
from typing import Final, NamedTuple
from nmc_mm_lib import graph
import operator, math, sys, copy
import logging

# Multiplier for shape-to-shape evaluations that happen while refining on a
# restart:
RESTART_PENALTY_MULT: Final[float] = 2.0

class PathEnd:
    """
    PathEnd is a single node used within the overall tree structure. This
    roughly equates to the "path_end" data structure outlined in Figure 2 of
    Perrine et al., 2015.
    """
    refPoint: graph.Map.Trackpoint # The geographic reference point
                                   # corresponding with this result
    pointOnLink: graph.Map.PointOnLink # Current point being considered
    totalCost: float # "s", the total score of the path represented
    prevTreeNode: PathEnd | None # "p", the previous PathEnd step in this path
    totalDist: float # Relates to total score of the path represented
    totalLinkCount: int # the number of links that had been traversed
    routeInfo: list[Hashable] # "l", a list of map links that have been
                              # traversed on the shortest path
    restart: bool # "r", a Boolean signifying a discontinuity

    def __init__(self,
                 refPoint: graph.Map.Trackpoint,
                 pointOnLink: graph.Map.PointOnLink):
        """
        Sets up values in this object, many of which need to be mutable
        """
        self.refPoint = refPoint
        self.pointOnLink = pointOnLink
        
        self.totalCost = 0.0
        self.prevTreeNode = None
        self.totalDist = 0.0
        self.totalLinkCount = 0
        self.routeInfo = []
        self.restart = False
        
    def cleanCopy(self) -> PathEnd:
        """
        cleanCopy initializes a new PathEnd object based on another one.
        """
        return PathEnd(self.refPoint, self.pointOnLink)

class PathEngine:
    """
    PathEngine contains constraints that guide the creation of a path.
    """
    class Params(NamedTuple):
        """
        Used for configuring the desired behavior of the path engine
        """
        searchRadius: float # TODO: Fill these in with the most common defaults.
        radiusPrimary: float
        radiusSecondary: float
        limitPathDist: float
        limitDirectDist: float
        limitDirectDistRev: float
        distFactor: float
        driftFactor: float
        nonPerpPenalty: float
        limitClosestPoints: int
        limitSimulPaths: int

    params: Params
    prevCosts: list[float] # A list of limitSimulPaths cost values that can be
                           # used to determine if proposed paths are worth
                           # traversing.
    maxHops: int = 12 # Limits number of nodes to be traversed in path-finding
    tossRatio: float = 1.0 # Disable the invalidation of short paths
    shapeScatterCache: list[graph.Map.PointOnLink] | None = None
    forceLinks: Sequence[Iterable[Hashable]] | None | None = None

    def __init__(self, params: Params):
        """
        Initializes the path engine with the given parameters; others are
        implicit from the defaults
        """
        self.params = params

    def scoreFunction(self,
                      prevGeoPoint: graph.Map.PointOnLink | None,
                      distance: float,
                      geoPoint: graph.Map.PointOnLink | None) -> float:
        """
        scoreFunction calculates a cost value given prior path distance, and deviation from the VISTA link.
        This corresponds with algorithm "ScoreFunction" in Perrine et al., 2015.

        @param prevGeoPoint: graph.PointOnLink
        @param distance: float
        @param geoPoint: graph.PointOnLink
        @return: Calculated score value
        """
        cost: float
        if prevGeoPoint is None:
            # We are starting anew.  Count the "black line distance" from the VISTA link to the GTFS point:
            if geoPoint is not None:
                cost = geoPoint.refDist * self.params.driftFactor
                if geoPoint.nonPerpPenalty:
                    cost = cost * self.params.nonPerpPenalty
            else:
                cost = 0.0
            return cost
        else:
            # We're jumping from one link to another, so add the "black line" distance to the total VISTA link distance:
            if geoPoint is not None:
                cost = geoPoint.refDist * self.params.driftFactor
                if geoPoint.nonPerpPenalty:
                    cost = cost * self.params.nonPerpPenalty
            else:
                cost = 0.0            
            return cost + abs(distance) * self.params.distFactor
            # Change from Perrine et al., 2015: Use absolute value of distance here because all movement
            # should be incrementing even in cases where a proposed path is moving back and forth on a link
            # because of shape point noise or tiny U-turns.
        
    def exceedsPreviousCosts(self, cost: float) -> bool:
        """
        Returns true if the given cost value exceeds the most expensive cost already recorded (if the list is
        limitSimulPaths elements long)
        """
        return len(self.prevCosts) >= self.params.limitSimulPaths \
            and cost > self.prevCosts[-1]

    def _findShortestPaths(self, pathProcessor, shapeEntry, gtfsPointsPrev, gtfsPoints, baseMap, avoidRestartCode = 0):
        """
        _findShortestPaths coordinates the creation of a list of new tree nodes for each of the reachable new points.
        This is mostly an implementation of "FindShortestPath" in Figure 2 of Perrine et al., 2015.

        @type pathProcessor: graph.WalkPathProcessor
        @type shapeEntry: ShapesEntry
        @type gtfsPointsPrev: list<PathEnd> 
        @type gtfsPoints: list<PathEnd>
        @type baseMap: graph.Map
        @param avoidRestartCode: 0 to allow restarts; 1 to allow restarts but suppress message; 2 to avoid restarts
        @type avoidRestartCode: int
        @rtype list<PathEnd>
        """
        # Initialize the list of costs that will be used to reduce the number of path-finding iterations:
        self.prevCosts = []
        
        # Then, for each previous GTFS tree entry, find the shortest path to each current GTFS tree entry:
        # (On the first time through, this loop will be skipped).
        iterList = gtfsPointsPrev if gtfsPointsPrev else [None]
        for gtfsPointPrev in iterList:
            "@type gtfsPointPrev: PathEnd"
            for gtfsPoint in gtfsPoints:
                "@type gtfsPoint: PathEnd"
                # Calculate path from gtfsPointPrev to vistaPoint.
                traversed, distance, cost, totalLinkCount = ([], 0.0, self.scoreFunction(None, 0.0, gtfsPoint.pointOnLink), 0) if not gtfsPointPrev \
                    else pathProcessor.walkPath(gtfsPointPrev.pointOnLink, gtfsPoint.pointOnLink, gtfsPointPrev.totalCost, gtfsPointPrev.totalLinkCount)
                    
                if traversed is not None:
                    # A valid path was found:
                    if (gtfsPoint.prevTreeNode is None) or ((gtfsPoint.prevTreeNode is not None) \
                                    and (gtfsPointPrev.totalCost + cost < gtfsPoint.totalCost)):
                        # This is the first proposed parent, or the proposed parent is cheaper than what
                        # is there.  Set it:
                        gtfsPoint.prevTreeNode = gtfsPointPrev
                        gtfsPoint.routeInfo = traversed
                        gtfsPoint.totalLinkCount = totalLinkCount
                        if gtfsPointPrev is not None:
                            gtfsPoint.totalCost = gtfsPointPrev.totalCost + cost
                            gtfsPoint.totalDist = gtfsPointPrev.totalDist + distance
                        else:
                            gtfsPoint.totalCost = cost
                            gtfsPoint.totalDist = 0
                        if len(self.prevCosts) < self.limitSimulPaths:
                            self.prevCosts.append(gtfsPoint.totalCost)
                        else:
                            self.prevCosts[-1] = gtfsPoint.totalCost
                        self.prevCosts[:] = sorted(self.prevCosts[:])
                                                
        # Clean up tree entries that didn't get assigned to a parent:
        if len(gtfsPointsPrev) > 0:
            gtfsPointsWork = []
            for gtfsPoint in gtfsPoints:
                if gtfsPoint.prevTreeNode is not None:
                    gtfsPointsWork.append(gtfsPoint)
        else:
            gtfsPointsWork = gtfsPoints
                                        
        # Warn if we ended up with nothing and move on to the next GTFS point:
        if (len(gtfsPointsWork) == 0) and (avoidRestartCode < 2):
            if (avoidRestartCode < 1) and (len(gtfsPointsPrev) > 0):
                # Warn if we are not at the start and we didn't find valid VISTA points.
                if self.logFile is not None:
                    shapeTypeStr = "GTFS shape"
                    print("WARNING: No VISTA paths were found for %s %s, sequence %d." \
                          % (shapeTypeStr, str(shapeEntry.shapeID), shapeEntry.shapeSeq), file=self.logFile)
            
            # Figure out which if the previous paths is the cheapest.
            gtfsPointRestart = None
            "@type gtfsPointRestart: PathEnd"
            if len(gtfsPointsPrev) > 0:
                for gtfsPointPrev in gtfsPointsPrev:
                    "@type gtfsPointPrev: PathEnd"
                    if (gtfsPointRestart is None) or (gtfsPointPrev.totalCost < gtfsPointRestart.totalCost):
                        gtfsPointRestart = gtfsPointPrev
            
            # Mark a "break" in the continuity and link up with the newer candidates.  Keep a limited set.
            gtfsPoints = gtfsPoints[0:self.limitSimulPaths]
            if gtfsPointRestart is not None:
                for gtfsPoint in gtfsPoints:
                    "@type gtfsPoint: PathEnd"
                    gtfsPoint.restart = True
                    gtfsPoint.prevTreeNode = gtfsPointRestart
                    
                    # Fake a distance and cost from the linear distance so that we something to report later.
                    distance = linear.getNorm(gtfsPointRestart.pointOnLink.pointX, gtfsPointRestart.pointOnLink.pointY,
                                       gtfsPoint.pointOnLink.pointX, gtfsPoint.pointOnLink.pointY)
                    gtfsPoint.totalCost = gtfsPointRestart.totalCost + self.scoreFunction(gtfsPointRestart.pointOnLink,
                        distance, gtfsPoint.pointOnLink)
                    gtfsPoint.totalDist = gtfsPointRestart.totalDist + distance
        else:
            # Trim off lowest-scoring paths:
            gtfsPointsWork.sort(key = operator.attrgetter('totalCost'))
            gtfsPoints = gtfsPointsWork[0:self.limitSimulPaths]
            
        return gtfsPoints            

    def constructPath(self,
                      trackpoints: tuple[graph.Map.Trackpoint],
                      baseMap: graph.Map,
                      linkList: list[Hashable] | None = None) -> list[PathEnd] | None:
        """
        constructPath goes through a list of trackpoints and finds the shortest path through the given baseMap.
        This roughly corresponds with algorithms "WalkTrack" and "TrackpointArrives" in Figure 2 of Perrine et al. 2015.
        """
        # TODO: Rename gtfsPointsPrev to something independent of GTFS.
        gtfsPointsPrev: list[PathEnd] = []

        pathProcessor: graph.WalkPathProcessor \
            = graph.WalkPathProcessor(self, baseMap,
                self.params.limitDirectDist, self.params.limitPathDist,
                self.params.limitDirectDistRev, self.maxHops, linkList)
        shapeCtr: int
        startInvalidCheckFlag: bool = True
        startValidIndex: int = 0
        lastValidIndex: int = -1
        invalidCtr: int = 0
        
        logging.info("Building path...")
        
        # TODO: Rename shapeEntry to "trackpoint".
        shapeEntry: graph.Map.Trackpoint
        for shapeCtr, shapeEntry in enumerate(trackpoints):
            
            if shapeCtr % 10 == 0:
                logging.info("   ... %d of %d", shapeCtr, len(trackpoints))

            # TODO: move the forceLinks stuff to baseMap.findPointsOnLinks().
            closestLinks: list[graph.Map.PointOnLink]
            if self.forceLinks and shapeCtr < len(self.forceLinks) \
                    and self.forceLinks[shapeCtr] is not None:
                # Custom behavior for forcing the use of a limited set of links:
                # TODO: Why not make a Map out of the subset, and it will be more versatile?
                closestLinks = []
                linkID: Hashable
                link: graph.Map.LinkRecord
                for linkID in self.forceLinks[shapeCtr]:
                    link = baseMap.getLinkByID(linkID)
                    dist, percentAlong, isPerpendicular, pointAlong = baseMap.pointDist(shapeEntry, link)
                    closestLinks.append(graph.Map.PointOnLink(link, percentAlong, not isPerpendicular, dist, pointAlong))
                closestLinks.sort(key = operator.attrgetter('refDist'))
            else:
                # Normal behavior: search among all links:
                closestLinks = baseMap.findPointsOnLinks(shapeEntry, self.params.searchRadius, self.params.radiusPrimary,
                                self.params.radiusSecondary, (gtfsPoint.pointOnLink for gtfsPoint in gtfsPointsPrev),
                                self.params.limitClosestPoints)
                        
            if not closestLinks:
                lastValidIndex = shapeCtr
                invalidCtr += 1
                logging.warning("No closest links were found for "
                                f"trackpoint {shapeEntry.id}, sequence {shapeEntry.seq}.")
                continue
            else:
                if startInvalidCheckFlag:
                    startInvalidCheckFlag = False
                    startValidIndex = lastValidIndex
                invalidCtr = 0
            
            # Initialize blank endpoint entries:
            endPoints: list[PathEnd] = [PathEnd(shapeEntry, endPoint) for endPoint in closestLinks]
            
            # Find the shortest paths from gtfsPointsPrev to the handful of closestVISTA points:
            # (We're adding another layer to the tree, and previous tree nodes can be found by accessing
            # PathEnd.prevTreeNode)
            gtfsPointsPrev = self._findShortestPaths(pathProcessor, shapeEntry, gtfsPointsPrev, endPoints, baseMap)

        if startInvalidCheckFlag:
            startValidIndex = len(trackpoints)

        # Additional reporting on points found and not found:
        reportStr: str = ""
        missingEnds: int = 0
        if startValidIndex > 0:
            reportStr = f"{startValidIndex} are missing from the start"
            missingEnds = startValidIndex
        if lastValidIndex == len(trackpoints) and startValidIndex < len(trackpoints):
            if startValidIndex > 0:
                reportStr += " and "
            reportStr += f"{invalidCtr} are missing from the end"
            missingEnds += invalidCtr
        if len(reportStr) > 0:
            logging.warning(f"Out of {len(trackpoints)} georeference points, {reportStr}.")
        if float(missingEnds) / len(trackpoints) > self.tossRatio:
            logging.warning(f"Aborting ID {shapeEntry.id}.")
            return None

        # Now, extract the shortest path. First, find the end that has the cheapest cost:
        logging.info("Finalizing path...")
        gtfsPoint = None
        "@type gtfsPoint: PathEnd"
        if len(gtfsPointsPrev) > 0:
            for gtfsPointPrev in gtfsPointsPrev:
                "@type gtfsPointPrev: PathEnd"
                if (gtfsPoint is None) or (gtfsPointPrev.totalCost < gtfsPoint.totalCost):
                    gtfsPoint = gtfsPointPrev
                    
        # Then, follow that end to the beginning:
        ret = []
        "@type ret: list<PathEnd>"
        while gtfsPoint is not None:
            ret.append(gtfsPoint)
            gtfsPoint = gtfsPoint.prevTreeNode
            
        # Reverse the order of the list to go from start to end.
        return ret[::-1]

    @staticmethod
    def _findNextRestart(gtfsPath, startIndex = 0):
        """
        Goes through the list of GTFS points and gives the index of the point before a restart.
        @type gtfsPath: list<PathEnd>
        @type startIndex: int
        @return The index before the next restart, or -1 if not found.
        @rtype int
        """
        while startIndex < len(gtfsPath) and not gtfsPath[startIndex].restart:
            startIndex += 1
        if startIndex >= len(gtfsPath):
            return -1
        return startIndex 

    def setRefineParams(self, termRefactorRadius):
        """
        setRefineParams() sets the parameters that are specific to refining paths.
        @param termRefactorRadius: The radius around restart points that cause tree points to be reevaluated.
        @type termRefactorRadius: float
        """
        self.termRefactorRadius = termRefactorRadius
        self.termRefactorRadiusSq = termRefactorRadius ** 2
        
    def setForceLinks(self, forceLinks: Sequence[Iterable[Hashable]] | None):
        """
        Forces refinePath() to use specific links.
        @param forceLinks: A list of sets of links or None values where each element corresponds with the
            oldGTFSPath list passed into refinePath(). Set to None to disable entirely (default).
        @type forceLinks: list<graph.GraphLink>
        """
        self.forceLinks = forceLinks

    def _tryTreeStack(self, pathProcessor, oldTreeNode, prevTreeNodes, baseMap, evalCode, firstFlag, pathIndex=None):
        """
        _tryTreeStack() is a potentially recursively called internal worker method that reevaluates tree indices and
        generates new tree nodes.
        @type pathProcessor: graph.WalkPathProcessor
        @type oldTreeNode: PathEnd
        @type prevTreeNodes: list<PathEnd>
        @type baseMap: graph.GraphLib
        @param evalCode: Use 0: no evaluation, 1: full evaluation, 2: wrap up loose ends
        @type evalCode: int
        @type firstFlag: bool
        @param pathIndex: This must be provided for path refining that uses forced links.
        @type pathIndex: int
        @return The list of new tree nodes, and then the eval code that was most recently used. 
        @rtype: list<PathEnd>, int        
        """
        prevPointsOnLinks = [prevTreeNode.pointOnLink for prevTreeNode in prevTreeNodes]
        curListAll = []
        "@type curListALl: list<PathEnd>"

        if firstFlag:        
            self.shapeScatterCache = None
        
        # Are we in an area that requires reevaluation (e.g. restart)?  Deal with shape points here:
        if evalCode > 0:
            if evalCode == 1:
                # Check if we had found all of the shape proximity points already:
                if self.shapeScatterCache is None:
                    if pathIndex is not None and self.forceLinks is not None and pathIndex < len(self.forceLinks) \
                            and self.forceLinks[pathIndex] is not None:
                        # Specialized operation: force the use of the given link:
                        self.shapeScatterCache = []
                        link = self.forceLinks[pathIndex]
                        distSq, linkDist, perpendicular = link.pointDistSq(oldTreeNode.shapeEntry.pointX, oldTreeNode.shapeEntry.pointY)
                        self.shapeScatterCache.append(graph.PointOnLink(link, linkDist, not perpendicular, math.sqrt(distSq)))
                    else:
                        # Normal operation: find closest limitClosestPoints points to the shape point among all links: 
                        self.shapeScatterCache = baseMap.findPointsOnLinks(oldTreeNode.shapeEntry.pointX,
                            oldTreeNode.shapeEntry.pointY, self.searchRadius, self.radiusPrimary,
                            self.radiusSecondary, prevPointsOnLinks, self.limitClosestPoints)
                
                # Create new PathEnd objects:
                gtfsPoints = len(self.shapeScatterCache) * []
                "@type gtfsPoints: list<PathEnd>"
                for vistaPoint in self.shapeScatterCache:
                    "@type vistaPoint: graph.PointOnLink"
                    gtfsPoint = PathEnd(oldTreeNode.shapeEntry, vistaPoint)
                    "@type gtfsPoint: PathEnd"
                    gtfsPoints.append(gtfsPoint)
            elif evalCode == 2:
                # We are getting all previous points to converge down on one point, preserving the best one:
                gtfsPoints = [oldTreeNode.cleanCopy()]
                "@type gtfsPoints: list<PathEnd>"
    
            if gtfsPoints:
                # Find the shortest paths from gtfsPointsPrev to the handful of closestVISTA points:
                curList = self._findShortestPaths(pathProcessor, oldTreeNode.shapeEntry, prevTreeNodes,
                    gtfsPoints, baseMap, 1 if firstFlag else 2)
                # Check for restarts and penalize.  Only keep the first (cheapest) restart.  Meanwhile, append to
                # the results list:
                restartFlag = False
                for gtfsPoint in curList:
                    if gtfsPoint.restart:
                        if not restartFlag:
                            gtfsPoint.totalCost *= RESTART_PENALTY_MULT
                            curListAll.append(gtfsPoint)
                            restartFlag = True
                    else:
                        curListAll.append(gtfsPoint)
            else:
                # Evidently we didn't find any candidate points. In this case, duplicate the previously matched link:
                for prevTreeNode in prevTreeNodes:
                    curTreeNode = copy.copy(oldTreeNode)
                    curTreeNode.prevTreeNode = prevTreeNode
                    dist = linear.getNorm(oldTreeNode.shapeEntry.pointX, oldTreeNode.shapeEntry.pointY, curTreeNode.pointOnLink.pointX, curTreeNode.pointOnLink.pointY)
                    curTreeNode.totalDist += dist * RESTART_PENALTY_MULT
                    curTreeNode.totalCost += dist * RESTART_PENALTY_MULT
                    curListAll.append(curTreeNode)
                        
        elif firstFlag:
            # This happens if we are not reevaluating the existing paths at all.
            # TODO: The total cost isn't being added up properly here.  Try combining the evalCode 0 and 2 parts to
            # get the system to retrace the steps that had been traversed before.
            curTreeNode = copy.copy(oldTreeNode)
            "@type curTreeNode: PathEnd"
            if not prevTreeNodes: # This happens on the first element of a path.
                prevTreeNodes.append(None)
            assert len(prevTreeNodes) == 1 # There should just be one of these because we're drawing from a final tree.
            curTreeNode.prevTreeNode = prevTreeNodes[0] 
            curListAll.append(curTreeNode)
                
        if firstFlag:
            if evalCode == 2:
                # Only keep the best result when we converge down to one point:
                curListAll.sort(key = operator.attrgetter('totalCost'))
                curListAll = [curListAll[0]]
        
        # Limit the number of results:
        # TODO: Re-enable if needed?
        #curListAll = curListAll[0:self.limitSimulPaths]
        
        #curListAll.sort(key = operator.attrgetter('totalCost'))
        
        return curListAll, evalCode

    def refinePath(self, oldGTFSPath, baseMap):
        """
        refinePath goes through existing GTFS points and  tries to route from a restart. Uses termRefactorRadius.
        @type oldGTFSPath: list<PathEnd>
        @type baseMap: graph.GraphLib
        @rtype: list<PathEnd>
        """
        if self.logFile is not None:
            print("INFO: Refining path...", file=self.logFile)
            
        treeNodes = []
        
        pathProcessor = graph.WalkPathProcessor(self, self.limitDirectDist, self.limitPathDist, self.limitDirectDistRev,
            self.maxHops)
        "@type pathProcessor: graph.WalkPathProcessor"

        oldTreeNodeIndex = 0
        nextRestartIndex = -1
        evalCode = 0 # 0 = not in restart zone; 1 = in restart zone; 2 = tidying up after restart zone.
        while oldTreeNodeIndex < len(oldGTFSPath):
            # Check to see if we need to find the next restart:
            if (oldTreeNodeIndex == 0) or ((evalCode != 1) and (nextRestartIndex != -1) and (nextRestartIndex < oldTreeNodeIndex)):
                nextRestartIndex = self._findNextRestart(oldGTFSPath, nextRestartIndex + 1)
                
            # Check for restart point. Check if we are in the radius of the last known good point or the restart point.
            if (nextRestartIndex >= 1 and linear.getNormSq(oldGTFSPath[oldTreeNodeIndex].pointOnLink.pointX,
                        oldGTFSPath[oldTreeNodeIndex].pointOnLink.pointY, oldGTFSPath[nextRestartIndex - 1].pointOnLink.pointX,
                        oldGTFSPath[nextRestartIndex - 1].pointOnLink.pointY) < self.termRefactorRadiusSq) \
                    or (nextRestartIndex >= 0 and linear.getNormSq(oldGTFSPath[oldTreeNodeIndex].pointOnLink.pointX,
                        oldGTFSPath[oldTreeNodeIndex].pointOnLink.pointY, oldGTFSPath[nextRestartIndex].pointOnLink.pointX,
                        oldGTFSPath[nextRestartIndex].pointOnLink.pointY) < self.termRefactorRadiusSq):
                if evalCode == 0:
                    evalCode = 1 # Full reevaluation
                    if self.logFile is not None:
                        print("INFO: Enter restart zone at shapeID %s, seq %d..." % (str(oldGTFSPath[oldTreeNodeIndex].shapeEntry.shapeID),
                            oldGTFSPath[oldTreeNodeIndex].shapeEntry.shapeSeq), file=self.logFile)
            else:
                # Tie up loose ends if the previous round had new points found.
                if evalCode == 1:
                    evalCode = 2
                    shapeTypeStr = "GTFS shape"
                    if self.logFile is not None:
                        print("INFO: Exiting zone at %s %s, seq %d..." % (shapeTypeStr, str(oldGTFSPath[oldTreeNodeIndex].shapeEntry.shapeID),
                            oldGTFSPath[oldTreeNodeIndex].shapeEntry.shapeSeq), file=self.logFile)
                        
            if evalCode == 1:
                if self.logFile is not None:
                    print("INFO:   ... shape seq. %d" % oldGTFSPath[oldTreeNodeIndex].shapeEntry.shapeSeq, file=self.logFile)

            # TODO: Also if a shape point is flagged to be reevaluated.
            
            # Visit this shape point further and figure out how to reevaluate it.
            treeNodes, evalCode = self._tryTreeStack(pathProcessor, oldGTFSPath[oldTreeNodeIndex], treeNodes, baseMap,
                evalCode, True, oldTreeNodeIndex)
            
            if evalCode == 2:
                # We have tied up loose ends; now reset.
                # TODO: To always trace current paths for sanity-check, don't set evalCode to 0.
                evalCode = 0
            
            # Check to see if we have a complete path:
            flag = False
            for treeNode in treeNodes:
                "@type treeNode: PathEnd"
                if not treeNode.restart:
                    flag = True
            if not flag and self.logFile is not None:
                print("WARNING: No VISTA path found into GTFS shpaeID %s, seq %d; restarting." \
                    % (str(oldGTFSPath[oldTreeNodeIndex].shapeEntry.shapeID), oldGTFSPath[oldTreeNodeIndex].shapeEntry.shapeSeq),
                    file = self.logFile)
            oldTreeNodeIndex += 1
            
        # Now, extract the shortest path.  First, find the end that has the cheapest cost:
        if self.logFile is not None:
            print("INFO: Finishing path...", file = self.logFile)
        gtfsPoint = None
        "@type gtfsPoint: PathEnd"
        if len(treeNodes) > 0:
            for treeNodeElem in treeNodes:
                "@type treeNodeElem: PathEnd"
                if (gtfsPoint is None) or (treeNodeElem.totalCost < gtfsPoint.totalCost):
                    gtfsPoint = treeNodeElem
                    
        # Then, follow that end to the beginning:
        ret = []
        "@type ret: list<PathEnd>"
        while gtfsPoint is not None:
            ret.append(gtfsPoint)
            gtfsPoint = gtfsPoint.prevTreeNode
            
        # Reverse the order of the list to go from start to end.
        return ret[::-1]
            
def dumpStandardHeader(outFile=sys.stdout):
    """
    Outputs the CSV header that precedes the info from dumpStandardInfo().
    """
    print("shapeID,shapeSeq,shapeType,linkID,linkDist,totalDist,numLinksTrav,linksTrav", file=outFile)

def dumpStandardInfo(treeNodes, outFile=sys.stdout):
    """
    Outputs the body of a CSV format of VISTA path information.
    @type treeNodes: list<PathEnd>
    """
    for gtfsNode in treeNodes:
        "@type gtfsNode: PathEnd"
        # TODO: Since we don't have hint stuff anymore, we can either remove or repurpose shapeType.
        shapeType = 0
        outStr = "%s,%d,%d,%d,%g,%g" % (str(gtfsNode.shapeEntry.shapeID), gtfsNode.shapeEntry.shapeSeq, shapeType,
                                     gtfsNode.pointOnLink.link.id, gtfsNode.pointOnLink.dist, gtfsNode.totalDist)
        if gtfsNode.restart:
            # A length of -1 shall be a special indication saying that we are restarting and the link list
            # does not exist.
            outStr = outStr + ",-1"
        else:
            outStr = outStr + ",%d" % len(gtfsNode.routeInfo)
            for routeTraverse in gtfsNode.routeInfo:
                "@type routeTraverse: graph.GraphLink"
                outStr = outStr + ",%d" % routeTraverse.id
        print(outStr, file=outFile)

def readStandardDump(baseMap, gtfsShapes, inFile, shapeIDMaker = lambda x: int(x)):
    """
    readStandardDump reconstructs the tree entries that PathEngine had created.
    @type baseMap: graph.GraphLib
    @type gtfsShapes: dict<int, list<gtfs.ShapesEntry>>
    @type inFile: file
    @type shapeIDMaker: function
    @return A dictionary of shapeID to a list of PathEnds
    @rtype dict<int, list<PathEnd>>
    """
    ret = {}
    "@type ret: dict<int, list<PathEnd>>"

    # Sanity check:
    fileLine = inFile.readline()
    if not fileLine.startswith("shapeID,shapeSeq,shapeType,linkID,linkDist,totalDist,numLinksTrav,linksTrav"):
        print("ERROR: The path match file doesn't have the expected header.", file=sys.stderr)
        return None
        
    # Storage place for sequence numbers:
    shapeSeqs = {}
    "@type shapeSeqs: dict<int, int>"
    
    # Go through the lines of the file:
    for fileLine in inFile:
        if len(fileLine) > 0:
            lineElems = fileLine.split(',')
            shapeID = shapeIDMaker(lineElems[0])
            shapeSeq = int(lineElems[1])
            linkID = int(lineElems[3])
            linkDist = float(lineElems[4])
            distTotal = float(lineElems[5])
            linksTravCount = int(lineElems[6])
            if linksTravCount >= 0:
                linksTrav = linksTravCount * [None]
            else:
                # If linksTravCount is -1, then that signifies that we are restarting.  Deal with it later.
                linksTrav = []
            "@type linksTrav: list<graph.GraphNode>"
                        
            # Get the variable-length link list that happens at the end:
            contFlag = False
            for index in range(0, len(linksTrav)):
                linksTravID = int(lineElems[index + 7])
                if linksTravID not in baseMap.linkMap:
                    print("WARNING: The path match file refers to a nonexistent link ID %d." % linksTravID, file=sys.stderr)
                    contFlag = True
                    break
                linksTrav[index] = baseMap.linkMap[linksTravID]
            if contFlag:
                # This is run if the break above is run.
                continue

            # Resolve the shapeID list:
            if shapeID not in gtfsShapes:
                print("WARNING: The path match file refers to a nonexistent shape ID %s." % str(shapeID), file=sys.stderr)
                continue
            shapeElems = gtfsShapes[shapeID]
            "@type shapeElems: list<gtfs.ShapesEntry>"
            
            # Set up the shape index cache to reduce linear searching later on:
            if shapeID not in shapeSeqs:
                shapeSeqs[shapeID] = -1
                
            # Resolve the link object:
            if linkID not in baseMap.linkMap:
                print("WARNING: The path match file refers to a nonexistent link ID %d." % linkID, file=sys.stderr)
                continue
            link = baseMap.linkMap[linkID]
            "@type link: graph.GraphLink"
            
            # Resolve the shape entry:
            shapeEntry = None
            "@type shapeEntry: gtfs.ShapesEntry"
            for index in range(shapeSeqs[shapeID] + 1, len(shapeElems)):
                if shapeElems[index].shapeSeq == shapeSeq:
                    shapeEntry = shapeElems[index]
                    shapeSeqs[shapeID] = index
                    break
            if shapeEntry is None:
                print("WARNING: No GTFS shape entry for shape ID: %s, seq: %d; check for out of order."
                      % (str(shapeID), shapeSeq), file = sys.stderr)
                continue
            
            # Recalculate parameters needed for the tree node:
            pointX, pointY = baseMap.gps.gps2feet(shapeEntry.lat, shapeEntry.lng)
            distRef, distLinear, perpFlag = baseMap.linkMap[linkID].pointDist(pointX, pointY)
            pointOnLink = graph.PointOnLink(link, distLinear, not perpFlag, distRef)
            newEntry = PathEnd(shapeEntry, pointOnLink)
            newEntry.totalCost = distTotal # TotalCost won't be available.
            newEntry.totalDist = distTotal
            if linksTravCount == -1:
                # We are restarting the path and don't have complete information up to this point.
                newEntry.restart = True
            
            # Reconstruct the node list:
            newEntry.routeInfo = linksTrav
            
            if shapeID not in ret:
                ret[shapeID] = []
            else:
                # Restore the previous tree entry linkage:
                newEntry.prevTreeNode = ret[shapeID][-1]
            ret[shapeID].append(newEntry)

    # Return the tree nodes:
    return ret

