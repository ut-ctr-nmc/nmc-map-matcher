"""
path_engine.py contains logic for matching up trackpoints to base map paths
@author: Kenneth Perrine
@contact: kperrine@utexas.edu
@organization: Network Modeling Center, Center for Transportation Research,
    Cockrell School of Engineering, The University of Texas at Austin
@version: 2.0

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
import operator, copy
import logging

# Multiplier for trackpoint-to-trackpoint evaluations that happen while refining
# on a restart:
RESTART_PENALTY_MULT: Final[float] = 2.0

# Trackpoint processing logging interval:
POINT_LOG_INTERVAL: Final[int] = 10


class PathEnd:
    """
    PathEnd is a single node used within the overall tree structure. This
    roughly equates to the "path_end" data structure outlined in Figure 2 of
    Perrine et al., 2015.
    """

    refPoint: graph.Trackpoint
    pointOnLink: graph.Map.PointOnLink
    totalCost: float  # "s", the total score of the path represented
    prevTreeNode: "PathEnd | None"  # "p", the previous PathEnd step in this path
    totalDist: float  # Relates to total score of the path represented
    totalLinkCount: int  # the number of links that had been traversed
    routeInfo: list[
        graph.Map.LinkRecord
    ]  # "l", a list of map links that have been traversed on the shortest path
    restart: bool  # "r", a Boolean signifying a discontinuity

    def __init__(self, refPoint: graph.Trackpoint, pointOnLink: graph.Map.PointOnLink):
        """
        Sets up values in this object, many of which need to be mutable

        @param refPoint: The trackpoint that this PathEnd corresponds with
        @param pointOnLink: Current point being considered
        """
        self.refPoint = refPoint
        self.pointOnLink = pointOnLink

        self.totalCost = 0.0
        self.prevTreeNode = None
        self.totalDist = 0.0
        self.totalLinkCount = 0
        self.routeInfo = []
        self.restart = False

    def cleanCopy(self) -> "PathEnd":
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
        Used for configuring the desired behavior of the path engine. Remarks
        for each parameter coincide with constants in Perrine et al., 2015.

        @param searchRadius: "k": Radius (m) to search from trackpoint to perpendicular basemap links (default: 100.0)
        @param radiusPrimary: "k_p": Radius (m) to search from trackpoint to new basemap links (default: 100.0)
        @param radiusSecondary: "k_s": Radius (m) to search from basemap perpendicular point to previous point (default: 50.0)
        @param limitPathDist: Path distance (m) to allow new proposed paths from one point to another (default: 500.0)
        @param limitDirectDist: Radius (m) to allow new proposed paths from one point to another (default: 500.0)
        @param limitDirectDistRev: Radius (m) to allow backtracking on a link (e.g. entering an off-map parking lot) (default: 160.0)
        @param distFactor: "f_d": Cost multiplier for linear path distance (default: 1.0)
        @param driftFactor: "f_r": Cost multiplier for distance from trackpoint to its basemap link (default: 2.0)
        @param nonPerpPenalty: "f_p": Penalty multiplier for trackpoints that aren't perpendicular to basemap links (default: 1.5)
        @param limitClosestPoints: "q_p": Number of close-proximity points that are considered for each trackpoint (default: 12)
        @param limitSimulPaths: "q_e": Number of proposed paths (hypotheses) to maintain during pathfinding stage (default: 8)
        @param maxHops: Maximum number of basemap links to pursue in a path-finding operation (default: 12)
        @param tossRatio: Disable the invalidation of short paths (default: 1.0)
        """

        searchRadius: float = 100.0
        radiusPrimary: float = 100.0
        radiusSecondary: float = 50.0
        limitPathDist: float = 500.0
        limitDirectDist: float = 500.0
        limitDirectDistRev: float = 160.0
        distFactor: float = 1.0
        driftFactor: float = 2.0
        nonPerpPenalty: float = 1.5
        limitClosestPoints: int = 12
        limitSimulPaths: int = 8
        maxHops: int = 12
        tossRatio: float = 1.0

    params: Params
    pathPointsPrev: Sequence[PathEnd | None]
    prevCosts: list[float] = []  # A list of limitSimulPaths cost values that can be
    # used to determine if proposed paths are worth traversing.
    shapeScatterCache: list[graph.Map.PointOnLink] | None = None
    forceLinks: Sequence[Iterable[Hashable]] | None = None
    termRefactorRadius: float

    def __init__(self, params: Params = Params()):
        """
        Initializes the path engine with the given parameters; those not
        specified are left to be default

        @param params: Configuration parameters
        """
        self.params = params

    def _gatherWPPParams(self, map: graph.Map) -> graph.WalkPathProcessor.Params:
        """
        _gatherWPPParams gathers the parameters that are needed to initialize a WalkPathProcessor.

        @return The parameters to initialize a WalkPathProcessor.
        """
        # Make local copies of variables for baking into functions:
        # TODO: Is this necessary?
        driftFactor: float = self.params.driftFactor
        nonPerpPenalty: float = self.params.nonPerpPenalty
        distFactor: float = self.params.distFactor
        prevCosts: list[float] = self.prevCosts
        limitSimulPaths: int = self.params.limitSimulPaths

        # Bake functions from parameters given in this PathEngine.
        # TODO: Create a class derived from an "interface" instead of one-off functions.
        def scoreFunction(
            prevGeoPoint: graph.Map.PointOnLink | None,
            distance: float,
            geoPoint: graph.Map.PointOnLink | None,
        ) -> float:
            """
            scoreFunction calculates a cost value given prior path distance, and deviation from the basemap link.
            This corresponds with algorithm "ScoreFunction" in Perrine et al., 2015.

            @param prevGeoPoint: graph.PointOnLink
            @param distance: float
            @param geoPoint: graph.PointOnLink
            @return: Calculated score value
            """
            cost: float
            if prevGeoPoint is None:
                # We are starting anew. Count the "black line distance" from the basemap link to the trackpoint:
                if geoPoint is not None:
                    cost = geoPoint.refDist * driftFactor
                    if geoPoint.nonPerpPenalty:
                        cost = cost * nonPerpPenalty
                else:
                    cost = 0.0
                return cost
            else:
                # We're jumping from one link to another, so add the "black line" distance to the total basemap link distance:
                if geoPoint is not None:
                    cost = geoPoint.refDist * driftFactor
                    if geoPoint.nonPerpPenalty:
                        cost = cost * nonPerpPenalty
                else:
                    cost = 0.0
                return cost + abs(distance) * distFactor
                # Change from Perrine et al., 2015: Use absolute value of distance here because all movement
                # should be incrementing even in cases where a proposed path is moving back and forth on a link
                # because of shape point noise or tiny U-turns.

        # Bake an exceeds checker:
        def exceedsPreviousCosts(cost: float) -> bool:
            """
            Returns true if the given cost value exceeds the most expensive cost already recorded (if the list is
            limitSimulPaths elements long)

            @param cost: The cost value to check
            """
            return len(prevCosts) >= limitSimulPaths and cost > prevCosts[-1]

        return graph.WalkPathProcessor.Params(
            map=map,
            scoreFunction=scoreFunction,
            exceedsPreviousCosts=exceedsPreviousCosts,
            limitPathDist=self.params.limitPathDist,
            limitDirectDist=self.params.limitDirectDist,
            limitDirectDistRev=self.params.limitDirectDistRev,
            limitSteps=self.params.maxHops,
        )

    def _findShortestPaths(
        self,
        wppParams: graph.WalkPathProcessor.Params,
        shapeEntry: graph.Trackpoint,
        pathPoints: list[PathEnd],
        avoidRestartCode: int = 0,
    ) -> list[PathEnd]:
        """
        _findShortestPaths coordinates the creation of a list of new tree nodes for each of the reachable new points.
        This is mostly an implementation of "FindShortestPath" in Figure 2 of Perrine et al., 2015.

        @param avoidRestartCode: 0 to allow restarts; 1 to allow restarts but suppress message; 2 to avoid restarts
        """
        # Initialize the list of costs that will be used to reduce the number of path-finding iterations:
        self.prevCosts.clear()

        # Then, for each previous tree entry, find the shortest path to each current tree entry:
        # (On the first time through, this loop will be skipped).
        iterList: Sequence[PathEnd | None] = (
            self.pathPointsPrev if self.pathPointsPrev else [None]
        )
        pathPointPrev: PathEnd | None
        # TODO: Make pathProcessors here while constructing.
        for pathPointPrev in iterList:
            pathPoint: PathEnd
            for (
                pathPoint
            ) in pathPoints:  # TODO: What if these were found simultaneously?
                pathProcessor: graph.WalkPathProcessor = graph.WalkPathProcessor(
                    wppParams, pathPoint.pointOnLink
                )

                # Calculate path from pathPointPrev to candidate points.
                walkResult = (
                    graph.WalkPathProcessor.PathResult(
                        [],
                        0.0,
                        wppParams.scoreFunction(None, 0.0, pathPoint.pointOnLink),
                        0,
                    )
                    if not pathPointPrev
                    else pathProcessor.walkPath(
                        pathPointPrev.pointOnLink,
                        pathPointPrev.totalCost,
                        pathPointPrev.totalLinkCount,
                    )
                )

                if walkResult.linkList is not None and pathPointPrev:
                    # A valid path was found:
                    if (pathPoint.prevTreeNode is None) or (
                        (pathPoint.prevTreeNode is not None)
                        and (
                            pathPointPrev.totalCost + walkResult.cost
                            < pathPoint.totalCost
                        )
                    ):
                        # This is the first proposed parent, or the proposed parent is
                        # less costly than what was found previously. Set it:
                        # TODO: If more costly parents are to be replaced, we can pass pathPoint's
                        #       recent cost to walkPath(), or pass pathPoint to exceedsPreviousCosts(),
                        #       and check there.
                        pathPoint.prevTreeNode = pathPointPrev
                        pathPoint.routeInfo = walkResult.linkList
                        pathPoint.totalLinkCount = walkResult.linkListIndex
                        if pathPointPrev is not None:
                            pathPoint.totalCost = (
                                pathPointPrev.totalCost + walkResult.cost
                            )
                            pathPoint.totalDist = (
                                pathPointPrev.totalDist + walkResult.distance
                            )
                        else:
                            pathPoint.totalCost = walkResult.cost
                            pathPoint.totalDist = 0
                        if len(self.prevCosts) < self.params.limitSimulPaths:
                            self.prevCosts.append(pathPoint.totalCost)
                        else:
                            self.prevCosts[-1] = pathPoint.totalCost
                        self.prevCosts.sort()

        # Clean up tree entries that didn't get assigned to a parent:
        if len(self.pathPointsPrev) > 0:
            pathPointsWork: list[PathEnd] = []
            for pathPoint in pathPoints:
                if pathPoint.prevTreeNode is not None:
                    pathPointsWork.append(pathPoint)
        else:
            pathPointsWork = pathPoints

        # Warn if we ended up with nothing and move on to the next point:
        if (len(pathPointsWork) == 0) and (avoidRestartCode < 2):
            if (avoidRestartCode < 1) and (len(self.pathPointsPrev) > 0):
                # Warn if we are not at the start and we didn't find valid map points.
                logging.warning(
                    f"No map paths were found for path {shapeEntry.id}, sequence {shapeEntry.seq}."
                )

            # Figure out which of the previous paths is the cheapest.
            pathPointRestart: PathEnd | None = None
            if len(self.pathPointsPrev) > 0:
                pathPointPrev: PathEnd | None
                for pathPointPrev in self.pathPointsPrev:
                    if (pathPointRestart is None) or (
                        pathPointPrev
                        and pathPointPrev.totalCost < pathPointRestart.totalCost
                    ):
                        pathPointRestart = pathPointPrev

            # Mark a "break" in the continuity and link up with the newer candidates.  Keep a limited set.
            pathPoints = pathPoints[0 : self.params.limitSimulPaths]
            if pathPointRestart is not None:
                pathPoint: PathEnd
                for pathPoint in pathPoints:
                    pathPoint.restart = True
                    pathPoint.prevTreeNode = pathPointRestart

                    # Fake a distance and cost from the linear distance so that we something to report later.
                    # TODO: To separate from Shapely, consider adding a distance method to PointOnLink or Map.
                    distance = pathPointRestart.pointOnLink.point.distance(
                        pathPointRestart.pointOnLink.point
                    )
                    pathPoint.totalCost = (
                        pathPointRestart.totalCost
                        + wppParams.scoreFunction(
                            pathPointRestart.pointOnLink,
                            distance,
                            pathPoint.pointOnLink,
                        )
                    )
                    pathPoint.totalDist = pathPointRestart.totalDist + distance
        else:
            # Trim off lowest-scoring paths:
            pathPointsWork.sort(key=operator.attrgetter("totalCost"))
            pathPoints = pathPointsWork[0 : self.params.limitSimulPaths]

        return pathPoints

    def constructPath(
        self,
        trackpoints: Iterable[graph.Trackpoint],
        baseMap: graph.Map,
    ) -> list[PathEnd] | None:
        """
        constructPath goes through a list of trackpoints and finds the shortest path through the given baseMap.
        This roughly corresponds with algorithms "WalkTrack" and "TrackpointArrives" in Figure 2 of Perrine et al. 2015.
        """
        # Various local variable initializations:
        self.pathPointsPrev = []
        startInvalidCheckFlag: bool = True
        startValidIndex: int = 0
        lastValidIndex: int = -1
        invalidCtr: int = 0

        # Create the parameter set for the WalkPathProcessor:
        wppParams: graph.WalkPathProcessor.Params = self._gatherWPPParams(baseMap)

        logging.info("Building path...")

        trackpoint: graph.Trackpoint
        trackCtr: int = -1
        for trackCtr, trackpoint in enumerate(trackpoints):

            if (trackCtr + 1) % POINT_LOG_INTERVAL == 0:
                if isinstance(trackpoints, Sequence):
                    logging.info(f"   ... {trackCtr + 1} of {len(trackpoints)}")
                else:
                    logging.info(f"   ... {trackCtr + 1}")

            # TODO: move the forceLinks stuff to baseMap.findPointsOnLinks().
            closestLinks: list[graph.Map.PointOnLink]
            if (
                self.forceLinks
                and trackCtr < len(self.forceLinks)
                and self.forceLinks[trackCtr] is not None
            ):
                # Custom behavior for forcing the use of a limited set of links:
                # TODO: Why not make a Map out of the subset, and it will be more versatile?
                closestLinks = []
                linkID: Hashable
                link: graph.Map.LinkRecord | None
                for linkID in self.forceLinks[trackCtr]:
                    link = baseMap.getLinkByID(linkID)
                    if link is None:
                        logging.warning(
                            f"forceLinks for index {trackCtr} contains link ID {linkID} that doesn't exist in the map."
                        )
                        continue
                    dist, percentAlong, isPerpendicular, pointAlong = baseMap.pointDist(
                        trackpoint, link
                    )
                    closestLinks.append(
                        graph.Map.PointOnLink(
                            link, percentAlong, not isPerpendicular, dist, pointAlong
                        )
                    )
                closestLinks.sort(key=operator.attrgetter("refDist"))
            else:
                # Normal behavior: search among all links:
                closestLinks = baseMap.findPointsOnLinks(
                    trackpoint,
                    self.params.searchRadius,
                    self.params.radiusPrimary,
                    self.params.radiusSecondary,
                    (
                        pathPoint.pointOnLink
                        for pathPoint in self.pathPointsPrev
                        if pathPoint is not None
                    ),
                    self.params.limitClosestPoints,
                )

            if not closestLinks:
                lastValidIndex = trackCtr
                invalidCtr += 1
                logging.warning(
                    f"No closest links found for trackpoint {trackpoint.id}, seq. {trackpoint.seq}."
                )
                continue
            else:
                if startInvalidCheckFlag:
                    startInvalidCheckFlag = False
                    startValidIndex = lastValidIndex
                invalidCtr = 0

            # Initialize blank endpoint entries:
            endPoints: list[PathEnd] = [
                PathEnd(trackpoint, endPoint) for endPoint in closestLinks
            ]

            # Find the shortest paths from pathPointsPrev to the handful of closest base map points:
            # (We're adding another layer to the tree, and previous tree nodes can be found by accessing
            # PathEnd.prevTreeNode)
            self.pathPointsPrev = self._findShortestPaths(
                wppParams, trackpoint, endPoints
            )

        if not isinstance(trackpoints, Sequence):
            logging.info(f"Finished with {trackCtr + 1} trackpoints.")

        if startInvalidCheckFlag:
            startValidIndex = trackCtr + 1

        # Additional reporting on points found and not found:
        reportStr = ""
        missingEnds = 0
        if startValidIndex > 0:
            reportStr = f"{startValidIndex} are missing from the start"
            missingEnds = startValidIndex
        if lastValidIndex == trackCtr + 1 and startValidIndex < trackCtr + 1:
            if startValidIndex > 0:
                reportStr += " and "
            reportStr += f"{invalidCtr} are missing from the end"
            missingEnds += invalidCtr
        if len(reportStr) > 0:
            logging.warning(f"Out of {trackCtr + 1} georeference points, {reportStr}.")
        if float(missingEnds) / (trackCtr + 1) > self.params.tossRatio:
            logging.warning(f"Aborting ID {trackpoint.id}.")
            return None

        # Now, extract the shortest path. First, find the end that has the
        # cheapest cost. Note that there could be cases where multiple ends
        # have a samilar cost (especially if processing incoming points from a
        # live stream), and it could be appropriate to report multiple
        # candidate paths.
        logging.info("Finalizing path...")
        pathPoint: PathEnd | None = None
        if len(self.pathPointsPrev) > 0:
            pathPointPrev: PathEnd | None
            for pathPointPrev in self.pathPointsPrev:
                if (pathPoint is None) or (
                    pathPointPrev is not None
                    and pathPointPrev.totalCost < pathPoint.totalCost
                ):
                    pathPoint = pathPointPrev

        # Then, follow that end to the beginning:
        ret: list[PathEnd] = []
        while pathPoint is not None:
            ret.append(pathPoint)
            pathPoint = pathPoint.prevTreeNode

        # Reverse the order of the list to go from start to end.
        return ret[::-1]
