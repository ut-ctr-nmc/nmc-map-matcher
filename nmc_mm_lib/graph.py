"""
graph.py: Links and nodes for graph models; also a breadth-first search.
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

from collections.abc import Hashable, Iterable, Sequence, Generator
from typing import NamedTuple, MutableMapping, Any, Self
import shapely
from shapely.ops import transform
import networkx
import pyproj
from pyproj.enums import TransformDirection

# from nmc_mm_lib.path_engine import PathEngine
# TODO: Avoid circular reference; bring in score functions through other means


def hasMoreThan(iterable: Iterable[Any], count: int = 0) -> bool:
    """
    Determines if the iterable has more than the specified count of items.

    @param iterable: The iterable to check.
    @param count: The count to compare against (default: 0).
    @return: True if the iterable has more than count items.
    """
    iterator = iter(iterable)
    for _ in range(count + 1):
        try:
            next(iterator)
        except StopIteration:
            return False
    return True


def hasExactly(iterable: Iterable[Any], count: int = 0) -> bool:
    """
    Determines if the iterable has exactly the specified count of items.

    @param iterable: The iterable to check.
    @param count: The count to compare against (default: 0).
    @return: True if the iterable has exactly count items.
    """
    iterator = iter(iterable)
    for _ in range(count):
        try:
            next(iterator)
        except StopIteration:
            return False
    try:
        next(iterator)
        return False
    except StopIteration:
        return True


class Trackpoint(NamedTuple):
    """
    A container for a geocoordinate, usually part of a series
    """

    lonHoriz: float
    latVert: float
    point: shapely.geometry.Point
    id: Hashable | None = None
    seq: int | float | None = None


class Map:
    """
    Map is a container for a node-link graph, maintained internally in
    a Shapely SRTree for spatial indexing.
    """

    fromCRS: pyproj.CRS
    workingCRS: pyproj.CRS
    transformer: pyproj.Transformer
    graph: networkx.DiGraph = networkx.DiGraph()
    edgeIndexLookup: tuple["LinkRecord", ...] = tuple()
    linkIDLookup: dict[Hashable, "LinkRecord"] = {}
    tree: shapely.strtree.STRtree | None = None
    trackpointCtr: int = 0  # For auto-assigning trackpoint sequence numbers
    eqCutoff: int = 3  # Decimal places for equality checks

    def __init__(
        self, fromCRS: str = "EPSG:4326", workingCRS: str = "EPSG:3857"  # GPS
    ):
        """
        Initializes an empty Map.

        @param fromCRS: The coordinate reference system of input geometries
            (Default: GPS).
        @param workingCRS: The coordinate reference system to use for internal
            calculations (Default: Web Mercator).
        """
        self.fromCRS = pyproj.CRS(fromCRS)
        self.workingCRS = pyproj.CRS(workingCRS)
        self.transformer = pyproj.Transformer.from_crs(
            self.fromCRS, self.workingCRS, always_xy=True
        )

    def addNode(
        self,
        nodeID: Hashable,
        lonHoriz: float,
        latVert: float,
        score: float = 0.0,
        metadata: dict = {},
    ):
        """
        Adds a node to the map.

        @param nodeID: The unique identifier for the node.
        @param lonHoriz: The longitude or horizontal measure of the node.
        @param latVert: The latitude or vertical measure of the node.
        @param score: An optional score (penalty)for the node.
        @param metadata: Optional additional metadata to associate with the node.
        """
        x, y = self.transformer.transform(lonHoriz, latVert)
        self.graph.add_node(
            nodeID,
            x=x,
            y=y,
            score=score,
            lonHoriz=lonHoriz,
            latVert=latVert,
            **metadata,
        )

    class LinkRecord(NamedTuple):
        """
        For recording individual links in edge lookups
        """

        origNodeID: Hashable
        destNodeID: Hashable
        data: MutableMapping[str, Any]
        id: Hashable

        def getLength(self) -> float:
            return self.data["geometry"].length

        def getFirstCoords(self) -> tuple[float, float]:
            return self.data["geometry"].coords[0]

        def getPointAlong(
            self, value: float, normalize: bool = False
        ) -> tuple[float, float]:
            pointAlong: shapely.geometry.Point = self.data["geometry"].interpolate(
                value, normalized=normalize
            )
            return pointAlong.x, pointAlong.y

    def addLink(
        self,
        origNodeID: Hashable,
        destNodeID: Hashable,
        controlPoints: (
            Iterable[Sequence[float]] | shapely.geometry.LineString | None
        ) = None,
        linkID: Hashable | None = None,
        metadata: dict = {},
        hasEndpoints: bool = True,
        flatScore: float = 0.0,
        lengthWeight: float = 1.0,
    ) -> Hashable | None:
        """
        Adds a directed curved link to the map using control points.

        @param origNodeID: The originating node ID.
        @param destNodeID: The destination node ID.
        @param controlPoints: A list of (lon, lat) tuples (or LineString) representing control points for the curve, or None for defaults
        @param linkID: The unique identifier for the link, or None (default) to auto-generate
        @param metadata: Optional additional metadata to associate with the link.
        @param hasEndpoints: Whether the control points include the endpoints (default: True).
        @param flatScore: An optional flat score (penalty) for the link, or 0 use the length as score (default: 0).
        @param lengthWeight: A multiplier to apply to the length when computing the score, or 0 to deactivate (default: 1.0).
        @param controlPoints: A list of (lon, lat) tuples representing control points for the curve.
        @return: The linkID that was incorporated into the link, or None if IDs weren't used.
        """
        # First, convert controlPoints if None or LineString:
        if isinstance(controlPoints, shapely.geometry.LineString):
            controlPoints = controlPoints.coords
        if controlPoints is None:
            controlPoints = []
            hasEndpoints = False
        else:
            # Transform control points to working projection:
            lonlats = [(lon, lat) for lon, lat in controlPoints]
            controlPoints = [
                self.transformer.transform(lon, lat) for lon, lat in lonlats
            ]

        # Next, add endpoints if needed:
        if not hasEndpoints:
            origNode = self.graph.nodes[origNodeID]
            destNode = self.graph.nodes[destNodeID]
            controlPoints = (
                [(origNode["x"], origNode["y"])]
                + controlPoints
                + [(destNode["x"], destNode["y"])]
            )

        # Next, convert to LineString:
        controlPoints = shapely.geometry.LineString(controlPoints)

        # Other attributes to store:
        if linkID is None:
            linkID = f"{origNodeID}->{destNodeID}"
        self.addLineStringLink(
            origNodeID,
            destNodeID,
            controlPoints,
            linkID,
            metadata,
            flatScore=flatScore,
            lengthWeight=lengthWeight,
            alreadyXformed=True,
        )
        # TODO: Consider requiring unique link IDs and using them to look up links,
        # instead of indices that update at the mercy of the tree structure. That
        # would allow dynamic updates of the underlying map, and also allow for
        # independence from the underlying libraries.
        return linkID

    def addLineStringLink(
        self,
        origNodeID: Hashable,
        destNodeID: Hashable,
        geometry: shapely.geometry.LineString,
        linkID: Hashable,
        metadata: dict = {},
        flatScore: float = 0.0,
        lengthWeight: float = 1.0,
        alreadyXformed: bool = False,
    ) -> None:
        """
        Adds a directed link to the map.

        @param origNodeID: The originating node ID.
        @param destNodeID: The destination node ID.
        @param geometry: The geometry of the link as a LineString, already projected.
        @param linkID: The unique identifier for the link.
        @param metadata: Optional additional metadata to associate with the link.
        @param flatScore: An optional flat score (penalty) for the link, or 0 use the length as score (default: 0).
        @param lengthWeight: A multiplier to apply to the length when computing the score, or 0 to deactivate (default: 1.0).
        @param alreadyXformed: Whether the geometry needs to be transformed to working CRS.
        """
        if not alreadyXformed:
            geometry = transform(self.transformer.transform, geometry)
        self.graph.add_edge(
            origNodeID,
            destNodeID,
            geometry=geometry,
            flatScore=flatScore,
            lengthWeight=lengthWeight,
            treeIndex=-1,
            id=linkID,
            **metadata,
        )

    def completeMap(self, treeNeeded=True) -> None:
        """
        Completes the map by building the spatial index.

        @param treeNeeded: Whether to build the spatial index tree (default: True).
        """
        self.edgeIndexLookup = tuple(
            Map.LinkRecord(origNodeID=u, destNodeID=v, data=data, id=data["id"])
            for u, v, data in self.graph.edges(data=True)
        )
        self.linkIDLookup = {}

        index: int
        link: Map.LinkRecord
        for index, link in enumerate(self.edgeIndexLookup):
            link.data["treeIndex"] = index
            self.linkIDLookup[link.id] = link
        if treeNeeded:
            self.tree = shapely.strtree.STRtree(
                [link.data["geometry"] for link in self.edgeIndexLookup]
            )
        else:
            self.tree = None

    def isReverseLink(self, linkA: LinkRecord | int, linkB: LinkRecord | int) -> bool:
        """
        Determines if linkB is the reverse of linkA.

        @param linkA: The first link expressed as an edge index or record.
        @param linkB: The second link, expressed as an edge index or record.
        @return: True if linkB is the reverse of linkA.
        """
        if isinstance(linkA, int):
            linkA = self.edgeIndexLookup[linkA]
        if isinstance(linkB, int):
            linkB = self.edgeIndexLookup[linkB]

        if linkA.origNodeID != linkB.destNodeID or linkA.destNodeID != linkB.origNodeID:
            return False

        geomA: shapely.geometry.LineString = linkA.data["geometry"]
        geomB: shapely.geometry.LineString = linkB.data["geometry"]
        if round(geomA.length, self.eqCutoff) != round(geomB.length, self.eqCutoff):
            return False
        tolerance = 1 / self.eqCutoff
        for (xA, yA), (xB, yB) in zip(geomA.coords, geomB.coords[::-1]):
            if abs(xA - xB) >= tolerance or abs(yA - yB) >= tolerance:
                return False
        return True

    def isSimilarLink(self, linkA: LinkRecord | int, linkB: LinkRecord | int) -> bool:
        """
        Determines if linkA is similar to linkB.

        @param linkA: The first link expressed as an edge index or record.
        @param linkB: The second link, expressed as an edge index or record.
        @return: True if linkB is the reverse of linkA.
        """
        if isinstance(linkA, int):
            linkA = self.edgeIndexLookup[linkA]
        if isinstance(linkB, int):
            linkB = self.edgeIndexLookup[linkB]

        if linkA.origNodeID != linkB.origNodeID or linkA.destNodeID != linkB.destNodeID:
            return False

        geomA: shapely.geometry.LineString = linkA.data["geometry"]
        geomB: shapely.geometry.LineString = linkB.data["geometry"]
        if round(geomA.length, self.eqCutoff) != round(geomB.length, self.eqCutoff):
            return False
        tolerance = 1 / self.eqCutoff
        for (xA, yA), (xB, yB) in zip(geomA.coords, geomB.coords):
            if abs(xA - xB) >= tolerance or abs(yA - yB) >= tolerance:
                return False
        return True

    def outgoingLinks(self, nodeID: Hashable) -> Generator[LinkRecord]:
        """
        Returns a generator of outgoing links from a given node.

        @param nodeID: The node ID to look up.
        @return: A tuple of LinkRecords for outgoing links.
        """
        return (
            self.edgeIndexLookup[treeIndex]
            for u, v, treeIndex in self.graph.edges(nodeID, data="treeIndex")
        )

    def getLinkByID(self, linkID: Hashable) -> LinkRecord | None:
        """
        Retrieves a link by its unique identifier.

        @param linkID: The unique identifier of the link.
        @return: The LinkRecord if found, otherwise None.
        """
        return self.linkIDLookup[linkID] if linkID in self.linkIDLookup else None

    def makeTrackpoint(
        self,
        lonHoriz: float,
        latVert: float,
        ident: Hashable | None = None,
        seq: int | None = None,
    ) -> Trackpoint:
        """
        Makes a Trackpoint from horizontal/vertical coordinates.

        @param lonHoriz: The longitude or horizontal measure of the point.
        @param latVert: The latitude or vertical measure of the point.
        @param ident: An optional unique identifier for the point.
        @param seq: An optional sequence number for the point (or auto-assigned).
        @return: The Trackpoint object.
        """
        x, y = self.transformer.transform(lonHoriz, latVert)
        point = shapely.geometry.Point(x, y)
        if seq is None:
            seq = self.trackpointCtr
            self.trackpointCtr += 1
        return Trackpoint(
            id=ident, lonHoriz=lonHoriz, latVert=latVert, point=point, seq=seq
        )

    def resetTrackpointCtr(self, value: int = 0) -> None:
        """
        Resets the internal trackpoint counter to a specified value.

        @param value: The value to set the counter to (default: 0).
        """
        self.trackpointCtr = value

    class PointOnLink(NamedTuple):
        """
        PointOnLink is a specific point on a link. This is documented in
        Figure 1 of Perrine, et al. 2015 as "point_on_link".
        """

        link: "Map.LinkRecord"  # The link that corresponds with this PointOnLink
        percentAlong: float  # Percentage of distance along the link
        nonPerpPenalty: (
            bool  # "not r", True if there is to be a non-perpendicular penalty applied
        )
        refDist: float  # "d_r", the reference distance, or the working radius from the original search point
        point: shapely.geometry.Point  # The point as it sits on the link

        def getDistanceAlong(self) -> float:
            """
            Calculates the distance along the link from the start to this point.

            @return: The distance along the link.
            """
            geometry: shapely.geometry.LineString = self.link.data["geometry"]
            return geometry.length * self.percentAlong

    def pointDist(
        self, trackPoint: Trackpoint, link: LinkRecord
    ) -> tuple[float, float, bool, shapely.geometry.Point]:
        """
        Determines the distance from a point to a link.

        @param trackPoint: The track point to evaluate.
        @param link: The link to evaluate.
        @return: A tuple of (distance from point to link, percent distance
            along link, whether the point is perpendicular to the link, and
            the point that sits on the link).
        """
        # TODO: Create a nicer return type for this. Use PointOnLink!

        geometry: shapely.geometry.LineString = link.data["geometry"]
        percentAlong: float = geometry.project(trackPoint.point, normalized=True)
        pointAlong: shapely.geometry.Point = geometry.interpolate(
            percentAlong, normalized=True
        )
        dist: float = trackPoint.point.distance(pointAlong)
        isPerpendicular: bool = percentAlong > 0.0 and percentAlong < 1.0
        return dist, percentAlong, isPerpendicular, pointAlong

    def revertPointOnLink(self, pointOnLink: PointOnLink) -> tuple[float, float]:
        """
        Reverts a PointOnLink back to longitude and latitude coordinates.

        @param pointOnLink: The PointOnLink to revert.
        @return: A tuple of (x or longitude, y or latitude) coordinates.
        """
        return self.revertPoint(pointOnLink.point.x, pointOnLink.point.y)

    def revertPoint(self, coordHoriz: float, coordVert: float) -> tuple[float, float]:
        """
        Reverts a projected working space point to the input space (e.g. lon/lat)

        @return: A tuple of (x or longitude, y or latitude) coordinates.
        """
        lon, lat = self.transformer.transform(
            coordHoriz, coordVert, direction=TransformDirection.INVERSE
        )
        return lon, lat

    def findPointsOnLinks(
        self,
        trackPoint: Trackpoint,
        radius: float,
        primaryRadius: float,
        secondaryRadius: float,
        prevPoints: Iterable[PointOnLink],
        limitClosestPoints: int | None = None,
    ) -> list[PointOnLink]:
        """
        findPointsOnLinks searches through the graph and finds all PointOnLinks
        that are within the radius. Then, eligible links are proposed
        primaryRadius distance around the considered point, or secondaryRadius
        distance from the previous map points. Returns an empty list if none
        are found. This corresponds with algorithm "FindPointsOnLinks" in
        Figure 1 of Perrine, et al. 2015. This expects that completeMap() has
        already been run.

        @param trackPoint: The point to search from
        @param radius: Maximum search radius distance
        @param primaryRadius: Maximum distance allowed from search point to link
        @param secondaryRadius: Maximum distance allowed from previous map points to link
        @param prevPoints: Previous PointOnLinks to consider for secondaryRadius matching
        @param limitClosestPoints: Maximum number of closest points to return
        @return: A list of PointOnLink objects found, sorted by increasing reference distance
        """
        if self.tree is None:
            raise RuntimeError(
                "Spatial index tree is not built. Call "
                + "completeMap() with treeNeeded=True before using "
                + "findPointsOnLinks()."
            )

        foundLinks = []

        # Find perpendicular and non-perpendicular PointOnLinks that are within radius.
        indices = self.tree.query(
            trackPoint.point, predicate="dwithin", distance=radius
        )
        for index in indices:
            linkRecord: Map.LinkRecord = self.edgeIndexLookup[index]
            refDist, percentAlong, isPerpendicular, pointAlong = self.pointDist(
                trackPoint, linkRecord
            )

            # Assume that if we are nonperpendicular with respect to the
            # end of the current link, and there are outgoing links, we may
            # in fact be perpendicular to one of those outgoing links, which
            # is more worthwhile.
            if not isPerpendicular and percentAlong >= 1.0:
                if self.graph.out_degree(linkRecord.destNodeID) > 0:
                    continue

            # A candidate:
            pointOnLink = Map.PointOnLink(
                linkRecord, percentAlong, not isPerpendicular, refDist, pointAlong
            )

            if refDist <= primaryRadius:
                # Immediate consideration if we are in the primary radius:
                foundLinks.append(pointOnLink)
            else:
                # Check to see if the point is close to a previous point. This allows candidate links to be tracked
                # that can possibly correspond with missing geometry, such as a bus going through a parking lot that
                # isn't represented in the underlying map.
                prevPoint: Map.PointOnLink
                for prevPoint in prevPoints:
                    dist = pointAlong.distance(prevPoint.point)
                    if dist < secondaryRadius:
                        # We have a winner:
                        # @TODO: Determine if we want to add a penalty for the second radius match.
                        foundLinks.append(pointOnLink)
                        break

        foundLinks.sort(key=lambda entry: entry.refDist)
        if limitClosestPoints is not None:
            foundLinks[:] = foundLinks[:limitClosestPoints]
        return foundLinks


class WalkPathProcessor:
    """
    WalkPathProcessor contains methods used to conduct the walkPath algorithm.
    It maintains a cache that persists in-between individual pathfinding
    operations.
    """

    map: Map
    uTurnInterPenalty: (
        float | None
    )  # Add this penalty to U-turns, or None if U-turns not allowed
    uTurnDeadEndPenalty: (
        float | None
    )  # Add this penalty to U-turns at dead-ends, or None for uTurnInterPenalty
    pathEngine: "PathEngine"  # The object that instanciates this class.
    backCache: dict[
        Hashable, dict[Hashable, Map.LinkRecord]
    ]  # Caches previous walkPathoperations to accelerate
    winner: "Next | None"  # Records the winning queue element
    processingQueue: list[
        "PathElement"
    ]  # Processing queue to facilitate the breadth-first search
    pointOnLinkOrig: Map.PointOnLink  # For internal record-keeping
    pointOnLinkDest: Map.PointOnLink  # For internal record-keeping
    linkList: (
        list[Hashable] | None
    )  # Constrains matching to this list of links IDs, if provided
    # TODO: Can the LinkList be more like a tree, or does it need to be?
    limitDistance: float
    limitRadiusRev: float
    limitSteps: int
    limitRadius: float
    backtrackScore: float
    queueCounter: int

    def __init__(
        self,
        pathEngine: "PathEngine",
        map: Map,
        limitRadius: float,
        limitDistance: float,
        limitRadiusRev: float,
        limitSteps: int,
        linkList: list[Hashable] | None = None,
    ):
        """
        This sets the parameters that are final for the entire walkPath
        algorithm execution:
        """
        self.pathEngine = pathEngine
        self.map = map
        self.limitDistance = limitDistance
        self.limitRadiusRev = limitRadiusRev
        self.limitSteps = limitSteps
        self.limitRadius = limitRadius

        self.uTurnInterPenalty = None  # Disable U-turns in intersections
        self.uTurnDeadEndPenalty = 50  # Allow U-turns at dead-ends

        # walkPath cache:
        self.backCache = {}

        # Keep the running score:
        self.backtrackScore = limitDistance

        # Record the winning queue element:
        self.winner = None

        # For tie-breaking when dealing with the priority queue.
        self.queueCounter = 0

        # List of required links for transit purposes.
        self.linkList = linkList

    class Next(NamedTuple):
        """
        Allows path match requests to be queued. Each of these represents a
        traversal from the start of incomingLink to the starts of the next
        possible links. The walkPath() method will create new Next instances
        for each of those possible links and enqueues them in the priority
        queue that coordinates the pathfinding operations.
        """

        # The previous Next structure that led to this one.
        prevStruct: Self | None
        incomingLink: Map.LinkRecord  # The link that we are to traverse.
        linkListIndex: int  # The count of how many links have been traversed
        # The total distance from the origin PointOnLink to the current location.
        distance: float
        # The total calculated cost from the origin PointOnLink to the current location.
        cost: float
        # The number of steps traversed from the origin PointOnLink to incomingLink.
        stepCount: int
        backtrackSet: set[
            Hashable
        ]  # A set of link IDs for all links that had already been traversed.

    def createNext(
        self,
        prevStruct: Next | None,
        incomingLink: Map.LinkRecord,
        startupCost: float = 0.0,
        linkListIndex: int = 0,
    ) -> Next:
        """
        Initializes elements that are stored within a new instance,
        performing necessary calculations

        @param prevStruct: The previous Next structure that led to this one
        @param incomingLink: The link that we are to traverse
        @param startupCost: The starting cost to add to this path segment
        @param linkListIndex: The count of how many links have been traversed
        @return: The newly created Next object
        """
        # linkDistPotential is the distance remaining on the link to be traversed plus prior journey.
        linkDistPotential: float
        stepCount: int
        if prevStruct is None:
            # First-time initialization:
            linkDistPotential = (
                1.0 - self.pointOnLinkOrig.percentAlong
            ) * self.pointOnLinkOrig.link.data["geometry"].length
            # TODO: Use incomingLink length for last term?
            stepCount = 0
        else:
            linkDistPotential = incomingLink.data["geometry"].length
            # linkDistPotential = (
            #     prevStruct.distance + incomingLink.data["geometry"].length
            # )
            stepCount = prevStruct.stepCount + 1

        cost: float
        if incomingLink is self.pointOnLinkDest.link:
            # Last-time initialization; we have hit the destination link:
            # We are stopping midway through this link.  So, subtract off the
            # distance from the end that we aren't traversing.
            linkDistPotential -= (
                1.0 - self.pointOnLinkDest.percentAlong
            ) * incomingLink.data["geometry"].length
            cost = startupCost + self.pathEngine.scoreFunction(
                self.pointOnLinkOrig, linkDistPotential, self.pointOnLinkDest
            )
        else:
            # Normal operation; we hadn't encountered the destination link yet:
            cost = startupCost + self.pathEngine.scoreFunction(
                self.pointOnLinkOrig, linkDistPotential, None
            )

        distance: float = (
            prevStruct.distance if prevStruct else 0.0
        ) + linkDistPotential

        # Make a copy of the set only if it is to change, and add in the new incoming link ID:
        oldBacktrackSet: set[Hashable] = (
            prevStruct.backtrackSet if prevStruct is not None else set()
        )
        if incomingLink.id not in oldBacktrackSet:
            self.backtrackSet = oldBacktrackSet | {incomingLink.id}  # This makes a copy
        else:
            self.backtrackSet = oldBacktrackSet

        return WalkPathProcessor.Next(
            prevStruct,
            incomingLink,
            linkListIndex,
            distance,
            cost,
            stepCount,
            self.backtrackSet,
        )

    class PathElement(NamedTuple):
        """
        PathElement is used to maintain the priority queue for walkPath.
        It contains the cost, a tie-breaker index, and the Next structure
        that is to be processed.
        """

        cost: float
        queueCounter: int
        nextStruct: "WalkPathProcessor.Next"

    # TODO: Create a return type for walkPath.

    def walkPath(
        self,
        pointOnLinkOrig: Map.PointOnLink,
        pointOnLinkDest: Map.PointOnLink,
        startupCost: float = 0.0,  # TODO: !!! Use link.data['flatScore'] !!!
        totalLinkCount: int = 0,
    ) -> tuple[list[Map.LinkRecord] | None, float, float, int]:
        """
        walkPath uses a breadth-first search to find the shortest distance from a given PointOnLink to another PointOnLink and
        returns a list of links representing nodes and following links encountered.  Specify a limiting radius for
        evaluating target nodes, and maximum distance traversed.  Also specify a smaller radius for small distances backwards.
        If nothing is found, then None is returned.  An empty list signifies that the destination is on the same link as the
        origin.
        """
        # Initializations:
        self.pointOnLinkOrig = (
            pointOnLinkOrig  # TOOD: Rearrange methods to keep these local
        )
        self.pointOnLinkDest = pointOnLinkDest
        self.winner = None
        self.backtrackScore = self.limitDistance

        # Are the points too far away to begin with?
        origDestDist = pointOnLinkOrig.point.distance(pointOnLinkDest.point)
        if origDestDist > self.limitRadius:
            return None, 0.0, 0.0, 0

        # Set a reasonable bound for the expected distance in this path search:
        self.backtrackScore = self.limitDistance

        # Set up a queue for the search. Preload the queue with the first starting location:
        self.processingQueue = [
            WalkPathProcessor.PathElement(
                cost=0.0,
                queueCounter=0,
                nextStruct=self.createNext(
                    None, pointOnLinkOrig.link, startupCost, totalLinkCount
                ),
            )
        ]
        self.queueCounter = 0

        # Do the breadth-first search:
        while self.processingQueue:
            self._walkPath(self.processingQueue.pop().nextStruct)

        # Set up the return:
        if self.winner is not None:
            # Iterate through all of the links we have traversed. (Ignore
            # first item because we hadn't technically traversed it).
            retList: list[Map.LinkRecord] = []
            element: WalkPathProcessor.Next = self.winner

            while element.prevStruct is not None:
                retList.append(element.incomingLink)
                element = element.prevStruct
            retList.reverse()
            return (
                retList,
                self.winner.distance,
                self.winner.cost - startupCost,
                self.winner.linkListIndex,
            )
        else:
            # We didn't find anything.
            return None, 0.0, 0.0, 0

    # _walkPath is called internally by walkPath().
    def _walkPath(self, walkPathElem: Next) -> None:
        """
        _walkPath is the internal processing element for the pathfinder.
        @type walkPathElem: _WalkPathNext
        """
        # Check maximum number of steps:
        if walkPathElem.stepCount >= self.limitSteps:
            return

        # Check total distance; we are not interested if we exceed our previous best score:
        if walkPathElem.distance >= self.backtrackScore:
            return

        # Do we exceed the worst cost in the list of simultaneous costs?
        if self.pathEngine.exceedsPreviousCosts(walkPathElem.cost):
            return

        # Are we at the destination?
        if walkPathElem.incomingLink is self.pointOnLinkDest.link:
            # We have a winner!
            self.winner = walkPathElem
            self.backtrackScore = walkPathElem.distance

            # Log the winner into the cache by looking at all of the parent elements:
            if self.pointOnLinkDest.link.id not in self.backCache:
                self.backCache[self.pointOnLinkDest.link.id] = {}
            mappings = self.backCache[self.pointOnLinkDest.link.id]
            "@type mappings: dict<int, GraphLink>"
            if walkPathElem.prevStruct is not None:
                element = walkPathElem.prevStruct
                "@type element: _WalkPathNext"
                while element.prevStruct is not None:
                    if (element.prevStruct.incomingLink.id in mappings) and (
                        mappings[element.prevStruct.incomingLink.id]
                        is element.incomingLink
                    ):
                        break
                    mappings[element.prevStruct.incomingLink.id] = element.incomingLink
                    element = element.prevStruct

            # Process the next queue element:
            return

        # Look at each link that comes out from the current node. First, see
        # if there is a shortcut to our destination already in the cache:
        myList: Iterable[Map.LinkRecord]
        if (self.pointOnLinkDest.link.id in self.backCache) and (
            walkPathElem.incomingLink.id in self.backCache[self.pointOnLinkDest.link.id]
        ):
            myList = (
                self.backCache[self.pointOnLinkDest.link.id][
                    walkPathElem.incomingLink.id
                ],
            )
        else:
            myList = self.map.outgoingLinks(walkPathElem.incomingLink.destNodeID)
        link: Map.LinkRecord
        for link in myList:
            # Filter out U-turns:
            penalty = 0.0
            if (
                self.uTurnDeadEndPenalty != 0 or self.uTurnInterPenalty != 0
            ) and self.map.isReverseLink(walkPathElem.incomingLink, link):
                # Is it a dead-end?
                if hasExactly(
                    self.map.outgoingLinks(walkPathElem.incomingLink.destNodeID), 1
                ):
                    if self.uTurnDeadEndPenalty is None:
                        if self.uTurnInterPenalty is None:
                            continue
                        else:
                            penalty = self.uTurnInterPenalty
                    else:
                        penalty = self.uTurnDeadEndPenalty
                else:
                    if self.uTurnInterPenalty is None:
                        continue
                    else:
                        penalty = self.uTurnInterPenalty
                penalty = self.pathEngine.scoreFunction(
                    None, penalty, None
                )  # TODO: !!! Use stuff in link.data !!!

            # Is this the next link we need to process according to the link list (transit)?
            if (
                self.linkList
                and walkPathElem.linkListIndex + 1 < len(self.linkList)
                and self.linkList[walkPathElem.linkListIndex + 1] != link.id
            ):
                continue
                # TODO: We want to eventually allow the path to be departed and then regained. How to do this? We can create a "path lost" state,
                # and as long as that state is True, then search forward in self.linkList to see if we regain the path. Or, add the indices into
                # the link list into the actual link graph objects (as sets). Departure from a set will incur a penalty, and encounter with a set
                # will allow the index to be reset to the last known value.

            # Had we visited this before?
            if link.data["id"] in walkPathElem.backtrackSet:
                continue
                # TODO: This won't work with park-and-rides where a path loops around on itself. This can possibly be fixed by adding a penalty
                # and allowing the path to be traversed. Turn this on with an option. Execution will probably be a bit slower.

            # Add to the queue for processing later:
            self.queueCounter += 1
            self.processingQueue.append(
                WalkPathProcessor.PathElement(
                    cost=walkPathElem.cost + penalty,
                    queueCounter=self.queueCounter,
                    nextStruct=self.createNext(
                        walkPathElem,
                        link,
                        walkPathElem.cost + penalty,
                        walkPathElem.linkListIndex + 1,
                    ),
                )
            )
