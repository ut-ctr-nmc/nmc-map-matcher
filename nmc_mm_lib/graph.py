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
import queue
from typing import Callable, MutableMapping, Any, NamedTuple, Self
from dataclasses import dataclass
import shapely
from shapely.ops import transform
import networkx
import pyproj
from pyproj.enums import TransformDirection


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


@dataclass(frozen=True)
class Trackpoint:
    """
    A container for a geocoordinate, usually part of a series
    """

    lonHoriz: float
    latVert: float
    point: shapely.geometry.Point
    id: Hashable | None
    seq: int | float | None

    def __rshift__(self, other: Self) -> float:
        """
        Calculates the distance between two Trackpoints.

        @param other: The other Trackpoint to measure distance to.
        @return: The distance between the two Trackpoints.
        """
        return self.point.distance(other.point)


class Map:
    """
    Map is a container for a node-link graph, maintained internally in
    a Shapely SRTree for spatial indexing.
    """

    fromCRS: pyproj.CRS
    workingCRS: pyproj.CRS
    transformer: pyproj.Transformer
    graph: networkx.MultiDiGraph = networkx.MultiDiGraph()
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

    @dataclass(frozen=True)
    class LinkRecord:
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

        def getOrigin(self) -> shapely.geometry.Point:
            return shapely.get_point(self.data["geometry"], 0)

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
            key=linkID,
            geometry=geometry,
            flatScore=flatScore,
            lengthWeight=lengthWeight,
            treeIndex=-1,
            **metadata,
        )

    def completeMap(self, treeNeeded=True) -> None:
        """
        Completes the map by building the spatial index.

        @param treeNeeded: Whether to build the spatial index tree (default: True).
        """
        self.edgeIndexLookup = tuple(
            Map.LinkRecord(origNodeID=u, destNodeID=v, data=data, id=key)
            for u, v, key, data in self.graph.edges(data=True, keys=True)
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
            for u, v, keys, treeIndex in self.graph.edges(
                nodeID, data="treeIndex", keys=True
            )
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

    @dataclass(frozen=True)
    class PointOnLink:
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
        origPoint: Trackpoint | None = (
            None  # The original trackpoint that led to this PointOnLink, if any
        )

        def getDistanceAlong(self) -> float:
            """
            Calculates the distance along the link from the start to this point.

            @return: The distance along the link.
            """
            geometry: shapely.geometry.LineString = self.link.data["geometry"]
            return geometry.length * self.percentAlong

        def findDistanceFrom(self, other: "Map.PointOnLink") -> float:
            """
            Calculates the distance from this PointOnLink to another PointOnLink.

            @param other: The other PointOnLink to measure distance to.
            @return: The distance between the two points.
            """
            return self.point.distance(other.point)

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
                linkRecord,
                percentAlong,
                not isPerpendicular,
                refDist,
                pointAlong,
                trackPoint,
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


class GlobalPathCache:
    """
    Logs the shortest path previously found ending at an end link
    """

    pathCache: dict[Hashable, dict[Hashable, Map.LinkRecord]] = {}
    # That's destination -> start -> next link in path

    def clear(self) -> None:
        self.pathCache.clear()

    def trimExcept(self, saveEndLinks: Iterable[Hashable]) -> None:
        """
        Trims the path cache to only keep paths for given end links.

        @param endLinks: An iterable of ending link identifiers to keep in the path cache
        """
        allowed = set(saveEndLinks)
        for k in self.pathCache.keys() - allowed:
            del self.pathCache[k]

    def update(
        self, endLink: Hashable, linkRecrds: Iterable[tuple[Hashable, Map.LinkRecord]]
    ) -> None:
        """
        Updates the path cache with a new path for a given start and end link.

        @param endLink: The ending link identifier
        @param linkRecrds: An iterable of tuples containing start and next link identifiers
        """
        if endLink not in self.pathCache:
            self.pathCache[endLink] = {}
        ref: dict[Hashable, Map.LinkRecord] = self.pathCache[endLink]
        for startLink, nextLink in linkRecrds:
            ref[startLink] = nextLink

    def check(self, startLink: Hashable, endLink: Hashable) -> Map.LinkRecord | None:
        """
        Checks the path cache for a given start and end link.

        @param startLink: The starting link identifier
        @param endLink: The ending link identifier
        @return: The next link in the path if it exists, otherwise None
        """
        return self.pathCache.get(endLink, {}).get(startLink)


class WalkPathProcessor:
    """
    WalkPathProcessor contains methods used to conduct the walkPath algorithm. This is the
    step to find the shortest path between parents and a destination PointOnLink. It
    maintains a cache that persists in-between individual pathfinding operations.
    """

    class Params(NamedTuple):
        """
        Used for configuring the desired behavior of WalkPathProcessor. Remarks
        for each parameter coincide with constants in Perrine et al., 2015.

        @param scoreFunction: Used for scoring, based upon costs derived from pair of PointOnLink objects
        @param exceedsPreviousCosts: A function that checks if a proposed path's cost exceeds the costs of previously found paths
        @param limitPathDist: Path distance (m) to allow new proposed paths from one point to another (default: 500.0)
        @param limitDirectDist: Radius (m) to allow new proposed paths from one point to another (default: 500.0)
        @param limitDirectDistRev: Radius (m) to allow backtracking on a link (e.g. entering an off-map parking lot) (default: 160.0)
        @param limitSteps: Maximum number of basemap links to pursue in a path-finding operation (default: 12)
        @param uTurnInterPenalty: Penalty to add to U-turns in intersections, or None to disable U-turns in intersections (default: None)
        @param uTurnDeadEndPenalty: Penalty to add to U-turns at dead-ends, or None for uTurnInterPenalty (default: 50)
        """

        map: Map
        scoreFunction: Callable[
            [Map.PointOnLink | None, float, Map.PointOnLink | None], float
        ]
        exceedsPreviousCosts: Callable[[float, float | None], bool]
        globalPathCache: GlobalPathCache | None = None
        allTargetLinkIDs: set[Hashable] = set()
        limitPathDist: float = 500.0
        limitDirectDist: float = 500.0
        limitDirectDistRev: float = 160.0
        limitSteps: int = 12
        uTurnInterPenalty: float | None = None
        uTurnDeadEndPenalty: float | None = 50.0

    pointOnLinkDest: Map.PointOnLink
    winner: "Next | None"  # Records the winning queue element
    processingQueue: queue.PriorityQueue[
        "PathElement"
    ]  # Processing queue to facilitate the breadth-first search
    pointOnLinkOrig: Map.PointOnLink  # For internal record-keeping
    backtrackLimit: float
    queueCounter: int

    def __init__(
        self,
        params: Params,
        pointOnLinkDest: Map.PointOnLink,
    ):
        """
        This sets the parameters that are final for the entire walkPath
        algorithm execution:
        """
        self.params = params
        self.pointOnLinkDest = pointOnLinkDest

        # Record the winning queue element:
        self.winner = None

        # For tie-breaking when dealing with the priority queue.
        self.queueCounter = 0

    @dataclass(frozen=True)
    class Next:
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
            stepCount = prevStruct.stepCount + 1

        cost: float
        if incomingLink is self.pointOnLinkDest.link:
            # Last-time initialization; we have hit the destination link:
            # We are stopping midway through this link. So, subtract off the
            # distance from the end that we aren't traversing.
            linkDistPotential -= (
                1.0 - self.pointOnLinkDest.percentAlong
            ) * incomingLink.data["geometry"].length
            cost = startupCost + self.params.scoreFunction(
                self.pointOnLinkOrig, linkDistPotential, self.pointOnLinkDest
            )
        else:
            # Normal operation; we hadn't encountered the destination link yet:
            cost = startupCost + self.params.scoreFunction(
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

        destDistance: float
        cost: float
        queueCounter: int
        nextStruct: "WalkPathProcessor.Next"

    @dataclass(frozen=True)
    class PathResult:
        """
        PathResult is a container for the results of a walkPath operation.
        """

        linkList: list[Map.LinkRecord] | None
        distance: float
        cost: float
        linkListIndex: int

    def walkPath(
        self,
        pointOnLinkOrig: Map.PointOnLink,
        startupCost: float = 0.0,  # TODO: !!! Use link.data['flatScore'] !!!
        totalLinkCount: int = 0,
    ) -> PathResult:
        """
        walkPath uses a breadth-first search to find the shortest distance from a given
        PointOnLink to another PointOnLink and returns a list of links representing nodes
        and following links encountered. If nothing is found, then None is returned. An
        empty list signifies that the destination is on the same link as the origin.
        """
        # Initializations:
        self.pointOnLinkOrig = (
            pointOnLinkOrig  # TOOD: Rearrange methods to keep these local
        )
        self.winner = None

        # Are the points too far away to begin with?
        origDestDist = pointOnLinkOrig.findDistanceFrom(self.pointOnLinkDest)
        if origDestDist > self.params.limitDirectDist:
            return WalkPathProcessor.PathResult(
                linkList=None, distance=0.0, cost=0.0, linkListIndex=0
            )

        # Set a reasonable bound for the expected distance in this path search:
        self.backtrackLimit = self.params.limitPathDist

        # Set up a queue for the search. Preload the queue with the first starting location:
        self.processingQueue = queue.PriorityQueue()
        self.processingQueue.put(
            WalkPathProcessor.PathElement(
                destDistance=origDestDist,
                cost=0.0,
                queueCounter=0,
                nextStruct=self.createNext(
                    None, pointOnLinkOrig.link, startupCost, totalLinkCount
                ),
            )
        )
        self.queueCounter = 0

        # Do the breadth-first search:
        while not self.processingQueue.empty():
            self._walkPath(self.processingQueue.get().nextStruct)

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
            return WalkPathProcessor.PathResult(
                linkList=retList,
                distance=self.winner.distance,
                cost=self.winner.cost - startupCost,
                linkListIndex=self.winner.linkListIndex,
            )
        else:
            # We didn't find anything.
            return WalkPathProcessor.PathResult(
                linkList=None, distance=0.0, cost=0.0, linkListIndex=0
            )

    def _walkPath(self, walkPathElem: Next) -> None:
        """
        _walkPath is the internal processing element for the pathfinder. It is
        called internally from walkPath(), and recursively.
        """
        # Check maximum number of steps:
        if walkPathElem.stepCount >= self.params.limitSteps:
            return

        # Check total distance; we are not interested if we exceed our previous best score:
        if walkPathElem.distance >= self.backtrackLimit:
            return

        # Do we exceed the worst cost in the list of simultaneous costs?
        currentPoint = (
            Map.PointOnLink(
                walkPathElem.incomingLink,
                0,
                False,
                0,
                walkPathElem.incomingLink.getOrigin(),  # TODO: Problem. May need to use current point. But where does dist score come from?
            )
            if walkPathElem.incomingLink is not self.pointOnLinkOrig.link
            else self.pointOnLinkOrig
        )

        # Get distance from the current point to the destination point:
        crowsDistance = currentPoint.findDistanceFrom(self.pointOnLinkDest)

        # What about costs from simultaneous paths?
        if self.params.exceedsPreviousCosts(walkPathElem.cost, crowsDistance):
            return

        # Are we at one of the destinations of interest?
        if (
            self.params.globalPathCache
            and walkPathElem.incomingLink.id in self.params.allTargetLinkIDs
            and not self.params.globalPathCache.check(
                self.pointOnLinkOrig.link.id, self.pointOnLinkDest.link.id
            )
        ):
            # Log the path into the cache by looking at all of the parent elements:
            prevLinks: list[tuple[Hashable, Map.LinkRecord]] = []
            if walkPathElem.prevStruct:
                current: Map.LinkRecord = walkPathElem.incomingLink
                prior: WalkPathProcessor.Next = walkPathElem.prevStruct
                while prior.prevStruct:
                    if prior.incomingLink is not current:
                        prevLinks.append((prior.incomingLink.id, current))
                    current = prior.incomingLink
                    prior = prior.prevStruct
            self.params.globalPathCache.update(walkPathElem.incomingLink.id, prevLinks)
            # TODO: Add in third parameter to prevLinks that is "expense", so we can
            # quickly compare below to see if a path is viable to pursue.

        # Are we at our currently analyzed destination?
        if walkPathElem.incomingLink is self.pointOnLinkDest.link:
            # Log winner for this round:
            self.winner = walkPathElem
            self.backtrackLimit = walkPathElem.distance

            # Process the next queue element:
            return

        # Can we possibly get back to the destination without exceeding the limit?
        if (
            walkPathElem.distance
            + crowsDistance
            - self.pointOnLinkOrig.link.getLength()
            > self.backtrackLimit
        ):
            return

        # Look at each link that comes out from the current node. First, see
        # if there is a shortcut to our destination already in the cache:
        myList: Iterable[Map.LinkRecord]
        shortcut: Map.LinkRecord | None
        if self.params.globalPathCache:
            shortcut = self.params.globalPathCache.check(
                walkPathElem.incomingLink.id, self.pointOnLinkDest.link.id
            )
        else:
            shortcut = None
        if shortcut:
            myList = (shortcut,)
        else:
            myList = self.params.map.outgoingLinks(walkPathElem.incomingLink.destNodeID)

        link: Map.LinkRecord
        for link in myList:
            # Filter out U-turns:
            penalty = 0.0
            if (
                self.params.uTurnDeadEndPenalty or self.params.uTurnInterPenalty
            ) and self.params.map.isReverseLink(walkPathElem.incomingLink, link):
                # Is it a dead-end?
                if hasExactly(
                    self.params.map.outgoingLinks(walkPathElem.incomingLink.destNodeID),
                    1,
                ):
                    if self.params.uTurnDeadEndPenalty is None:
                        if self.params.uTurnInterPenalty is None:
                            continue
                        else:
                            penalty = self.params.uTurnInterPenalty
                    else:
                        penalty = self.params.uTurnDeadEndPenalty
                else:
                    if self.params.uTurnInterPenalty is None:
                        continue
                    else:
                        penalty = self.params.uTurnInterPenalty
                penalty = self.params.scoreFunction(
                    None, penalty, None
                )  # TODO: !!! Use stuff in link.data !!!

            # Had we visited this before?
            if link.id in walkPathElem.backtrackSet:
                continue
                # TODO: This won't work with park-and-rides where a path loops around on itself. This can possibly be fixed by adding a penalty
                # and allowing the path to be traversed. Turn this on with an option. Execution will probably be a bit slower.

            # Add to the queue for processing later:
            self.queueCounter += 1
            self.processingQueue.put(
                WalkPathProcessor.PathElement(
                    destDistance=crowsDistance,
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

        # If for some reason depth-first is worth trying, uncomment these lines:
        # while not self.processingQueue.empty():
        #    self._walkPath(self.processingQueue.get().nextStruct)
