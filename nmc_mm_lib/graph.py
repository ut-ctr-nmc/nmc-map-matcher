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
from operator import index
import sys
from typing import Hashable, Iterable, NamedTuple, Sequence, MutableMapping, Generator, Any
import shapely
from shapely.ops import transform
import networkx
import pyproj

class Map:
    """
    Map is a container for a node-link graph, maintained internally in
    a Shapely SRTree for spatial indexing.
    """
    class LinkRecord(NamedTuple):
        """
        For recording individual links in edge lookups
        """
        origNodeID: Hashable
        destNodeID: Hashable
        data: MutableMapping[str, Any]

    fromCRS: pyproj.CRS
    workingCRS: pyproj.CRS
    transformer: pyproj.Transformer
    graph: networkx.DiGraph
    edgeIndexLookup: tuple[LinkRecord, ...]
    tree: shapely.strtree.STRtree
    eqCutoff: int = 2 # Decimal places for equality checks

    def __init__(self,
                 fromCRS: str = "EPSG:4326", # GPS
                 workingCRS: str = "EPSG:3857"):
        """
        Initializes an empty Map.

        @param fromCRS: The coordinate reference system of input geometries
            (Default: GPS).
        @param workingCRS: The coordinate reference system to use for internal
            calculations (Default: Web Mercator).
        """
        self.fromCRS = pyproj.CRS(fromCRS)
        self.workingCRS = pyproj.CRS(workingCRS)
        self.transformer = pyproj.Transformer.from_crs(self.fromCRS,
                                    self.workingCRS, always_xy=True)

        self.graph = networkx.DiGraph()
        self.edgeIndexLookup = tuple()
        self.tree = shapely.strtree.STRtree([])

    def addNode(self,
                nodeID: Hashable,
                lonHoriz: float, latVert: float,
                score: float = 0.0,
                metadata: dict = {}):
        """
        Adds a node to the map.

        @param nodeID: The unique identifier for the node.
        @param lonHoriz: The longitude or horizontal measure of the node.
        @param latVert: The latitude or vertical measure of the node.
        @param score: An optional score (penalty)for the node.
        @param metadata: Optional additional metadata to associate with the node.
        """
        x, y = self.transformer.transform(lonHoriz, latVert)
        self.graph.add_node(nodeID, x=x, y=y, score=score,
                            lonHoriz=lonHoriz, latVert=latVert, **metadata)

    def addLink(self,
                origNodeID: Hashable,
                destNodeID: Hashable, 
                controlPoints: Iterable[Sequence[float]] \
                    | shapely.geometry.LineString | None,
                linkID: Hashable | None = "",
                metadata: dict = {},
                hasEndpoints: bool = True,
                flatScore: float = 0.0,
                lengthWeight: float = 1.0) -> None:
        """
        Adds a directed curved link to the map using control points.

        @param origNodeID: The originating node ID.
        @param destNodeID: The destination node ID.
        @param controlPoints: A list of (lon, lat) tuples (or LineString) representing control points for the curve, or None for defaults
        @param linkID: The unique identifier for the link (optional), "" for default, or None.
        @param metadata: Optional additional metadata to associate with the link.
        @param hasEndpoints: Whether the control points include the endpoints (default: True).
        @param flatScore: An optional flat score (penalty) for the link, or 0 use the length as score (default: 0).
        @param lengthWeight: A multiplier to apply to the length when computing the score, or 0 to deactivate (default: 1.0).
        @param controlPoints: A list of (lon, lat) tuples representing control points for the curve.
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
            controlPoints = [self.transformer.transform(lon, lat) \
                             for lon, lat in lonlats]

        # Next, add endpoints if needed:
        if not hasEndpoints:
            origNode = self.graph.nodes[origNodeID]
            destNode = self.graph.nodes[destNodeID]
            controlPoints = [(origNode["x"], origNode["y"])] + controlPoints \
                + [(destNode["x"], destNode["y"])]

        # Next, convert to LineString:
        controlPoints = shapely.geometry.LineString(controlPoints)

        # Other attributes to store:
        if linkID == "":
            linkID = f"{origNodeID}->{destNodeID}"
        if linkID:
            metadata['id'] = linkID

        # TODO: Need progressive score multiplier for partial scoring!
        self.addLineStringLink(origNodeID, destNodeID, controlPoints, metadata,
                               flatScore=flatScore, lengthWeight=lengthWeight,
                               alreadyXformed=True)
        
    def addLineStringLink(self,
                          origNodeID: Hashable,
                          destNodeID: Hashable,
                          geometry: shapely.geometry.LineString,
                          metadata: dict = {},
                          flatScore: float = 0.0,
                          lengthWeight: float = 1.0,
                          alreadyXformed: bool = False) -> None:
        """
        Adds a directed link to the map.

        @param origNodeID: The originating node ID.
        @param destNodeID: The destination node ID.
        @param geometry: The geometry of the link as a LineString, already projected.
        @param metadata: Optional additional metadata to associate with the link.
        @param flatScore: An optional flat score (penalty) for the link, or 0 use the length as score (default: 0).
        @param lengthWeight: A multiplier to apply to the length when computing the score, or 0 to deactivate (default: 1.0).
        @param alreadyXformed: Whether the geometry needs to be transformed to working CRS.
        """
        if not alreadyXformed:
            geometry = transform(self.transformer.transform, geometry)
        self.graph.add_edge(origNodeID, destNodeID, geometry=geometry,
                            flatScore=flatScore, lengthWeight=lengthWeight,
                            treeIndex=-1, **metadata)

    def completeMap(self) -> None:
        """
        Completes the map by building the spatial index.
        """
        self.edgeIndexLookup = tuple(Map.LinkRecord(origNodeID=u, destNodeID=v,
                    data=data) for u, v, data in self.graph.edges(data=True))
        for index, element in enumerate(self.edgeIndexLookup):
            element.data['treeIndex'] = index
        self.tree = shapely.strtree.STRtree([element.data["geometry"] \
                                        for element in self.edgeIndexLookup])

    def isReverseLink(self,
                      linkA: LinkRecord | int,
                      linkB: LinkRecord | int) -> bool:
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
        
        if linkA.origNodeID != linkB.destNodeID \
                or linkA.destNodeID != linkB.origNodeID:
            return False

        geomA: shapely.geometry.LineString = linkA.data['geometry']
        geomB: shapely.geometry.LineString = linkB.data['geometry']
        if round(geomA.length, self.eqCutoff) != round(geomB.length, self.eqCutoff):
            return False
        tolerance = 1 / self.eqCutoff
        for (xA, yA), (xB, yB) in zip(geomA.coords, geomB.coords[::-1]):
            if abs(xA - xB) >= tolerance or abs(yA - yB) >= tolerance:
                return False
        return True

    def isSimilarLink(self,
                      linkA: LinkRecord | int,
                      linkB: LinkRecord | int) -> bool:
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
        
        if linkA.origNodeID != linkB.origNodeID \
                or linkA.destNodeID != linkB.destNodeID:
            return False

        geomA: shapely.geometry.LineString = linkA.data['geometry']
        geomB: shapely.geometry.LineString = linkB.data['geometry']
        if round(geomA.length, self.eqCutoff) != round(geomB.length, self.eqCutoff):
            return False
        tolerance = 1 / self.eqCutoff
        for (xA, yA), (xB, yB) in zip(geomA.coords, geomB.coords):
            if abs(xA - xB) >= tolerance or abs(yA - yB) >= tolerance:
                return False
        return True

    def outgoingEdges(self,
                      nodeID: Hashable) -> Generator[LinkRecord]:
        """
        Returns a generator of outgoing edges from a given node.

        @param nodeID: The node ID to look up.
        @return: A tuple of LinkRecords for outgoing edges.
        """
        return (self.edgeIndexLookup[treeIndex] \
            for u, v, treeIndex in self.graph.edges(nodeID, data="treeIndex"))

    class Trackpoint(NamedTuple):
        """
        Keeps an original track point and transformed coordinates.
        """
        lonHoriz: float
        latVert: float
        point: shapely.geometry.Point
        id: Hashable | None = None

    def makeTrackpoint(self,
                       lonHoriz: float,
                       latVert: float,
                       ident: Hashable | None = None) -> Trackpoint:
        """
        Makes a Trackpoint from horizontal/vertical coordinates.

        @param ident: The unique identifier for the point that is written to the Trackpoint object.
        @param lonHoriz: The longitude or horizontal measure of the point.
        @param latVert: The latitude or vertical measure of the point.
        @return: The Trackpoint object.
        """
        x, y = self.transformer.transform(lonHoriz, latVert)
        point = shapely.geometry.Point(x, y)
        return Map.Trackpoint(id=ident, lonHoriz=lonHoriz, latVert=latVert,
                              point=point)

    class PointOnLink(NamedTuple):
        """
        PointOnLink is a specific point on a link. This is documented in
        Figure 1 of Perrine, et al. 2015 as "point_on_link".
        """
        link: Map.LinkRecord # The link that corresponds with this PointOnLink
        distPercent: float # Percentage of distance along the link
        nonPerpPenalty: bool # "not r", True if there is to be a non-perpendicular penalty applied
        refDist: float # "d_r", the reference distance, or the working radius from the original search point
        point: shapely.geometry.Point # The point as it sits on the link

    def findPointsOnLinks(self,
                          trackPoint: Trackpoint,
                          radius: float,
                          primaryRadius: float,
                          secondaryRadius: float,
                          prevPoints: Iterable[PointOnLink],
                          limitClosestPoints: int | None = None) -> list[PointOnLink]:
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
        foundLinks = []

        # Find perpendicular and non-perpendicular PointOnLinks that are within radius.
        indices = self.tree.query(trackPoint.point, predicate="dwithin", distance=radius)
        for index in indices:
            linkRecord: Map.LinkRecord = self.edgeIndexLookup[index]
            geometry: shapely.geometry.LineString = linkRecord.data["geometry"]
            percentAlong: float = geometry.project(trackPoint.point, normalized=True)
            pointAlong: shapely.geometry.Point = geometry.interpolate(percentAlong, normalized=True)
            refDist: float = trackPoint.point.distance(pointAlong)

            # Determine if perpendicular:
            isPerpendicular = percentAlong > 0.0 and percentAlong < 1.0

            # Assume that if we are nonperpendicular with respect to the
            # end of the current link, and there are outgoing links, we may
            # in fact be perpendicular to one of those outgoing links, which
            # is more worthwhile.
            if not isPerpendicular and percentAlong >= 1.0:
                if self.graph.out_degree(linkRecord.destNodeID) > 0:
                    continue

            # A candidate:
            pointOnLink = Map.PointOnLink(linkRecord, percentAlong, not isPerpendicular, refDist, pointAlong)

            if refDist <= primaryRadius:
                # Immediate consideration if we are in the primary radius:
                foundLinks.append(pointOnLink)
            else:
                # Check to see if the point is close to a previous point. This allows candidate links to be tracked
                # that can possibly correspond with missing geometry, such as a bus going through a parking lot that
                # isn't represented in the underlying map.
                for prevPoint in prevPoints:
                    dist = pointAlong.distance(prevPoint.point)
                    if (dist < secondaryRadius):
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
    WalkPathProcessor contains methods used to conduct the walkPath algorithm.  It maintains a cache that
    persists in-between individual pathfinding operations.
    
    @ivar uTurnInterPenalty: Set this to none if U-turns are not allowed in intersections; otherwise, this
        number of feet are added at intersection U-turns.
    @type uTurnInterPenalty: float
    @ivar uTurnDeadEndPenalty: Set this to none to use the penalty value in uTurnInterEnable; otherwise,
        this number of feet are added at U-turns at dead-ends.
    @type uTurnDeadEndPenalty: float
    @ivar pathEngine: A reference to the object that instanciates this class.
    @type pathEngine: path_engine.PathEngine
    @ivar backCache: Caches previous walkPath operations to accelerate processing a little bit 
    @type backCache: dict<int, dict<int, GraphLink>>
    @ivar winner: Records the winning queue element 
    @type winner: _WalkPathNext
    @ivar processingQueue: Processing queue to facilitate the breadth-first search
    @type processingQueue: []
    @ivar pointOnLinkOrig: For internal record-keeping    
    @type pointOnLinkOrig: PointOnLink
    @ivar pointOnLinkDest: For internal record-keeping
    @type pointOnLinkDest: PointOnLink
    @ivar linkList: A list of link objects that are to be used for matching, or None if no list.
    @type linkList: list<GraphLink>
    """        
    def __init__(self, pathEngine, limitRadius, limitDistance, limitRadiusRev, limitSteps, linkList=None):
        """
        This sets the parameters that are final for the entire walkPath algorithm execution:
        @type pathEngine: path_engine.PathEngine
        @type limitRadius: float
        @type limitDistance: float
        @type limitRadiusRev: float
        @type limitSteps: int
        """
        self.pathEngine = pathEngine
        self.limitDistance = limitDistance
        self.limitRadiusRev = limitRadiusRev
        self.limitSteps = limitSteps

        self.limitRadius = limitRadius
        self.limitRadiusSq = (limitRadius ** 2) if limitRadius < sys.float_info.max else sys.float_info.max

        self.uTurnInterPenalty = None # Disable U-turns in intersections
        self.uTurnDeadEndPenalty = 50 # Allow U-turns at dead-ends
        
        # walkPath cache:
        self.backCache = {}
        
        # Keep the running score:
        self.backtrackScore = limitDistance
        
        # Record the winning queue element:
        self.winner = None
        
        # Other variables that exist throughout pathfinding iterations:
        self.processingQueue = None
        self.pointOnLinkOrig = None
        self.pointOnLinkDest = None
        
        # For tie-breaking when dealing with the priority queue.
        self.queueCounter = 0
        
        # List of required links for transit purposes.
        self.linkList = linkList
        
    class _WalkPathNext:
        """
        _WalkPathNext allows path match requests to be queued. Each of these represents a traversal from the
        start of incomingLink to the starts of the next possible links. The walkPath() method will create new
        _WalkPathNext instances for each of those possible links and enqueues them in the priority queue that
        coordinates the pathfinding operations.
        @ivar prevStruct: The previous _WalkPathNext object that represents the link traversal for the previous
            link in the path.
        @type prevStruct: _WalkPathNext
        @ivar incomingLink: The link that we are to traverse.
        @type incomingLink: GraphLink
        @ivar linkListIndex: The count of how many links have been traversed
        @type linkListIndex: int
        @ivar distance: The total distance from the origin PointOnLink to the current location.
        @type distance: float
        @ivar cost: The total calculated cost from the origin PointOnLink to the current location.
        @type cost: float
        @ivar stepCount: The number of steps traversed from the origin PointOnLink to incomingLink.
        @type stepCount: int
        @ivar backtrackSet: A set of link unique IDs for all links that had already been traversed.
        @ivar backtrackSet: set<int>
        """
        def __init__(self, processor, prevStruct, incomingLink, startupCost=0.0, linkListIndex=0):
            """
            This initializes the elements that are stored within this object.
            @type processor: WalkPathProcessor
            @type prevStruct: _WalkPathNext
            @type incomingLink: GraphLink
            @type startupCost: float
            """
            self.prevStruct = prevStruct
            self.incomingLink = incomingLink
            self.linkListIndex = linkListIndex
            
            
            if prevStruct is None:
                # First-time initialization:
                linkDistance = processor.pointOnLinkOrig.link.distance - processor.pointOnLinkOrig.dist
                self.stepCount = 0
            else:
                linkDistance = incomingLink.distance
                self.stepCount = prevStruct.stepCount + 1

            if incomingLink is processor.pointOnLinkDest.link:
                # Last-time initialization; we have hit the destination link:
                # We are stopping midway through this link.  So, subtract off the distance from the
                # end that we aren't traversing.
                linkDistance -= processor.pointOnLinkDest.link.distance - processor.pointOnLinkDest.dist
                self.cost = startupCost + processor.pathEngine.scoreFunction(processor.pointOnLinkOrig, linkDistance, processor.pointOnLinkDest)                
            else:
                # Normal operation; we hadn't encountered the destination link yet:
                self.cost = startupCost + processor.pathEngine.scoreFunction(processor.pointOnLinkOrig, linkDistance, None)
                
            self.distance = (prevStruct.distance if prevStruct else 0.0) + linkDistance

            # Make a copy of the set only if it is to change, and add in the new incoming link ID:
            oldBacktrackSet = prevStruct.backtrackSet if prevStruct is not None else set()
            "@type oldBacktrackSet: set<int>"
            if incomingLink.id not in oldBacktrackSet: 
                self.backtrackSet = set(oldBacktrackSet)
                self.backtrackSet.add(incomingLink.id)
            else:
                self.backtrackSet = oldBacktrackSet
    
    def walkPath(self, pointOnLinkOrig, pointOnLinkDest, startupCost=0.0, totalLinkCount=0):
        """
        walkPath uses a breadth-first search to find the shortest distance from a given PointOnLink to another PointOnLink and
        returns a list of links representing nodes and following links encountered.  Specify a limiting radius for
        evaluating target nodes, and maximum distance traversed.  Also specify a smaller radius for small distances backwards.
        If nothing is found, then None is returned.  An empty list signifies that the destination is on the same link as the
        origin.
        @type pointOnLinkOrig: PointOnLink
        @type pointOnLinkDest: PointOnLink
        @return List of new GraphLinks traversed, distance, and cost 
        @rtype list<GraphLink>, float, float
        """
        # Initializations:
        self.pointOnLinkOrig = pointOnLinkOrig
        self.pointOnLinkDest = pointOnLinkDest
        self.winner = None
        self.backtrackScore = self.limitDistance
        
        # Are the points too far away to begin with?
        origDestDistSq = linear.getNormSq(self.pointOnLinkDest.pointX, self.pointOnLinkDest.pointY, self.pointOnLinkOrig.pointX, self.pointOnLinkOrig.pointY)
        if origDestDistSq > self.limitRadiusSq:
            return None, 0.0, 0.0, 0
        
        # Set a reasonable bound for the expected distance in this path search:
        self.backtrackScore = self.limitDistance

        # Set up a queue for the search.  Preload the queue with the first starting location:
        self.processingQueue = []
        heappush(self.processingQueue, (0.0, 0, self._WalkPathNext(self, None, self.pointOnLinkOrig.link, startupCost, totalLinkCount)))
        self.queueCounter = 0
        
        # Do the breadth-first search:
        while self.processingQueue:
            self._walkPath(heappop(self.processingQueue)[-1])
  
        # Set up the return:
        if self.winner is not None:
            # Iterate through all of the links we have traversed. (Ignore first item because we
            # hadn't technically traversed it).
            retList = []
            "@type retList: list<GraphLink>"
            element = self.winner
            "@type element: _WalkPathNext"
            while element.prevStruct is not None:
                retList.append(element.incomingLink)
                element = element.prevStruct
            retList.reverse()
            return retList, self.winner.distance, self.winner.cost - startupCost, self.winner.linkListIndex
        else:
            # We didn't find anything.
            return None, 0.0, 0.0, 0
        
    # _walkPath is called internally by walkPath().
    def _walkPath(self, walkPathElem):
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
                    if (element.prevStruct.incomingLink.id in mappings) \
                            and (mappings[element.prevStruct.incomingLink.id] is element.incomingLink):
                        break
                    mappings[element.prevStruct.incomingLink.id] = element.incomingLink
                    element = element.prevStruct
                
            # Process the next queue element:
            return
        
        # Look at each link that comes out from the current node.
        # First, see if there is a shortcut to our destination already in the cache:
        if (self.pointOnLinkDest.link.id in self.backCache) and \
                (walkPathElem.incomingLink.id in self.backCache[self.pointOnLinkDest.link.id]):
            myList = [self.backCache[self.pointOnLinkDest.link.id][walkPathElem.incomingLink.id]]
        else:
            myList = walkPathElem.incomingLink.destNode.outgoingLinkMap.values()
        for link in myList:
            # Filter out U-turns:
            penalty = 0.0            
            if (self.uTurnDeadEndPenalty != 0 or self.uTurnInterPenalty != 0) and walkPathElem.incomingLink.isComplementary(link):
                # Is it a dead-end?
                if len(walkPathElem.incomingLink.destNode.outgoingLinkMap) == 1:
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
                penalty = self.pathEngine.scoreFunction(None, penalty, None)
                
            # Is this the next link we need to process according to the link list (transit)?
            if self.linkList and walkPathElem.linkListIndex + 1 < len(self.linkList) and self.linkList[walkPathElem.linkListIndex + 1].id != link.id:
                continue
                # TODO: We want to eventually allow the path to be departed and then regained. How to do this? We can create a "path lost" state,
                # and as long as that state is True, then search forward in self.linkList to see if we regain the path. Or, add the indices into
                # the link list into the actual link graph objects (as sets). Departure from a set will incur a penalty, and encounter with a set
                # will allow the index to be reset to the last known value.  
                                    
            # Had we visited this before?
            if link.id in walkPathElem.backtrackSet:
                continue
                # TODO: This won't work with park-and-rides where a path loops around on itself. This can possibly be fixed by adding a penalty
                # and allowing the path to be traversed. Turn this on with an option. Execution will probably be a bit slower.
            
            # Add to the queue for processing later:
            self.queueCounter += 1
            heappush(self.processingQueue, (walkPathElem.cost + penalty, self.queueCounter, self._WalkPathNext(self, walkPathElem, link, walkPathElem.cost + penalty,
                walkPathElem.linkListIndex + 1)))
