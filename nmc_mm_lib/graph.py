"""
graph.py: Links and nodes for graph models; also a breadth-first search.
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
from nmc_mm_lib import linear, gps
import sys, math
from heapq import heappush, heappop

"The maximum number of lines allowed in a quad for optimized QuadSet lookup."
DEFAULT_QUAD_LIMIT = 200

"Tolerance in the closest-point finder for dealing with corner arcs."
EPSILON = 0.001

class GraphLink:
    """
    GraphLink is a link that connects one node to another.
    @ivar id: ?
    @ivar origNode: GraphNode
    @ivar destNode: GraphNode
    @ivar vertices: list<GraphLinkVertex>
    @ivar distance: float
    @ivar graphLib: GraphLib
    """
    def __init__(self, ident, origNode, destNode, graphLib):
        """
        @type ident: int
        @type origNode: GraphNode
        @type destNode: GraphNode 
        """
        self.id = ident
        self.origNode = origNode
        self.destNode = destNode
        
        # The head of the linked-list of member vertices within this Link
        self.vertices = None
        
        # Placeholders for efficiency of measurement:
        self.distance = 0.0
        
        # Reference to coordinate system and other links:
        self.graphLib = graphLib
        
    def makeVertices(self):
        """
        For links that are straight and have no curvature, creates a pair of vertices that correspond
        with the given link, drawing a line from the origin node to the destination node.
        """
        linkVertices = [GraphLinkVertex(self.origNode.gpsLat, self.origNode.gpsLng), 
                        GraphLinkVertex(self.destNode.gpsLat, self.destNode.gpsLng)]
        self.addVertices(linkVertices)
        
    def addVertices(self, linkVertices):
        """
        Adds the list of vertices to the given link, computing coordinates and distances along the way.
        @type linkVertices: list<GraphLinkVertex>
        """
        """
        # TEST!
        newList = [None] * 4
        newList[0] = linkVertices[0]
        newList[1] = GraphLinkVertex((linkVertices[1].lat - linkVertices[0].lat) / 3.0 + linkVertices[0].lat,
                                     (linkVertices[1].lng - linkVertices[0].lng) / 3.0 + linkVertices[0].lng)
        newList[2] = GraphLinkVertex((linkVertices[1].lat - linkVertices[0].lat) * 2 / 3.0 + linkVertices[0].lat,
                                     (linkVertices[1].lng - linkVertices[0].lng) * 2 / 3.0 + linkVertices[0].lng)
        newList[3] = linkVertices[1]
        linkVertices = newList
        """
        
        gps2feet = self.graphLib.gps.gps2feet
        prevVertex = linkVertices[0]
        prevVertex.pointX, prevVertex.pointY = gps2feet(prevVertex.lat, prevVertex.lng)
        prevVertex.distance = 0.0
        self.vertices = linkVertices
        for nextVertex in linkVertices[1:]:
            nextVertex.pointX, nextVertex.pointY = gps2feet(nextVertex.lat, nextVertex.lng)
            nextVertex.distance = prevVertex.distance + linear.getNorm(prevVertex.pointX, prevVertex.pointY,
                nextVertex.pointX, nextVertex.pointY)
            prevVertex = nextVertex            
            
    def pointDistSq(self, pointX, pointY):
        """
        Finds the minimal point distance of all of the link segments between the vertices.
        @return Tuple of the minimum squared distance of a line segment from the point, the distance of the link traversed,
            and also returns whether that point is effectively perpendicular from the entire segment.
        @rtype float, float, bool
        """
        minDistSq = sys.float_info.max
        minLinkDist = sys.float_info.max
        minPerpendicular = False
        minIndex = -1
        prevVertex = self.vertices[0]
        for prevIndex, nextVertex in enumerate(self.vertices[1:]):
            distSq, linkDist, perpendicular = linear.pointDistSq(pointX, pointY, prevVertex.pointX, prevVertex.pointY,
                nextVertex.pointX, nextVertex.pointY, nextVertex.distance - prevVertex.distance)
            if distSq < minDistSq - EPSILON:
                minDistSq, minLinkDist, minPerpendicular = distSq, linkDist, perpendicular
                minIndex = prevIndex
            prevVertex = nextVertex
        
        # If we are nonperpendicular, see if the point falls within the arc that sits between the neighboring
        # segments.
        if not minPerpendicular:
            if (minIndex > 0 or minIndex == 0 and minLinkDist > EPSILON) and (minIndex < len(self.vertices) - 1 \
                    or minIndex == len(self.vertices) - 1 and minLinkDist + self.vertices[-2].distance < self.vertices[-1].distance - EPSILON): 
                minPerpendicular = True
                
        """
        # TEST!
        print("md: %g; ld: %g; p: %d" % (math.sqrt(minDistSq), minLinkDist + self.vertices[minIndex].distance, minPerpendicular))
        """

        return minDistSq, minLinkDist + self.vertices[minIndex].distance, minPerpendicular 
            
    def isComplementary(self, otherLink):
        """
        Return True if the given link directly flows in the opposite direction of this link.
        @type otherLink: GraphLink
        """
        # Test 1: Check that the nodes are shared in a complementary way:
        if not (otherLink.destNode is self.origNode and otherLink.origNode is self.destNode \
                and self.origNode is not self.destNode):
            return False

        # Test 2: Compare vertices in between the nodes:
        index = 1
        limit = math.ceil(len(self.vertices) / 2)
        while index < limit:
            if abs(self.vertices[-(index + 2)].pointX - self.vertices[index].pointX) >= EPSILON \
                    or abs(self.vertices[-(index + 2)].pointY - self.vertices[index].pointY) >= EPSILON:
                return False
            index += 1
        
        return True
    
    def __lt__(self, other):
        return self.id < other.id
    
class GraphLinkVertex:
    """
    GraphLinkVertex is a vertex possibly among several for a GraphLink. Use GraphLink.addVertices()
    to associate this with the parent link to compute distance and X-Y coordinates.
    @ivar lat: float
    @ivar lng: float
    @ivar pointX: float
    @ivar pointY: float
    @ivar distance: float
    """
    def __init__(self, lat, lng):
        """
        Sets the next vertex to allow for a segment to be expressed between this object and nextVertex.
        @type gpsLat: float
        @type gpsLng: float
        """
        self.lat = lat
        self.lng = lng
        self.pointX = 0.0
        self.pointY = 0.0
        self.distance = 0.0

class GraphNode:
    """
    GraphNode is a node that connects multiple links together.
    """
    def __init__(self, ident, gpsLat, gpsLng):
        """
        @type ident: int
        @type gpsLat: float
        @type gpsLng: float
        """
        self.id = ident
        self.gpsLat = gpsLat
        self.gpsLng = gpsLng
        self.outgoingLinkMap = {}
        "@type self.outgoingLinkMap: dict<int, GraphLink>" 
        
        # Placeholders for coordinates in feet; these won't be filled out until this is added to the GraphLib.
        self.coordX = 0.0
        self.coordY = 0.0

class PointOnLink:
    """
    PointOnLink is a specific point on a link.  This is documented in Figure 1 of Perrine, et al. 2015
    as "point_on_link".
    
    @ivar link: "L", the link that corresponds with this PointOnLink
    @type link: GraphLink
    @ivar dist: "d", the distance along the link from the origin in feet
    @type dist: float
    @ivar nonPerpPenalty: "not r", true if there is to be a non-perpendicular penalty applied
    @type nonPerpPenalty: bool
    @ivar refDist: "d_r", the reference distance, or the working radius from the original search point
    @type refDist: float
    @ivar pointX: The point x-coordinate
    @type pointX: float
    @ivar pointY: The point y-coordinate
    @type pointY: float
    """
    def __init__(self, link, dist, nonPerpPenalty = False, refDist = 0.0):
        """
        @type link: GraphLink
        @type dist: float
        @type nonPerpPenalty: bool
        @type refDist: float
        """
        self.link = link
        self.dist = dist
        self.nonPerpPenalty = nonPerpPenalty
        self.refDist = refDist
        
        if link is not None:
            if link.distance == 0:
                dist = 0
                norm = 1
            else:
                norm = link.distance
            
            self.pointX = link.origNode.coordX + (link.destNode.coordX - link.origNode.coordX) * dist / norm
            self.pointY = link.origNode.coordY + (link.destNode.coordY - link.origNode.coordY) * dist / norm
        else:
            self.pointX = 0
            self.pointY = 0
                
class GraphLib:
    """
    GraphLib is the container that holds an entire graph.
    
    @ivar gps: Reference GPS center coordinates plus calculator
    @type gps: gps.GPS
    @ivar nodeMap: Collection of nodes that are in this graph.
    @type nodeMap: dict<int, GraphNode>
    @ivar linkMap: Collection of links that are in this graph.
    @type linkMap: dict<int, GraphLink>
    @ivar prevLinkID: Previous link ID for cases where we are dealing with single-paths
    @ivar quadLimit: The maximum number of points allowed at a QuadSet layer.
    @ivar quadSet: The linear.QuadSet object that assists in finding lines of closest perpendicular distances
    """
    def __init__(self, gpsCtrLat, gpsCtrLng, quadLimit=DEFAULT_QUAD_LIMIT):
        """
        @type gpsCtrLat: float
        @type gpsCtrLng: float
        """
        self.gps = gps.GPS(gpsCtrLat, gpsCtrLng)
        self.nodeMap = {}
        self.linkMap = {}
        self.prevLinkID = 0
        self.quadLimit = quadLimit
        self.quadSet = None

    def addNode(self, node):
        """
        addNode adds a node to the GraphLib and translates its coordinates to feet. Not supported for single-path.
        @type node: GraphNode
        """
        node.coordX, node.coordY = self.gps.gps2feet(node.gpsLat, node.gpsLng)
        self.nodeMap[node.id] = node
        
    def addLink(self, link):
        """
        addLink adds a link to the GraphLib and updates its respective nodes.  Call 
        addNode first.
        @type link: GraphLink
        """
        if link.origNode.id not in self.nodeMap:
            print('WARNING: Node %d is not present.' % link.origNode.id, file = sys.stderr)
            return
        link.distance = linear.getNorm(link.origNode.coordX, link.origNode.coordY, link.destNode.coordX, link.destNode.coordY)
        ourID = link.id
        self.prevLinkID = link.id
        self.linkMap[ourID] = link
        link.origNode.outgoingLinkMap[link.id] = link
            
    def generateQuadSet(self):
        """
        This performs the task of generating the quadtree for this GraphLib. Calls to GraphLink.addVertices()
        should have been made, or if there are no vertices, makeVertices() will be called to create straight
        segments between nodes.
        """
        minX = sys.float_info.max
        minY = sys.float_info.max
        maxX = sys.float_info.min
        maxY = sys.float_info.min
        for link in self.linkMap.values():
            if not link.vertices:
                link.makeVertices()
            for vertex in link.vertices:
                minX = min(minX, vertex.pointX)
                minY = min(minY, vertex.pointY)
                maxX = max(maxX, vertex.pointX)
                maxY = max(maxY, vertex.pointY)
        
        self.quadSet = linear.QuadSet(self.quadLimit, minX, minY, maxX, maxY)
        for link in self.linkMap.values():
            self.quadSet.storeLink(link)
        
    def findPointsOnLinks(self, pointX, pointY, radius, primaryRadius, secondaryRadius, prevPoints, limitClosestPoints=sys.maxsize):
        """
        findPointsOnLinks searches through the graph and finds all PointOnLinks that are within the radius.
        Then, eligible links are proposed primaryRadius distance around the GTFS point, or secondaryRadius
        distance from the previous VISTA points.  Returns an empty list if none are found.  This corresponds
        with algorithm "FindPointsOnLinks" in Figure 1 of Perrine, et al. 2015. This expects that
        generateQuadSet has already been run.
        @type pointX: float
        @type pointY: float
        @type radius: float
        @type primaryRadius: float
        @type secondaryRadius: float
        @type prevPoints: list<PointOnLink>
        @type limitClosestPoints: int
        @rtype list<PointOnLink>
        """
        retList = []
        secondaryRadiusSq = secondaryRadius ** 2

        # Find perpendicular and non-perpendicular PointOnLinks that are within radius.
        for refDist, linkDist, perpendicular, link in self.quadSet.retrieveLinks(pointX, pointY, radius):
            # Everything coming back from retrieveLinks is sorted according to the distance from point to
            # line, and is limited to the given radius. Are we done?
            if len(retList) >= limitClosestPoints:
                break
            
            # Filter out duplicate locations represented by a nonperpendicular match to the end of one link and a
            # nonperpendicular match to the start of the following link. Keep the downstream one:                
            if not perpendicular and linkDist > 0 and len(link.destNode.outgoingLinkMap) > 0:
                continue
            
            """
            # TEST!
            print("id: %d, ld: %g, rd: %g" % (link.id, linkDist, refDist))
            """
            
            # Here is a candidate.
            pointOnLink = PointOnLink(link, linkDist, not perpendicular, refDist)
            if refDist <= primaryRadius:
                retList.append(pointOnLink)
            else:
                # Check to see if the point is close to a previous point. This allows candidate links to be tracked
                # that can possibly correspond with missing geometry, such as a bus going through a parking lot that
                # isn't represented in the underlying map.
                for prevPoint in prevPoints:
                    "@type prevPoint: PointOnLink"
                    distSq = linear.getNormSq(pointOnLink.pointX, pointOnLink.pointY, prevPoint.pointX, prevPoint.pointY)
                    if (distSq < secondaryRadiusSq):
                        # We have a winner:
                        retList.append(pointOnLink)
                        break

        # Return the limitClosestPoints number of points: 
        return retList

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
    """        
    def __init__(self, pathEngine, limitRadius, limitDistance, limitRadiusRev, limitSteps):
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
        self.uTurnDeadEndPenalty = None # Disable U-turns at dead-ends
        
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
        @ivar distance: The total distance from the origin PointOnLink to the current location.
        @type distance: float
        @ivar cost: The total calculated cost from the origin PointOnLink to the current location.
        @type cost: float
        @ivar stepCount: The number of steps traversed from the origin PointOnLink to incomingLink.
        @type stepCount: int
        @ivar backtrackSet: A set of link unique IDs for all links that had already been traversed.
        @ivar backtrackSet: set<int>
        """
        def __init__(self, processor, prevStruct, incomingLink, startupCost=0.0):
            """
            This initializes the elements that are stored within this object.
            @type processor: WalkPathProcessor
            @type prevStruct: _WalkPathNext
            @type incomingLink: GraphLink
            @type startupCost: float
            """
            self.prevStruct = prevStruct
            self.incomingLink = incomingLink
            
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
                
                """
                # TEST!                
                print("B: pil=%d; il=%d; sc=%g; d=%g; c=%g" % (prevStruct.incomingLink.id if prevStruct else -1, incomingLink.id, startupCost,
                                                               linkDistance, self.cost))
                """
                
            else:
                # Normal operation; we hadn't encountered the destination link yet:
                self.cost = startupCost + processor.pathEngine.scoreFunction(processor.pointOnLinkOrig, linkDistance, None)
                
                """
                # TEST!                
                print("A: pil=%d; il=%d; sc=%g; d=%g; c=%g" % (prevStruct.incomingLink.id if prevStruct else -1, incomingLink.id, startupCost,
                                                               linkDistance, self.cost))
                """
                
            self.distance = prevStruct.distance if prevStruct else 0.0 + linkDistance

            # Make a copy of the set only if it is to change, and add in the new incoming link ID:
            oldBacktrackSet = prevStruct.backtrackSet if prevStruct is not None else set()
            "@type oldBacktrackSet: set<int>"
            if incomingLink.id not in oldBacktrackSet: 
                self.backtrackSet = set(oldBacktrackSet)
                self.backtrackSet.add(incomingLink.id)
            else:
                self.backtrackSet = oldBacktrackSet
    
    def walkPath(self, pointOnLinkOrig, pointOnLinkDest, startupCost=0.0):
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
            return (None, 0.0, 0.0)
        
        # Set a reasonable bound for the expected distance in this path search:
        self.backtrackScore = self.limitDistance

        # Set up a queue for the search.  Preload the queue with the first starting location:
        self.processingQueue = []
        heappush(self.processingQueue, (0.0, 0, self._WalkPathNext(self, None, self.pointOnLinkOrig.link, startupCost)))
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
            return (retList, self.winner.distance, self.winner.cost - startupCost)
        else:
            # We didn't find anything.
            return (None, 0.0, 0.0)
        
    # _walkPath is called internally by walkPath().
    def _walkPath(self, walkPathElem):
        """
        _walkPath is the internal processing element for the pathfinder.
        @type walkPathElem: _WalkPathNext
        """
        """
        # TEST!
        if self.pointOnLinkOrig.link.id == 31 and self.pointOnLinkDest.link.id == 31:
            j = 0
            j += 1
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
            if (self.uTurnDeadEndPenalty != 0 or self.uTurnDeadEndPenalty != 0) and walkPathElem.incomingLink.isComplementary(link):
                # Is it a dead-end?
                if len(link.destNode.outgoingLinkMap) == 1:
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
                                    
            # Had we visited this before?
            if link.id in walkPathElem.backtrackSet:
                continue
            
            # Add to the queue for processing later:
            self.queueCounter += 1
            heappush(self.processingQueue, (walkPathElem.cost + penalty, self.queueCounter, self._WalkPathNext(self, walkPathElem, link, walkPathElem.cost + penalty)))
