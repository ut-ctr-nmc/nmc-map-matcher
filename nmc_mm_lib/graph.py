"""
graph.py: Links and nodes for graph models; also a rudimentary breadth-first search.
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
from collections import deque
import linear, gps, sys, math, operator

class GraphLink:
    """
    GraphLink is a link that connects one node to another.
    """
    def __init__(self, ident, origNode, destNode):
        """
        @type ident: int
        @type origNode: GraphNode
        @type destNode: GraphNode 
        """
        self.id = ident
        self.origNode = origNode
        self.destNode = destNode
        
        # Placeholders for efficiency of measurement:
        self.distance = 0.0
        
    def isComplementary(self, otherLink):
        """
        Return True if the given link directly flows in the opposite direction of this link.
        @type otherLink: GraphLink
        """
        return otherLink.destNode is self.origNode and otherLink.origNode is self.destNode

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
    @ivar nodeMap: Collection of nodes
    @type nodeMap: dict<int, GraphNode>
    @ivar linkMap: Collection of links
    @type linkMap: dict<int, GraphLink>
    """
    def __init__(self, gpsCtrLat, gpsCtrLng):
        """
        @type gpsCtrLat: float
        @type gpsCtrLng: float
        """
        self.gps = gps.GPS(gpsCtrLat, gpsCtrLng)
        self.nodeMap = {}
        self.linkMap = {}

    def addNode(self, node):
        """
        addNode adds a node to the GraphLib and translates its coordinates to feet.
        @type node: GraphNode
        """
        (node.coordX, node.coordY) = self.gps.gps2feet(node.gpsLat, node.gpsLng)
        self.nodeMap[node.id] = node
        
    def addLink(self, link):
        """
        addLink adds a link to the GraphLib and updates its respective nodes.  Call 
        addNode first.
        @type link: GraphLink
        """
        if link.origNode.id not in self.nodeMap:
            print('WARNING: Node %d is not present.' % link.origNode, file = sys.stderr)
            return
        link.distance = linear.getNorm(link.origNode.coordX, link.origNode.coordY, link.destNode.coordX, link.destNode.coordY)
        self.linkMap[link.id] = link
        self.nodeMap[link.origNode.id].outgoingLinkMap[link.id] = link
        
    def findPointsOnLinks(self, pointX, pointY, radius, primaryRadius, secondaryRadius, prevPoints, limitClosestPoints = sys.maxint):
        """
        findPointsOnLinks searches through the graph and finds all PointOnLinks that are within the radius.
        Then, eligible links are proposed primaryRadius distance around the GTFS point, or secondaryRadius
        distance from the previous VISTA points.  Returns an empty list if none are found.  This corresponds
        with algorithm "FindPointsOnLinks" in Figure 1 of Perrine, et al. 2015.
        @type pointX: float
        @type pointY: float
        @type radius: float
        @type primaryRadius: float
        @type secondaryRadius: float
        @type prevPoints: list<PointOnLink>
        @type limitClosestPoints: int
        @rtype list<PointOnLink>
        """
        # TODO: This brute-force implementation can be more efficient with quad trees, etc. rather than
        # scanning through all elements.
        radiusSq = radius ** 2
        primaryRadiusSq = primaryRadius ** 2
        secondaryRadiusSq = secondaryRadius ** 2
        retSet = set()
        
        # Find perpendicular and non-perpendicular PointOnLinks that are within radius.
        for link in self.linkMap.values():
            "@type link: graph.GraphLink"
            (distSq, linkDist, perpendicular) = linear.pointDistSq(pointX, pointY, link.origNode.coordX, link.origNode.coordY,
                                                                   link.destNode.coordX, link.destNode.coordY, link.distance)
            if distSq <= radiusSq:
                pointOnLink = PointOnLink(link, linkDist, not perpendicular, math.sqrt(distSq))
                
                # We are within the initial search radius.  Are we then within the primary radius?
                if distSq <= primaryRadiusSq:
                    # Yes, easy.  Add to the list:
                    retSet.add(pointOnLink)
                else:
                    # Check to see if the point is close to a previous point:
                    for prevPoint in prevPoints:
                        "@type prevPoint: PointOnLink"
                        distSq = linear.getNormSq(pointOnLink.pointX, pointOnLink.pointY, prevPoint.pointX, prevPoint.pointY)
                        if (distSq < secondaryRadiusSq):
                            # We have a winner:
                            retSet.add(pointOnLink)
                            break
                    
        ret = list(retSet)
        
        # TODO: If there is a nonperpendicular link and distance = 0, and there also exists in the set a link
        # that leads to the first link's parent node, then get rid of that first link.
        
        # Keep limited number of closest values 
        ret.sort(key = operator.attrgetter('refDist'))
        return ret[0:limitClosestPoints]

class WalkPathProcessor:
    """
    WalkPathProcessor contains methods used to conduct the walkPath algorithm.  It maintains a cache that
    persists in-between individual pathfinding operations.
    
    @ivar backCache: Caches previous walkPath operations to accelerate processing a little bit 
    @type backCache: dict<int, dict<int, GraphLink>>
    @ivar winner: Records the winning queue element 
    @type winner: _WalkPathNext
    @ivar processingQueue: Processing queue to facilitate the breadth-first search
    @type processingQueue: deque
    @ivar pointOnLinkOrig: For internal record-keeping    
    @type pointOnLinkOrig: PointOnLink
    @ivar pointOnLinkDest: For internal record-keeping
    @type pointOnLinkDest: PointOnLink
    """        
    def __init__(self, limitRadius, limitDistance, limitRadiusRev, limitSteps):
        """
        This sets the parameters that are final for the entire walkPath algorithm execution:
        @type limitRadius: float
        @type limitDistance: float
        @type limitRadiusRev: float
        @type limitSteps: int
        """
        self.limitDistance = limitDistance
        self.limitRadiusRev = limitRadiusRev
        self.limitSteps = limitSteps

        self.limitRadius = limitRadius
        self.limitRadiusSq = (limitRadius ** 2) if limitRadius < sys.float_info.max else sys.float_info.max

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
        
    class _WalkPathNext:
        """
        _WalkPathNext allows walk_path requests to be queued.
        """
        def __init__(self, processor, prevStruct, incomingLink):
            """
            This initializes the elements that are stored within this object.
            @type processor: WalkPathProcessor
            @type prevStruct: _WalkPathNext
            @type incomingLink: GraphLink
            """
            self.prevStruct = prevStruct
            self.incomingLink = incomingLink
            if prevStruct is None:
                # First-time initialization:
                self.distance = processor.pointOnLinkOrig.link.distance - processor.pointOnLinkOrig.dist
                self.stepCount = 0
            else:
                self.distance = prevStruct.distance + incomingLink.distance
                self.stepCount = prevStruct.stepCount + 1
                
            # Last-time initialization; we have hit the destination link:
            if incomingLink is processor.pointOnLinkDest.link:
                # We are stopping midway through this link.  So, subtract off the distance from the
                # end that we aren't traversing.
                self.distance -= processor.pointOnLinkDest.link.distance - processor.pointOnLinkDest.dist
            
            # Make a copy of the set only if it is to change, and add in the new incoming link ID:
            oldBacktrackSet = prevStruct.backtrackSet if prevStruct is not None else set()
            "@type oldBacktrackSet: set<int>"
            if incomingLink.id not in oldBacktrackSet: 
                self.backtrackSet = set(oldBacktrackSet)
                self.backtrackSet.add(incomingLink.id)
            else:
                self.backtrackSet = oldBacktrackSet
    
    def walkPath(self, pointOnLinkOrig, pointOnLinkDest):
        """
        walkPath uses a breadth-first search to find the shortest distance from a given PointOnLink to another PointOnLink and
        returns a list of links representing nodes and following links encountered.  Specify a limiting radius for
        evaluating target nodes, and maximum distance traversed.  Also specify a smaller radius for small distances backwards.
        If nothing is found, then None is returned.  An empty list signifies that the destination is on the same link as the
        origin.
        @type pointOnLinkOrig: PointOnLink
        @type pointOnLinkDest: PointOnLink
        @rtype list<GraphLink>, float
        """
        # Initializations:
        self.pointOnLinkOrig = pointOnLinkOrig
        self.pointOnLinkDest = pointOnLinkDest
        self.winner = None
        self.backtrackScore = self.limitDistance
        
        # Are the points too far away to begin with?
        if (self.pointOnLinkDest.pointX - self.pointOnLinkOrig.pointX) ** 2 \
                + (self.pointOnLinkDest.pointY - self.pointOnLinkOrig.pointY) ** 2 > self.limitRadiusSq:
            return (None, 0)

        # Set up a queue for the search.  Preload the queue with the first starting location:
        self.processingQueue = deque()
        self.processingQueue.append(self._walk_path_next(self, None, self.pointOnLinkOrig.link))
        
        # Do the breadth-first search:
        while len(self.processingQueue) > 0:
            self._walkPath(self.processingQueue.popleft())
        
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
            return (retList, self.winner.distance)
        else:
            # We didn't find anything.
            return (None, 0)
        
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
            """
            # Filter out U-turns:
            if walkPathElem.incomingLink.isComplementary(link):
                continue
            """
            
            # Had we visited this before?
            if link.id in walkPathElem.backtrackSet:
                continue
            
            # Add to the queue for processing later:
            self.processingQueue.append(self._WalkPathNext(self, walkPathElem, link))
