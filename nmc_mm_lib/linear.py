"""
linear.py: Methods that have to do with geometry
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
from nmc_mm_lib import compat
import sys, math, unittest
from heapq import heappush, heappop

"Tolerance in the point-distance finder for finding ill-defined lines"
EPSILON = 0.001

def pointDistSq(pointX, pointY, lineX1, lineY1, lineX2, lineY2, norm):
    """
    pointDistSq returns the squared distance of a line segment from a point, the distance of the segment traversed,
    and also returns whether that point is perpendicularly distant from the segment (true/false).
    @type pointX: float
    @type pointY: float
    @type lineX1: float
    @type lineY1: float
    @type lineX2: float
    @type lineY2: float
    @type norm: float
    @return Tuple of the three factors
    @rtype float, float, bool
    """
    perpendicular = False
    
    # Case 1: We have an ill-defined line:
    if abs(lineX2 - lineX1) < EPSILON and abs(lineY2 - lineY1) < EPSILON:
        distSq = (pointX - lineX1) ** 2 + (pointY - lineY1) ** 2
        linkDist = 0
    
    # Case 2: We have an arbitrary line:
    else:
        # Do a vector projection:
        a = (lineX2 - lineX1)
        b = (lineY2 - lineY1)
        scale = (a * (pointX - lineX1) + b * (pointY - lineY1)) / norm ** 2
            
        # Is the projection actually on the line segment?
        if scale < 0:
            distSq = (pointX - lineX1) ** 2 + (pointY - lineY1) ** 2
            linkDist = 0
        elif scale > 1:
            distSq = (pointX - lineX2) ** 2 + (pointY - lineY2) ** 2
            linkDist = norm
        else:
            # ProjX and Y is the point on the line that's closest to the point.
            projX = lineX1 + scale * a
            projY = lineY1 + scale * b
            distSq = (pointX - projX) ** 2 + (pointY - projY) ** 2
            linkDist = scale * norm
            perpendicular = True
            
    return (distSq, linkDist, perpendicular)

def pointDist(pointX, pointY, lineX1, lineY1, lineX2, lineY2):
    """
    pointDist returns the distance of a line segment from a point, the distance of the segment traversed,
    and also returns whether that point is perpendicularly distant from the segment (true/false).
    @type pointX: float
    @type pointY: float
    @type lineX1: float
    @type lineY1: float
    @type lineX2: float
    @type lineY2: float
    @rtype float, float, bool
    """
    # Calculate the length of the segment for normalization:
    norm = getNorm(lineX1, lineY1, lineX2, lineY2)
    
    # Get the distance from the point to the segment:
    (dist, linkDist, perpendicular) = pointDistSq(pointX, pointY, lineX1, lineY1, lineX2, lineY2, norm)
    
    # Prepare the result:
    return (math.sqrt(dist), linkDist, perpendicular)

def getNorm(lineX1, lineY1, lineX2, lineY2):
    """
    getNorm is a convenience function to get the normalization factor (distance) of a line segment
    @type lineX1: float
    @type lineY1: float
    @type lineX2: float
    @type lineY2: float
    @rtype float
    """
    norm = math.sqrt((lineX2 - lineX1) ** 2 + (lineY2 - lineY1) ** 2)
    return norm

def getNormSq(lineX1, lineY1, lineX2, lineY2):
    """
    getNorm is a convenience function to get the normalization factor (distance) of a line segment
    @type lineX1: float
    @type lineY1: float
    @type lineX2: float
    @type lineY2: float
    @rtype float
    """
    norm = (lineX2 - lineX1) ** 2 + (lineY2 - lineY1) ** 2
    return norm

def lineIntersectsRectangle(pointX1, pointY1, pointX2, pointY2, uCornerX, uCornerY, lCornerX, lCornerY):
    """
    Returns true if the given line intersects the given rectangle. Lifted from 
    http://stackoverflow.com/questions/5514366/how-to-know-if-a-line-intersects-a-rectangle
    """
    # Part 1: See if beginning or end of line is in the rectangle. Then we intersect.
    if pointX1 >= uCornerX and pointX1 <= lCornerX and pointY1 >= uCornerY and pointY1 <= lCornerY \
            or pointX2 >= uCornerX and pointX2 <= lCornerX and pointY2 >= uCornerY and pointY2 <= lCornerY:
        return True
    
    # Part 2: Check if we intersect the rectangle:
    if pointX2 > pointX1:
        if uCornerX > pointX2 or lCornerX < pointX1:
            return False
    else:
        if uCornerX > pointX1 or lCornerX < pointX2:
            return False
        
    if pointY2 > pointY1:
        if uCornerY > pointY2 or lCornerY < pointY1:
            return False
    else:
        if uCornerY > pointY1 or lCornerY < pointY2:
            return False
        
    if pointX2 == pointX1:
        lineYuCornerX = pointY1
        lineYlCornerX = pointY2
    else:
        slope = (pointY2 - pointY1) / (pointX2 - pointX1)
        lineYuCornerX = slope * (uCornerX - pointX1) + pointY1
        lineYlCornerX = slope * (lCornerX - pointX1) + pointY1
    if lineYlCornerX < lineYuCornerX:
        lineYuCornerX, lineYlCornerX = (lineYlCornerX, lineYuCornerX)
    return uCornerY > lineYuCornerX and uCornerY < lineYlCornerX \
            or lCornerY > lineYuCornerX and lCornerY < lineYlCornerX \
            or uCornerY < lineYuCornerX and lCornerY > lineYlCornerX
            
class QuadSet:
    """
    A quasi-quad tree implementation to accelerate the searching of points-on-lines.
    @ivar pointLimit: The maximum allowed points within a quad before a subdivision happens.
    @ivar quadElement: A _QuadElement that is the topmost member of this quadtree
    @ivar prevQuadElement: For optimization purposes in adding lines to this quadtree, a shortcut for small links that
        are entirely contained is facilitated by storing the previously updated _QuadElement.
    """
    def __init__(self, pointLimit, uCornerX, uCornerY, lCornerX, lCornerY):
        """
        Sets up a QuadSet quadtree such that the maximum number of points stored at a given location is limited to pointLimit.
        @type pointLimit: int
        @type uCornerX: float
        @type uCornerY: float
        @type lCornerX: float
        @type lCornerY: float
        """
        self.pointLimit = pointLimit
        dimension = max(abs(lCornerX - uCornerX), abs(lCornerY - uCornerY)) / 2
        centerX = (uCornerX + lCornerX) / 2
        centerY = (uCornerY + lCornerY) / 2
        self.quadElement = _QuadElement(self, 1, centerX - dimension, centerY - dimension, centerX + dimension, centerY + dimension)
        self.prevQuadElement = None

    def storeLink(self, link):
        """
        Drills down and stores the given link in all _QuadElements that intersect it.
        @type pointX1: float
        @type pointY1: float
        @type pointX2: float
        @type pointY2: float
        @type link: graph.GraphLink
        """
        
        if link.id == 5778:
            j = 0
            j += 1
        
        
        # Do a quick test first to see if the line we're attempting to add is entirely contained by the
        # previously written quad:
        if self.prevQuadElement is not None:
            target = self.prevQuadElement 
            # Local copies for speed:
            vertices = link.vertices
            for vertIndex, vertex in enumerate(vertices[:-1]):
                if not (vertex.pointX >= self.prevQuadElement.uCornerX and vertex.pointX <= self.prevQuadElement.lCornerX \
                        and vertex.pointY >= self.prevQuadElement.uCornerY and vertex.pointY <= self.prevQuadElement.lCornerY \
                        and vertices[vertIndex + 1].pointX >= self.prevQuadElement.uCornerX and vertices[vertIndex + 1].pointX <= self.prevQuadElement.lCornerX \
                        and vertices[vertIndex + 1].pointY >= self.prevQuadElement.uCornerY and vertices[vertIndex + 1].pointY <= self.prevQuadElement.lCornerY):
                    target = self.quadElement
                    break
        else:
            target = self.quadElement
        
        # Nope, do the normal searching and storing:
        target.storeLink(link, 0)        

    def retrieveLinks(self, pointX, pointY, maxRadius=None):
        """
        Sets up a generator that returns all line objects in order of perpendicular distance from the
        point. The yield is the tuple (refDistance, minDistance, lineDist, perpendicular, link) 
        """
        maxRadiusSq = maxRadius ** 2 if maxRadius is not None else sys.float_info.max
        heap = []
        traversedSet = set()
        heappush(heap, (0.0, self.quadElement))
        
        while len(heap) > 0:
            distSq, element = heappop(heap)
            "@type element: tuple<float, boolean, graph.GraphLink>"
            if isinstance(element, _QuadElement):
                if not element.bottomLayer:
                    # Process the next layer
                    for quadElement in element.members:
                        if quadElement is not None:
                            distSq = quadElement.minimalDistanceSq(pointX, pointY)
                            if distSq <= maxRadiusSq:
                                heappush(heap, (distSq, quadElement))
                else:
                    # We got to the bottom; evaluate all links that intersect the rectangle:
                    for link in element.memberMap:
                        if link in traversedSet:
                            continue
                        distSq, lineDist, perpendicular = link.pointDistSq(pointX, pointY)
                        if distSq <= maxRadiusSq:
                            heappush(heap, (distSq, (lineDist, perpendicular, link)))
                            traversedSet.add(link)
            else: # It's a link, annotated with the stuff that was computed earlier:
                linkDist, perpendicular, link = element
                yield math.sqrt(distSq), linkDist, perpendicular, link
    
class _QuadElement:
    """
    Represents a rectangle within this QuadSet.
    @type members: list<_QuadElement>
    @type memberMap: dict<graph.GraphLink, set<int>>
    @type quadSet: QuadSet
    @type bottomLayer: bool
    @type uCornerX: float
    @type uCornerY: float
    @type lCornerX: float
    @type lCornerY: float
    """
    def __init__(self, quadSet, layer, uCornerX, uCornerY, lCornerX, lCornerY):
        """
        @type quadSet: QuadSet
        @type layer: int
        @type uCornerX: float
        @type uCornerY: float
        @type lCornerX: float
        @type lCornerY: float
        @type centerX: float
        @type centerY: float
        """
        self.members = None
        self.memberMap = {}
        self.bottomLayer = True
        self.quadSet = quadSet
        self.layer = layer
        self.uCornerX = uCornerX
        self.uCornerY = uCornerY
        self.lCornerX = lCornerX
        self.lCornerY = lCornerY
        self.centerX = (uCornerX + lCornerX) / 2
        self.centerY = (uCornerY + lCornerY) / 2
        
    def minimalDistanceSq(self, pointX, pointY, pad=0.01):
        """
        Returns the squared distance from the given point to the nearest side or corner of this rectangle, adjusted
        for padding. If the point is in the rectangle, then 0.0 is returned.
        """
        if pointX >= self.uCornerX and pointX <= self.lCornerX and pointY >= self.uCornerY and pointY <= self.lCornerY: 
            # We are inside the rectangle. Return 0.0.
            return 0.0
        else:
            edgeNorm = self.lCornerX - self.uCornerX - 2 * pad
            distSq1, _, _ = pointDistSq(pointX, pointY, self.uCornerX + pad, self.uCornerY + pad,
                self.lCornerX - pad, self.uCornerY + pad, edgeNorm)
            distSq2, _, _ = pointDistSq(pointX, pointY, self.lCornerX - pad, self.uCornerY + pad,
                self.lCornerX - pad, self.lCornerY - pad, edgeNorm)
            distSq3, _, _ = pointDistSq(pointX, pointY, self.lCornerX - pad, self.lCornerY - pad,
                self.uCornerX + pad, self.lCornerY - pad, edgeNorm)
            distSq4, _, _ = pointDistSq(pointX, pointY, self.uCornerX + pad, self.lCornerY - pad,
                self.uCornerX + pad, self.uCornerY + pad, edgeNorm)
            return min(distSq1, distSq2, distSq3, distSq4)    

    def storeLink(self, link, vertIndex):
        """
        Drills down and stores the given link in all _QuadElements that intersect it.
        @type link graph.GraphLink
        """
        if self.bottomLayer:
            if link not in self.memberMap:
                self.memberMap[link] = set()
            if len(self.memberMap) < self.quadSet.pointLimit:
                self.memberMap[link].add(vertIndex)
                self.quadSet.prevQuadElement = self
            else:
                # Promote this layer to be above the bottom:
                self._promote()
            
        if not self.bottomLayer:
            # Get local copies for optimization purposes:
            centerX, centerY = (self.centerX, self.centerY)
            uCornerX, uCornerY, lCornerX, lCornerY = (self.uCornerX, self.uCornerY, self.lCornerX, self.lCornerY)
            _prepareBelow = self._prepareBelow
            vertices = link.vertices
            pointX1, pointY1, pointX2, pointY2 = (vertices[vertIndex].pointX, vertices[vertIndex].pointY,
                vertices[vertIndex + 1].pointX, vertices[vertIndex + 1].pointY)
    
            # Check if there is intersection in any of these rectangles:
            if lineIntersectsRectangle(pointX1, pointY1, pointX2, pointY2, uCornerX, uCornerY, centerX, centerY):
                quadElement = _prepareBelow(0, uCornerX, uCornerY, centerX, centerY)
                if quadElement is not None:
                    quadElement.storeLink(link, vertIndex)
            if lineIntersectsRectangle(pointX1, pointY1, pointX2, pointY2, centerX, uCornerY, lCornerX, centerY):
                quadElement = _prepareBelow(1, centerX, uCornerY, lCornerX, centerY)
                if quadElement is not None:
                    quadElement.storeLink(link, vertIndex)
            if lineIntersectsRectangle(pointX1, pointY1, pointX2, pointY2, centerX, centerY, lCornerX, lCornerY):
                quadElement = _prepareBelow(2, centerX, centerY, lCornerX, lCornerY)
                if quadElement is not None:
                    quadElement.storeLink(link, vertIndex)
            if lineIntersectsRectangle(pointX1, pointY1, pointX2, pointY2, uCornerX, centerY, centerX, lCornerY):
                quadElement = _prepareBelow(3, uCornerX, centerY, centerX, lCornerY)
                if quadElement is not None:
                    quadElement.storeLink(link, vertIndex)

    def _prepareBelow(self, index, uCornerX, uCornerY, lCornerX, lCornerY):
        """
        Checks to see if the quad below exists, and if it doesn't, creates it.
        """
        if self.members[index] is None:
            self.members[index] = _QuadElement(self.quadSet, self.layer + 1, uCornerX, uCornerY, lCornerX, lCornerY)
        return self.members[index]
    
    def _promote(self):
        """
        Moves all member lines to new elements below.
        """
        self.members = 4 * [None]
        memberMap = self.memberMap
        self.memberMap = None
        self.bottomLayer = False
        for link, vertexSet in compat.iteritems(memberMap):
            for vertexIndex in vertexSet:
                self.storeLink(link, vertexIndex)
        
class TestLinear(unittest.TestCase):

    def test_horizontalLine(self):
        """
        Test 1: Horizontal line
        """
        self.assertEqual(pointDist(1, 3, -3, 2, 4, 2), (1, 4, True), "Line (-3, 2)-(4, 2) and Point (1, 3)") 
        self.assertEqual(pointDist(-4, 3, -3, 2, 4, 2), (math.sqrt(2), 0, False), "Line (-3, 2)-(4, 2) and Point (-4, 3)") 
        self.assertEqual(pointDist(5, 1, -3, 2, 4, 2), (math.sqrt(2), 7, False), "Line (-3, 2)-(4, 2) and Point (5, 1)") 
        self.assertEqual(pointDist(0, 2, -3, 2, 4, 2), (0, 3, True), "Line (-3, 2)-(4, 2) and Point (0, 2)") 
    
    def test_verticalLine(self):
        """
        Test 2: Vertical line
        """
        self.assertEqual(pointDist(3, 1, 2, -3, 2, 4), (1, 4, True), "Line (2, -3)-(2, 4) and Point (3, 1)")
        self.assertEqual(pointDist(3, -4, 2, -3, 2, 4), (math.sqrt(2), 0, False), "Line (2, -3)-(2, 4) and Point (3, -4)")
        self.assertEqual(pointDist(1, 5, 2, -3, 2, 4), (math.sqrt(2), 7, False), "Line (2, -3)-(2, 4) and Point (1, 5)")
        self.assertEqual(pointDist(2, 0, 2, -3, 2, 4), (0, 3, True), "Line (2, -3)-(2, 4) and Point (2, 0)")
    
    def test_arbitraryLine(self):
        """
        Test 3: Arbitrary line
        """
        self.assertEqual(pointDist(1, 1, -2, -1, 1, 2), (math.sqrt(2) / 2, math.sqrt(2) * 5 / 2, True), "Line (-2, -1)-(1, 2) and Point (1, 1)")
        
        (distRef, distLen, perpFlag) = pointDist(2, 1, -2, -1, 1, 2) 
        self.assertAlmostEqual(distRef, math.sqrt(2), 5, "Line (-2, -1)-(1, 2) and Point (2, 1) ref")
        self.assertAlmostEqual(distLen, math.sqrt(2) * 3, 5, "Line (-2, -1)-(1, 2) and Point (2, 1) len")
        # perpFlag is arbitrary in the case because of precision error.
        
        (distRef, distLen, perpFlag) = pointDist(2, 2, -2, -1, 1, 2)
        self.assertAlmostEqual(distRef, 1.0, 5, "Line (-2, -1)-(1, 2) and Point (2, 2) ref")
        self.assertAlmostEqual(distLen, math.sqrt(2) * 3, 5, "Line (-2, -1)-(1, 2) and Point (2, 2) len")
        self.assertEqual(perpFlag, False, "Line (-2, -1)-(1, 2) and Point (2, 2) perp")
        
        self.assertEqual(pointDist(-2, -2, -2, -1, 1, 2), (1, 0, False), "Line (-2, -1)-(1, 2) and Point (-2, -2)")

if __name__ == '__main__':
    unittest.main()
    