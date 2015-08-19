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
import sys, math, unittest
from collections import deque
from heapq import heappush, heappop

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
    @rtype float, float, bool
    """
    perpendicular = False
    
    # Case 1: We have an ill-defined line:
    if lineX1 == lineX2 and lineY1 == lineY2:
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

"SQRT2 = math.sqrt(2)"

def _calcY(xVal, pointX1, pointY1, pointX2, pointY2):
    if pointX2 == pointX1:
        return None
    return pointY1 + (xVal - pointX1) * (pointY2 - pointY1) / (pointX2 - pointX1)

def _calcX(yVal, pointX1, pointY1, pointX2, pointY2):
    if pointY2 == pointY1:
        return None
    return pointX1 + (yVal - pointY1) * (pointX2 - pointX1) / (pointY2 - pointY1)

def lineIntersectsRectangle(pointX1, pointY1, pointX2, pointY2, uCornerX, uCornerY, lCornerX, lCornerY):
    """
    Returns true if the given line intersects the given rectangle. Lifted from 
    http://stackoverflow.com/questions/1354472/detect-if-line-segment-intersects-square
    """
    # Part 1: See if beginning or end of line is in the rectangle. Then we intersect.
    if pointInRectangle(pointX1, pointY1, uCornerX, uCornerY, lCornerX, lCornerY) \
            or pointInRectangle(pointX2, pointY2, uCornerX, uCornerY, lCornerX, lCornerY):
        return True
    
    # Part 2: Check if we intersect the rectangle:
    if lineIntersectsLine(pointX1, pointY1, pointX2, pointY2, uCornerX, uCornerY, lCornerX, uCornerY) \
            or lineIntersectsLine(pointX1, pointY1, pointX2, pointY2, lCornerX, uCornerY, lCornerX, lCornerY) \
            or lineIntersectsLine(pointX1, pointY1, pointX2, pointY2, lCornerX, lCornerY, uCornerX, lCornerY) \
            or lineIntersectsLine(pointX1, pointY1, pointX2, pointY2, uCornerX, lCornerY, uCornerX, uCornerY):
        return True
    """
    lY = _calcX(lCornerY, pointX1, pointY1, pointX2, pointY2)
    uY = _calcX(uCornerY, pointX1, pointY1, pointX2, pointY2)
    uX = _calcY(uCornerX, pointX1, pointY1, pointX2, pointY2)
    lX = _calcY(lCornerX, pointX1, pointY1, pointX2, pointY2)    
    return lY is not None and lY < lCornerX and lY >= uCornerX \
        or uY is not None and uY < lCornerX and uY >= uCornerX \
        or uX is not None and uX < lCornerY and uX >= uCornerY \
        or lX is not None and lX < lCornerY and lX >= uCornerY
    """
    
    
    """
    # Part 1: See if beginning or end of line is in the rectangle. Then we intersect.
    if pointInRectangle(pointX1, pointY1, uCornerX, uCornerY, lCornerX, lCornerY) \
            or pointInRectangle(pointX2, pointY2, uCornerX, uCornerY, lCornerX, lCornerY):
        return True
    
    # Part 2: Check if we intersect the rectangle:
    """

    """
        # Now check for line-to-rectangle intersection 
        if lineIntersectsLine(pointX1, pointY1, pointX2, pointY2, uCornerX, uCornerY, lCornerX, uCornerY) \
                or lineIntersectsLine(pointX1, pointY1, pointX2, pointY2, lCornerX, uCornerY, lCornerX, lCornerY) \
                or lineIntersectsLine(pointX1, pointY1, pointX2, pointY2, lCornerX, lCornerY, uCornerX, lCornerY) \
                or lineIntersectsLine(pointX1, pointY1, pointX2, pointY2, uCornerX, lCornerY, uCornerX, uCornerY):
            return True
    """
        
    "return False"
        
def lineIntersectsLine(pointAX1, pointAY1, pointAX2, pointAY2, pointBX1, pointBY1, pointBX2, pointBY2):
    """
    Returns true if the two given lines intersects. Lifted from
    http://stackoverflow.com/questions/5514366/how-to-know-if-a-line-intersects-a-rectangle
    """ 
    q = (pointAY1 - pointBY1) * (pointBX2 - pointBX1) - (pointAX1 - pointBX1) * (pointBY2 - pointBY1)
    d = (pointAX2 - pointAX1) * (pointBY2 - pointBY1) - (pointAY2 - pointAY1) * (pointBX2 - pointBX1)

    if d == 0:
        return False

    r = q / d
    
    q = (pointAY1 - pointBY1) * (pointAX2 - pointAX1) - (pointAX1 - pointBX1) * (pointAY2 - pointAY1)
    s = q / d

    return not (r < 0 or r > 1 or s < 0 or s > 1)

def pointInRectangle(pointX, pointY, uCornerX, uCornerY, lCornerX, lCornerY):
    """
    Returns true if the given point is inside the given rectangle.
    """
    return pointX >= uCornerX and pointX <= lCornerX and pointY >= uCornerY and pointY <= lCornerY

class QuadSet:
    """
    A quasi-quad tree implementation to accelerate the searching of points-on-lines.
    """
    def __init__(self, resolution, uCornerX, uCornerY, lCornerX, lCornerY):
        """
        Sets up a QuadSet such that the size of the resulting rectangle will be finer than the resolution value given.
        @type resolution: float
        @type uCornerX: float
        @type uCornerY: float
        @type lCornerX: float
        @type lCornerY: float
        """
        self.layers = int(math.ceil(max(math.log(abs(lCornerX - uCornerX) / resolution, 2), math.log(abs(lCornerY - uCornerY) / resolution, 2)))) + 1
        centerX = (uCornerX + lCornerX) / 2
        centerY = (uCornerY + lCornerY) / 2
        if self.layers > 2:
            self.quadElement = _QuadElement(self, 1, centerX - (resolution * 2 ** (self.layers - 1)) / 2, centerY - (resolution * 2 ** (self.layers - 1)) / 2,
                                        centerX + (resolution * 2 ** (self.layers - 1)) / 2, centerY + (resolution * 2 ** (self.layers - 1)) / 2)
        else:
            self.quadElement = _QuadBottom(self, 1, centerX - (resolution * 2 ** (self.layers - 1)) / 2, centerY - (resolution * 2 ** (self.layers - 1)) / 2,
                                        centerX + (resolution * 2 ** (self.layers - 1)) / 2, centerY + (resolution * 2 ** (self.layers - 1)) / 2)

    def storeLine(self, pointX1, pointY1, pointX2, pointY2, obj):
        """
        Drills down and stores the given line in all _QuadElements that intersect it.
        @type pointX1: float
        @type pointY1: float
        @type pointX2: float
        @type pointY2: float
        @type obj: ?
        """
        self.quadElement.storeLine(_QuadElementLine(pointX1, pointY1, pointX2, pointY2, obj))

    def retrieveLines(self, pointX, pointY, maxRadius=None):
        """
        Sets up a generator that returns all line objects in order of perpendicular distance from the
        point. The yield is the tuple (refDistance, minDistance, lineDist, perpendicular, obj) 
        """
        maxRadiusSq = maxRadius ** 2 if maxRadius is not None else sys.float_info.max
        heap = []
        traversedSet = set()
        heappush(heap, (0.0, self.quadElement))
        
        while len(heap) > 0:
            distSq, element = heappop(heap)
            if isinstance(element, _QuadElement):
                # Process the next layer
                for quadElement in element.members:
                    if quadElement is not None:
                        distSq = quadElement.minimalDistanceSq(pointX, pointY)
                        if distSq <= maxRadiusSq:
                            heappush(heap, (distSq, quadElement))
            elif isinstance(element, _QuadBottom):
                # We got to the bottom; evaluate all lines that intersect the rectangle:
                for line in element.members:
                    distSq, lineDist, perpendicular = pointDistSq(pointX, pointY, line.pointX1, line.pointY1, line.pointX2,
                                                                  line.pointY2, line.norm)
                    if distSq <= maxRadiusSq and line not in traversedSet:
                        heappush(heap, (distSq, (lineDist, perpendicular, line.obj)))
                        traversedSet.add(line)
            else: # It's a line, annotated with the stuff that was computed earlier:
                lineDist, perpendicular, lineObj = element
                yield (math.sqrt(distSq), lineDist, perpendicular, lineObj)

class _QuadElementLine:
    """
    Storage object for a line to be used by _QuadElement.
    """
    def __init__(self, pointX1, pointY1, pointX2, pointY2, obj):
        self.pointX1 = pointX1 
        self.pointY1 = pointY1 
        self.pointX2 = pointX2 
        self.pointY2 = pointY2
        self.norm = getNorm(pointX1, pointY1, pointX2, pointY2)
        self.obj = obj
    
class _QuadBase:
    """
    Base class for _QuadElement and _QuadBottom.
    @type quadSet: QuadSet
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
        """
        self.quadSet = quadSet
        self.layer = layer
        self.uCornerX = uCornerX
        self.uCornerY = uCornerY
        self.lCornerX = lCornerX
        self.lCornerY = lCornerY

    """
    def centralRadiusSq(self, pointX, pointY):
        ""
        Returns the distance from the given point to the center of this rectangle. 
        ""
        return getNormSq(pointX, pointY, (self.lCornerX - self.uCornerX) / 2, (self.lCornerY - self.uCornerY) / 2)
    """
    
    def minimalDistanceSq(self, pointX, pointY, pad=0.01):
        """
        Returns the squared distance from the given point to the nearest side or corner of this rectangle, adjusted
        for padding. If the point is in the rectangle, then 0.0 is returned.
        """
        if pointInRectangle(pointX, pointY, self.uCornerX, self.uCornerY, self.lCornerX, self.lCornerY):
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
    
    """
    def minimalRadiusSq(self, pointX, pointY, pad=1.0):
        ""
        Returns the squared distance from the given point to the nearest corner of this rectangle, adjusted for padding.
        If the point is in the rectangle, then 0.0 is returned.
        ""
        if pointX < self.lCornerX and pointX >= self.uCornerX and pointY < self.lCornerY and pointY >= self.uCornerX:
            # We are inside the rectangle. Return 0.0.
            ourRadius = 0.0
        else:
            if abs(self.uCornerX - pointX) < abs(self.lCornerX - pointX):
                if abs(self.uCornerY - pointY) < abs(self.lCornerY - pointY):
                    ourRadius = getNormSq(pointX, pointY, self.uCornerX + pad, self.uCornerY + pad)
                else:
                    ourRadius = getNormSq(pointX, pointY, self.uCornerX + pad, self.lCornerY - pad)
            else:
                if abs(self.uCornerY - pointY) < abs(self.lCornerY - pointY):
                    ourRadius = getNormSq(pointX, pointY, self.lCornerX - pad, self.uCornerY + pad)
                else:
                    ourRadius = getNormSq(pointX, pointY, self.lCornerX - pad, self.lCornerY - pad)
        return ourRadius
    """

    """
    def intersectsRadiusSq(self, pointX, pointY, radiusSq):
        ""
        Returns true if this _QuadElement intersects the circle depicted by the given radius squared. 
        ""
        return self.minimalRadiusSq(pointX, pointY) <= radiusSq
    """

    """
    def intersectsLine(self, pointX1, pointY1, pointX2, pointY2):
        ""
        Returns true if this _QuadElement intersects the given line.
        ""
        return lineIntersectsRectangle(pointX1, pointY1, pointX2, pointY2, self.uCornerX, self.uCornerY, self.lCornerX, self.lCornerY)
    """
        
class _QuadElement(_QuadBase):
    """
    Represents a rectangle within this QuadSet.
    @ivar members: set<_QuadElement>
    """
    def __init__(self, quadSet, layer, uCornerX, uCornerY, lCornerX, lCornerY):
        """
        @type quadSet: QuadSet
        @type layer: int
        @type uCornerX: float
        @type uCornerY: float
        @type lCornerX: float
        @type lCornerY: float
        """
        _QuadBase.__init__(self, quadSet, layer, uCornerX, uCornerY, lCornerX, lCornerY)
        self.members = 4 * [None]
            
    def _prepareBelow(self, index, uCornerX, uCornerY, lCornerX, lCornerY):
        """
        Checks to see if the quad below exists, and if it doesn't, creates it, which may be a _QuadElement
        or a _QuadBottom.
        """
        if self.members[index] is None:
            if self.layer + 1 < self.quadSet.layers:
                self.members[index] = _QuadElement(self.quadSet, self.layer + 1, uCornerX, uCornerY, lCornerX, lCornerY)
            else:
                self.members[index] = _QuadBottom(self.quadSet, self.layer + 1, uCornerX, uCornerY, lCornerX, lCornerY)
        return self.members[index]
    
    def storeLine(self, line):
        """
        Drills down and stores the given line in all _QuadElements that intersect it.
        @type line _QuadElementLine
        """
        centerX = (self.uCornerX + self.lCornerX) / 2
        centerY = (self.uCornerY + self.lCornerY) / 2
        for index in range(0, 4):
            quadElement = None
            if index == 0:
                if lineIntersectsRectangle(line.pointX1, line.pointY1, line.pointX2, line.pointY2, self.uCornerX, self.uCornerY, centerX, centerY):
                    quadElement = self._prepareBelow(0, self.uCornerX, self.uCornerY, centerX, centerY)
            elif index == 1:
                if lineIntersectsRectangle(line.pointX1, line.pointY1, line.pointX2, line.pointY2, centerX, self.uCornerY, self.lCornerX, centerY):
                    quadElement = self._prepareBelow(1, centerX, self.uCornerY, self.lCornerX, centerY)
            elif index == 2:
                if lineIntersectsRectangle(line.pointX1, line.pointY1, line.pointX2, line.pointY2, centerX, centerY, self.lCornerX, self.lCornerY):
                    quadElement = self._prepareBelow(2, centerX, centerY, self.lCornerX, self.lCornerY)
            elif index == 3:
                if lineIntersectsRectangle(line.pointX1, line.pointY1, line.pointX2, line.pointY2, self.uCornerX, centerY, centerX, self.lCornerY):
                    quadElement = self._prepareBelow(3, self.uCornerX, centerY, centerX, self.lCornerY)
            
            if quadElement is not None:
                quadElement.storeLine(line)

    """
    def retrieveLines(self, pointX, pointY, maxRadiusSq, traversedSet=set()):
        ""
        Continues the generation of perpendicular points. The yield is the tuple (refDistance, minDistance, perpendicular, obj) 
        ""
        for quadElement in self.members:
            heappush(self.quadSet.heap, (quadElement.minimalRadiusSq(pointX, pointY), quadElement))
        
        # Process the lowest-scoring item in the heap.
        # *** WHILE?
        _, bestQuadElement = heappop(self.quadSet.heap) 
        for refDistance, minDistance, lineDist, perpendicular, obj in bestQuadElement.retrieveLines(pointX, pointY, maxRadiusSq, traversedSet):
            yield (refDistance, minDistance, lineDist, perpendicular, obj)
        
        ""
        radii = []
        for quadElement in self.members:
            if quadElement is not None:
                radii.append((quadElement.minimalRadiusSq(pointX, pointY), quadElement))
        radii[:] = sorted(radii)
        for _, quadElement in radii:
            for refDistance, minDistance, lineDist, perpendicular, obj in quadElement.retrieveLines(pointX, pointY, maxRadiusSq, traversedSet):
                yield (refDistance, minDistance, lineDist, perpendicular, obj)
        ""
    """
            
class _QuadBottom(_QuadBase):
    """
    Represents a bottom-most rectangle within the set which has references to actual lines.
    """
    def __init__(self, quadSet, layer, uCornerX, uCornerY, lCornerX, lCornerY):
        """
        @type quadSet: QuadSet
        @type layer: int
        @type uCornerX: float
        @type uCornerY: float
        @type lCornerX: float
        @type lCornerY: float
        """
        _QuadBase.__init__(self, quadSet, layer, uCornerX, uCornerY, lCornerX, lCornerY)
        self.members = []
    
    def storeLine(self, line):
        """
        Stores the given line in this _QuadBottom.
        """
        self.members.append(line)

    """
    def retrieveLines(self, pointX, pointY, maxRadiusSq, traversedSet=set()):
        ""
        Continues the generation of perpendicular points. The yield is the tuple (refDistance, minDistance, perpendicular, obj) 
        ""
        refDists = []
        for line in self.members:
            distSq, lineDist, perpendicular = pointDistSq(pointX, pointY, line.pointX1, line.pointY1, line.pointX2, line.pointY2, line.norm)
            if distSq <= maxRadiusSq and line not in traversedSet:
                refDists.append((distSq, lineDist, perpendicular, line))
                traversedSet.add(line)
        if len(refDists) > 0:
            refDists[:] = sorted(refDists)
            minDistance = math.sqrt(self.minimalRadiusSq(pointX, pointY))
            for distSq, lineDist, perpendicular, line in refDists:
                yield (math.sqrt(distSq), minDistance, lineDist, perpendicular, line.obj)
    """

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
    