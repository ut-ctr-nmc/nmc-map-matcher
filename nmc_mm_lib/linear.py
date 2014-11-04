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
import math, unittest

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
    