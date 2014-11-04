"""
gps.py: GPS-related methods
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
import math

FT_PER_DEGREE = 69.172 * 5280

class GPS:
    """
    GPS class calculates distances from GPS coordinates using a center as a reference to allow for rectangular
    approximations in calculations. At this time, because of the rectangular approximation, this is not
    recommended for geographic areas that are larger than a metropolitan region, for example.
    
    TODO: Insert a better GPS coordinate approximator or a library.
    """
    def  __init__(self, latCtr, lngCtr):
        """
        In the constructor, provide the reference GPS center, so that the neighborhood can be treated as
        a rectangular grid.
        @type latCtr: float
        @type lngCtr: float
        """
        self.latCtr = latCtr
        self.latCtrCos = math.cos(latCtr)
        self.oneDegreeLngFt = self.latCtrCos * FT_PER_DEGREE
        self.lngCtr = lngCtr

    def gps2feet(self, lat, lng):
        """
        gps2feet converts the GPS coordinate pair to feet relative to the given center.
        @type lat: float
        @type lng: float
        @rtype float, float
        """
        # Length of 1 degree of Longitude = cosine (latitude) * length of degree (miles) at equator
        coordX = (lng - self.lngCtr) * self.oneDegreeLngFt
        
        # length of 1 degree of latitude = 69.172 miles
        coordY = (lat - self.latCtr) * FT_PER_DEGREE
        return (coordX, coordY)
        
    def feet2gps(self, coordX, coordY):
        """
        feet2gps converts the rectangular pair centered around the assumed center to GPS coordinates.
        @type coordX: float
        @type coordY: float
        @rtype float, float
        """
        lng = coordX / self.oneDegreeLngFt + self.lngCtr
        lat = coordY / FT_PER_DEGREE + self.latCtr
        return (lat, lng)
        
    def gps2dist(self, lat1, lng1, lat2, lng2):
        """
        gps2dist finds the distance between the GPS coordinate pair
        @type lat1: float
        @type lng1: float
        @type lat2: float
        @type lng2: float
        @rtype float
        """
        (coord1X, coord1Y) = self.gps2feet(lat1, lng1)
        (coord2X, coord2Y) = self.gps2feet(lat2, lng2)
        return math.sqrt((coord1Y - coord2Y) ** 2 + (coord1X - coord2X) ** 2)
        