"""
shape.py: Definitions for entities that exist within shape datasets
@author: Kenneth Perrine
@contact: kperrine@utexas.edu
@organization: Network Modeling Center, Center for Transportation Research,
    Cockrell School of Engineering, The University of Texas at Austin 
@version: 1.1

@copyright: (C) 2015, The University of Texas at Austin
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

class Shapes(dict):
    """
    Shapes is a dictionary of Shape objects, organized according to shapeID -> Shape
    """
    pass

class Shape(list):
    """
    Shape is a single shape that contains multiple ShapesEntry objects arranged by sequence number.
    @ivar shapeID: an undefined type for this shape's unique identifier
    """
    def __init__(self, shapeID):
        """
        @type shapeID: ?
        """
        self.shapeID = shapeID
    
class ShapesEntry:
    """
    ShapesEntry is a single shape dataset entry.
    @ivar shapeID: ? signifying the shape ID of this ShapesEntry
    @ivar shapeSeq: int signifying the ordering of this shape relative to others
    @ivar lat: float signifying latitude
    @ivar lng: float signifying longitude
    @ivar time: datetime object representing a timestamp, or None if undefined.
    @ivar hintFlag: bool signifying whether this ShapesEntry is a hint or not.
    @ivar pointX: float signifying the translated rectangular coordinate X value.
    @ivar pointY: float signifying the translated rectangular coordinate Y value.
    """
    def __init__(self, shape, shapeSeq, lat, lng, hintFlag=False):
        """
        @type shape: Shape
        @type shapeSeq: int
        @type lat: float
        @type lng: float
        @type hintFlag: bool
        """
        self.shape = shape
        self.shapeSeq = shapeSeq
        self.lat = lat
        self.lng = lng
        self.hintFlag = hintFlag
        self.time = None
        
        self.pointX = 0
        self.pointY = 0

class Routes(dict):
    """
    Routes is a dictionary of routeID -> Route.
    """
    pass

class Route:
    """
    Route is a single route with name.
    @ivar routeID: ? signifying a route ID
    @ivar shortName: ? that is the short name of the route
    @ivar name: str that is the full name of the route
    """
    def __init__(self, routeID, shortName, name):
        """
        @type routeID: ?
        @type shortName: ?
        @type name: str
        """
        self.routeID = routeID
        self.shortName = shortName
        self.name = name
     
class Trips(dict):
    """
    Trips is a dictionary of tripID -> Trip.
    """
    pass
        
class Trip:
    """
    Trip is a single trip entry.
    @ival tripID: ? signifying a trip ID
    @ival route: RoutesEntry that the trip is a member of
    @ival tripHeadsign: str that signifies the direction that the trip is taking
    @ival shape: Shape that captures all of the shape elements that comprise this trip, or None if not defined
    """
    def __init__(self, tripID, route, tripHeadsign, shape=None):
        """
        @type tripID: ?
        @type route: RoutesEntry
        @type tripHeadsign: str
        @type shape: Shape
        """
        self.tripID = tripID
        self.route = route
        self.tripHeadsign = tripHeadsign
        self.shape = shape

class Stops(dict):
    """
    Stops is a dictionary of stopID -> Stop
    """
    pass

class Stop:
    """
    Stop is a single stops file entry.
    @ivar stopID: ? representing the stop ID
    @ivar stopName: str representing the name of the stop
    @ivar gpsLat: float representing the latitude of the stop 
    @ivar gpsLng: float representing the longitude of the stop
    @ivar pointX: float representing the stop according to rectangular X coordinates around a reference point
    @ivar pointY: float representing the stop according to rectangular Y coordinates around a reference point
    """
    def __init__(self, stopID, stopName, gpsLat, gpsLng):
        """
        @type stopID: ?
        @type stopName: str
        @type gpsLat: float
        @type gpsLng: float
        """
        self.stopID = stopID
        self.stopName = stopName
        self.gpsLat = gpsLat
        self.gpsLng = gpsLng
        
        self.pointX = 0
        self.pointY = 0

class StopTimes(dict):
    """
    StopTimes represents a set of StopTime objects arranged as a dictionary of tripID -> list<StopTime>.
    """
    pass

class StopTime:
    """
    StopTime is a single stoptimes file entry.
    @ivar trip: Trip that this StopTime is a member of
    @ivar stop: Stop that this StopTime is a member of
    @ivar stopSeq: int that is the stop sequence number
    @ivar arrivalTime: datetime that is the arrival time of this StopTime
    @ivar departureTime: datetime that is the departure time of this StopTime
    """
    def __init__(self, trip, stop, stopSeq):
        """
        @type trip: TripsEntry
        @type stop: StopsEntry
        @type stopSeq: int
        """
        self.trip = trip
        self.stop = stop
        self.stopSeq = stopSeq
        
        self.arrivalTime = None
        self.departureTime = None
