"""
track.py: Definitions for entities that exist within shape datasets
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

class Tracks(dict):
    """
    Tracks is a dictionary of Track objects, organized according to trackID -> Track
    """
    pass

class Track(list):
    """
    Track is a single shape that contains multiple TrackEntry objects arranged by sequence number.
    @ivar trackID: an undefined type for this shape's unique identifier
    """
    def __init__(self, trackID):
        """
        @type trackID: str
        """
        self.trackID = trackID
    
class Trackpoint:
    """
    Trackpoint is a single shape dataset entry.
    @ivar trackID: Identifier that signifies the parent of this TrackEntry
    @ivar trackSeq: int signifying the ordering of this element relative to others.
        The pair trackID, trackSeq are to be unique.
    @ivar lat: float signifying latitude
    @ivar lng: float signifying longitude
    @ivar hintFlag: bool signifying whether this TrackEntry is a hint or not.
    @ivar pointX: float signifying the translated rectangular coordinate X value.
    @ivar pointY: float signifying the translated rectangular coordinate Y value.
    @ivar name: str that is a name attached to this Trackpoint
    """
    def __init__(self, trackID, trackSeq, lat, lng, hintFlag=False):
        """
        @type trackID: str
        @type trackSeq: int
        @type lat: float
        @type lng: float
        @type hintFlag: bool
        """
        self.trackID = trackID
        self.trackSeq = trackSeq
        self.lat = lat
        self.lng = lng
        self.hintFlag = hintFlag
        self.name = str + "-" + str(trackSeq)
        
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
    @ivar routeID: str signifying a route ID
    @ivar shortName: str that is the short name of the route
    @ivar name: str that is the full name of the route
    """
    def __init__(self, routeID, shortName, name):
        """
        @type routeID: str
        @type shortName: str
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
    @ival tripID: str signifying a trip ID
    @ival route: RoutesEntry that the trip is a member of
    @ival tripHeadsign: str that signifies the direction that the trip is taking
    @ival track: Shape that captures all of the track elements that comprise this trip, or None if not defined
    """
    def __init__(self, tripID, route, tripHeadsign, track=None):
        """
        @type tripID: str
        @type route: RoutesEntry
        @type tripHeadsign: str
        @type track: Track
        """
        self.tripID = tripID
        self.route = route
        self.tripHeadsign = tripHeadsign
        self.track = track

class Timepoints(dict):
    """
    Timepoints represents Trackpoints that have timestamps associated with them, arranged as a 
    dictionary of tripID -> list<StopTime>.
    """
    pass

class Timepoint:
    """
    Timepoint ties a time to a Trackpoint.
    @ivar trip: Trip that this StopTime is a member of
    @ivar trackpoint: The reference to the individual trackpoint
    @ivar timeSeq: int that is a sequence number to arrange Timepoints
    @ivar arrivalTime: datetime that is the arrival time of this StopTime
    @ivar departureTime: datetime that is the departure time of this StopTime
    """
    def __init__(self, trip, trackpoint, timeSeq):
        """
        @type trip: Trip
        @type trackpoint: Trackpoint
        @type timeSeq: int
        """
        self.trip = trip
        self.trackpoint = trackpoint
        self.timeSeq = timeSeq
        
        self.arrivalTime = None
        self.departureTime = None
