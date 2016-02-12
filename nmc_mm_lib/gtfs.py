"""
gtfs.py: Definitions for entities that exist within GTFS datasets
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
import os, operator, sys
from datetime import datetime, timedelta

class ShapesEntry:
    """
    ShapesEntry is a single GTFS shape file entry.
    """
    def __init__(self, shapeID, shapeSeq, lat, lng, hintFlag = False):
        """
        @type shapeSeq: int
        @type lat: float
        @type lng: float
        @type hintFlag: bool
        @ivar typeID: int
        """
        self.shapeID = shapeID
        self.shapeSeq = shapeSeq
        self.lat = lat
        self.lng = lng
        self.hintFlag = hintFlag
        
        self.pointX = 0
        self.pointY = 0
        
        self.typeID = 0

def fillShapes(filePath, gps):
    """
    fillShapes retrieves the shape information from a shape file and returns a list of shape entries.
    @type filePath: str
    @type gps: gps.GPS
    @return A map of shape_id to a list of shape entries
    @rtype dict<int, list<ShapesEntry>>
    """
    ret = {}
    "@type ret: dict<int, list<ShapesEntry>>"
    filename = os.path.join(filePath, "shapes.txt") 
    with open(filename, 'r') as inFile:
        # Sanity check:
        fileLine = inFile.readline()
        if not fileLine.startswith("shape_id,shape_pt_lat,shape_pt_lon"):
            print("ERROR: The shapes.txt file doesn't have the expected header.", file = sys.stderr)
            return None
        
        # Go through the lines of the file:
        for fileLine in inFile:
            if len(fileLine) > 0:
                lineElems = fileLine.split(',')
                newEntry = ShapesEntry(int(lineElems[0]), int(lineElems[3]), float(lineElems[1]),
                                        float(lineElems[2]), False)
                (newEntry.pointX, newEntry.pointY) = gps.gps2feet(newEntry.lat, newEntry.lng)
                if newEntry.shapeID not in ret:
                    ret[newEntry.shapeID] = []
                ret[newEntry.shapeID].append(newEntry)

    # Ensure that the lists are sorted:
    for shapesEntries in ret.values():
        "@type shapesEntries: list<ShapesEntry>"
        shapesEntries.sort(key = operator.attrgetter('shapeSeq'))
                    
    # Return the shapes file contents:
    return ret

class RoutesEntry:
    """
    RoutesEntry is a single GTFS route with name.
    """
    def __init__(self, routeID, shortName, name):
        """
        @type routeID: int
        @type shortName: str
        @type name: str
        """
        self.routeID = routeID
        self.shortName = shortName
        self.name = name
        
def fillRoutes(filePath):
    """
    fillRoutes retrieves route name information from a GTFS repository.
    @type filepath: str
    @return A map from routeID to a RoutesEntry object.
    @rtype dict<int, RoutesEntry>
    """
    ret = {}
    "@type ret: dict<int, RoutesEntry>"
    filename = os.path.join(filePath, "routes.txt") 
    with open(filename, 'r') as inFile:
        # Sanity check:
        fileLine = inFile.readline()
        if not fileLine.startswith("route_id,agency_id,route_short_name,route_long_name"):
            print("ERROR: The routes.txt file doesn't have the expected header.", file = sys.stderr)
            return None
        
        # Go through the lines of the file:
        for fileLine in inFile:
            if len(fileLine) > 0:
                lineElems = fileLine.split(',')
                newEntry = RoutesEntry(int(lineElems[0]), lineElems[2], lineElems[3])
                ret[newEntry.routeID] = newEntry
                    
    # Return the stops file contents:
    return ret    

class TripsEntry:
    """
    TripsEntry is a single GTFS trip file entry.  Its key is tripID.
    """
    def __init__(self, tripID, route, tripHeadsign, shapeEntries):
        """
        @type tripID: int
        @type route: RoutesEntry
        @type tripHeadsign: str
        @type shapeEntries: list<ShapesEntry>
        """
        self.tripID = tripID
        self.route = route
        self.tripHeadsign = tripHeadsign
        self.shapeEntries = shapeEntries
        
    def __hash__(self):
        return self.tripID
    
    def __eq__(self, other):
        return self.tripID == other.tripID

def fillTrips(filePath, shapes, routes, unusedShapeIDs=None, restrictService=None):
    """
    fillTrips retrieves the trip information from a GTFS repository.
    @type filePath: str
    @type shapes: dict<int, list<ShapesEntry>>
    @type routes: dict<int, RoutesEntry>
    @type unusedShapeIDs: set<int>
    @type restrictService: set<string>
    @return A map of trip_id to TripsEntry records, as well as a list of unused trip IDs
    @rtype (dict<int, TripsEntry>, set<int>)
    """
    if not unusedShapeIDs:
        unusedShapeIDs = set()
    if not restrictService:
        restrictService = set()
    ret = {}
    "@type ret: dict<int, TripsEntry>"
    unusedTripIDs = set()
    "@type unusedTripIDs: set<int>"
    filename = os.path.join(filePath, "trips.txt") 
    with open(filename, 'r') as inFile:
        # Sanity check:
        fileLine = inFile.readline()
        if not fileLine.startswith("route_id,service_id,trip_id,trip_headsign,trip_short_name,direction_id,block_id,shape_id"):
            print("ERROR: The trips.txt file doesn't have the expected header.", file = sys.stderr)
            return None
        
        # Go through the lines of the file:
        for fileLine in inFile:
            if len(fileLine) > 0:
                lineElems = fileLine.split(',')
                shapeID = int(lineElems[7])
                routeID = int(lineElems[0])
                serviceID = lineElems[1]
                tripID = int(lineElems[2])
                if (shapeID in unusedShapeIDs) or ((len(restrictService) > 0) and (serviceID not in restrictService)):
                    unusedTripIDs.add(tripID)
                else:
                    if shapeID not in shapes:
                        print("WARNING: GTFS Trips file expects undefined shape ID %d" % shapeID, file = sys.stderr)
                        unusedTripIDs.add(tripID)                        
                    else:
                        if routeID not in routes:
                            print("WARNING: GTFS Trips file expects undefined route ID %d" % routeID, file = sys.stderr)
                            unusedTripIDs.add(tripID)                        
                        else:
                            newEntry = TripsEntry(tripID, routes[routeID], lineElems[3], shapes[shapeID])
                            ret[newEntry.tripID] = newEntry
                            
    # Return the trips file contents:
    return (ret, unusedTripIDs)

class StopsEntry:
    """
    StopsEntry is a single GTFS stops file entry.
    """
    def __init__(self, stopID, stopName, gpsLat, gpsLng):
        """
        @type stopID: int
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

def fillStops(filePath, gps):
    """
    fillStops retrieves the stop information from a GTFS repository.
    @type filePath: str
    @type GPS: gps.GPS
    @return A map of stop_id to a StopsEntry.
    @rtype dict<int, StopsEntry>
    """
    ret = {}
    "@type ret: dict<int, StopsEntry>"
    filename = os.path.join(filePath, "stops.txt") 
    with open(filename, 'r') as inFile:
        # Sanity check:
        fileLine = inFile.readline()
        if not fileLine.startswith("stop_id,stop_code,stop_name,stop_desc,stop_lat,stop_lon"):
            print("ERROR: The stops.txt file doesn't have the expected header.", file = sys.stderr)
            return None
        
        # Go through the lines of the file:
        for fileLine in inFile:
            if len(fileLine) > 0:
                lineElems = fileLine.split(',')
                newEntry = StopsEntry(int(lineElems[0]), lineElems[2], float(lineElems[4]), float(lineElems[5]))
                (newEntry.pointX, newEntry.pointY) = gps.gps2feet(newEntry.gpsLat, newEntry.gpsLng)
                ret[newEntry.stopID] = newEntry
                    
    # Return the stops file contents:
    return ret

class StopTimesEntry:
    """
    StopTimesEntry is a single GTFS stoptimes file entry.
    """
    def __init__(self, trip, stop, stopSeq):
        """
        @type trip: TripsEntry
        @type stop: StopsEntry
        @type stopSeq: int
        @type arrivalTime: datetime
        @type departureTime: datetime
        """
        self.trip = trip
        self.stop = stop
        self.stopSeq = stopSeq
        
        self.arrivalTime = None
        self.departureTime = None

def fillStopTimes(filePath, trips, stops, unusedTripIDs):
    """
    fillStopTimes retrieves the stoptime information from a GTFS repository.
    @type filePath: str
    @type trips: dict<int, TripsEntry>
    @type stops: dict<int, StopsEntry>
    @type unusedTripIDs: set<int>
    @return A map of TripsEntry to a list of stop entries plus the start and end times
    @rtype dict<TripsEntry, list<StopTimesEntry>>
    """
    stopTimes = {}
    "@type stopTimes: dict<TripsEntry, list<StopTimesEntry>>"
    
    filename = os.path.join(filePath, "stop_times.txt")
    with open(filename, 'r') as inFile:
        # Sanity check:
        fileLine = inFile.readline()
        if not fileLine.startswith("trip_id,arrival_time,departure_time,stop_id"):
            print("ERROR: The stop_times.txt file doesn't have the expected header.", file = sys.stderr)
            return None
        
        # Go through the lines of the file:
        for fileLine in inFile:
            if len(fileLine) > 0:
                lineElems = fileLine.split(',')
                tripID = int(lineElems[0])
                if tripID not in unusedTripIDs:
                    if not tripID in trips:
                        print("WARNING: GTFS Stop Times file expects undefined trip ID %d" % tripID, file = sys.stderr)
                        continue
                    
                    # Split apart time string because hours can roll around.
                    timeElems = lineElems[1].split(':')
                    timeHour = int(timeElems[0])
                    timeDays = timeHour / 24
                    timeHour = timeHour % 24
                    arrivalTime = datetime(1900, 1, 1, timeHour, int(timeElems[1]), int(timeElems[2]))
                    arrivalTime += timedelta(days = timeDays)

                    timeElems = lineElems[2].split(':')
                    timeHour = int(timeElems[0])
                    timeDays = timeHour / 24
                    timeHour = timeHour % 24
                    departureTime = datetime(1900, 1, 1, timeHour, int(timeElems[1]), int(timeElems[2]))
                    departureTime += timedelta(days = timeDays)
                    
                    stopID = int(lineElems[3])
                    if not stopID in stops:
                        print("WARNING: GTFS Stop Times file expects undefined stop ID %d" % stopID, file = sys.stderr)
                        continue
                    newEntry = StopTimesEntry(trips[tripID], stops[stopID], int(lineElems[4]))
                    newEntry.arrivalTime = arrivalTime
                    newEntry.departureTime = departureTime
                    if newEntry.trip not in stopTimes:
                        stopTimes[newEntry.trip] = []
                    stopTimes[newEntry.trip].append(newEntry)
                    del newEntry
                del tripID

    # Return the stop times file contents:
    return stopTimes
