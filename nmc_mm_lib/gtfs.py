"""
gtfs.py: Definitions for entities that exist within GTFS datasets
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
from collections.abc import Iterable
import csv
import logging
import os, operator
from datetime import datetime, timedelta
from typing import Hashable, NamedTuple, Self

# TODO: This would be more compact with Pandas.

class ShapesEntry(NamedTuple):
    """
    ShapesEntry is a single GTFS shape file entry.
    """
    shapeID: int
    shapeSeq: int
    lat: float
    lng: float

    def __hash__(self) -> int:
        return hash((self.shapeID, self.shapeSeq))
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ShapesEntry):
            return NotImplemented
        return self.shapeID == other.shapeID and self.shapeSeq == other.shapeSeq

def fillShapes(filePath: str) -> dict[int, list[ShapesEntry]]:
    """
    fillShapes retrieves the shape information from a shape file and returns a list of shape entries.

    @return A map of shape_id to a list of shape entries
    """
    ret: dict[int, list[ShapesEntry]] = {}
    filename = os.path.join(filePath, "shapes.txt") 
    with open(filename, mode='r', newline='') as inFile:
        csvReader = csv.DictReader(inFile)
        for fileLine in csvReader:
            newEntry = ShapesEntry(shapeID=int(fileLine['shape_id']),
                                   shapeSeq=int(fileLine['shape_pt_sequence']),
                                   lat=float(fileLine['shape_pt_lat']),
                                   lng=float(fileLine['shape_pt_lon']))
            if newEntry.shapeID not in ret:
                ret[newEntry.shapeID] = []
            ret[newEntry.shapeID].append(newEntry)

    # Ensure that the lists are sorted:
    shapesEntries: list[ShapesEntry]
    for shapesEntries in ret.values():
        shapesEntries.sort(key=operator.attrgetter('shapeSeq'))

    # Return the shapes file contents:
    return ret

class RoutesEntry(NamedTuple):
    """
    RoutesEntry is a single GTFS route with name.
    """
    routeID: int
    shortName: str
    name: str

    def __hash__(self) -> int:
        return self.routeID
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, RoutesEntry):
            return NotImplemented
        return self.routeID == other.routeID
        
def fillRoutes(filePath: str) -> dict[int, RoutesEntry]:
    """
    fillRoutes retrieves route name information from a GTFS repository.

    @return A map from routeID to a RoutesEntry object.
    """
    ret: dict[int, RoutesEntry] = {}
    filename = os.path.join(filePath, "routes.txt") 
    with open(filename, mode='r', newline='') as inFile:
        csvReader = csv.DictReader(inFile)
        
        # Go through the lines of the file:
        for fileLine in csvReader:
            newEntry = RoutesEntry(routeID=int(fileLine['route_id']),
                                   shortName=fileLine['route_short_name'],
                                   name=fileLine['route_long_name'])
            ret[newEntry.routeID] = newEntry
                    
    # Return the routes file contents:
    return ret    

class TripsEntry(NamedTuple):
    """
    TripsEntry is a single GTFS trip file entry.  Its key is tripID.
    """
    tripID: int
    route: RoutesEntry
    tripHeadsign: str
    shapeEntries: list[ShapesEntry]
        
    def __hash__(self) -> int:
        return self.tripID
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, TripsEntry):
            return NotImplemented
        return self.tripID == other.tripID

def fillTrips(filePath: str,
              shapes: dict[int, list[ShapesEntry]],
              routes: dict[int, RoutesEntry],
              unusedShapeIDs: Iterable[int] | set[int] = {},
              restrictService: Iterable[str] | set[str] = {}) -> tuple[dict[int, TripsEntry], set[int]]:
    """
    fillTrips retrieves the trip information from a GTFS repository.

    @return A map of trip_id to TripsEntry records, as well as a list of unused trip IDs
    """
    ret: dict[int, TripsEntry] = {}
    unusedTripIDs: set[int] = set()
    shapeErrorIDs: set[int] = set()
    unusedShapeIDs = set(unusedShapeIDs)
    restrictService = set(restrictService)

    filename = os.path.join(filePath, "trips.txt") 
    with open(filename, mode='r', newline='') as inFile:
        csvReader = csv.DictReader(inFile)

        # Go through the lines of the file:
        for fileLine in csvReader:
            shapeID = int(fileLine['shape_id'])
            routeID = int(fileLine['route_id'])
            serviceID = fileLine['service_id']
            tripID = int(fileLine['trip_id'])
            if shapeID in unusedShapeIDs or len(restrictService) > 0 and serviceID not in restrictService:
                unusedTripIDs.add(tripID)
            else:
                if shapeID not in shapes:
                    if shapeID not in shapeErrorIDs:
                        logging.warning(f"GTFS Trip {tripID} expects undefined shape ID {shapeID}; skipping")
                        shapeErrorIDs.add(shapeID)                        
                    unusedTripIDs.add(tripID)
                else:
                    if routeID not in routes:
                        logging.warning(f"GTFS Trip {tripID} expects undefined route ID {routeID}; skipping")
                        unusedTripIDs.add(tripID)                        
                    else:
                        newEntry = TripsEntry(tripID=tripID,
                                              route=routes[routeID],
                                              tripHeadsign=fileLine['trip_headsign'],
                                              shapeEntries=shapes[shapeID])
                        ret[newEntry.tripID] = newEntry
                            
    # Return the trips file contents:
    return ret, unusedTripIDs

class StopsEntry(NamedTuple):
    """
    StopsEntry is a single GTFS stops file entry.
    """
    stopID: int
    stopName: str
    gpsLat: float
    gpsLng: float

    def __hash__(self) -> int:
        return self.stopID
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, StopsEntry):
            return NotImplemented
        return self.stopID == other.stopID

def fillStops(filePath: str) -> dict[int, StopsEntry]:
    """
    fillStops retrieves the stop information from a GTFS repository.

    @return A map of stop_id to a StopsEntry.
    """
    ret: dict[int, StopsEntry] = {}
    filename = os.path.join(filePath, "stops.txt") 
    with open(filename, mode='r', newline='') as inFile:
        csvReader = csv.DictReader(inFile)

        # Go through the lines of the file:
        for fileLine in csvReader:
            newEntry = StopsEntry(stopID=int(fileLine['stop_id']),
                                  stopName=fileLine['stop_name'],
                                  gpsLat=float(fileLine['stop_lat']),
                                  gpsLng=float(fileLine['stop_lon']))
            ret[newEntry.stopID] = newEntry
    
    # Return the stops file contents:
    return ret

class StopTimesEntry(NamedTuple):
    """
    StopTimesEntry is a single GTFS stoptimes file entry.
    """
    trip: TripsEntry
    stop: StopsEntry
    stopSeq: int
    arrivalTime: datetime
    departureTime: datetime

    def __hash__(self) -> int:
        return hash((self.trip, self.stopSeq))
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, StopTimesEntry):
            return NotImplemented
        return self.trip == other.trip and self.stopSeq == other.stopSeq

def parseGTFSTime(timeStr: str) -> datetime:
    """
    parseGTFSTime parses a GTFS time string into a datetime object.

    @return A datetime object representing the given time string, with respect to 1/1/1900
    """
    # Split apart time string this way and count from epoch because GTFS stops may express times for
    # early morning service with hours being greater than 23.
    timeElems = timeStr.split(':')
    timeHour = int(timeElems[0])
    timeDays = timeHour // 24
    timeHour = timeHour % 24
    timeObj = datetime(1900, 1, 1, timeHour, int(timeElems[1]), int(timeElems[2]))
    timeObj += timedelta(days = timeDays)
    return timeObj

def fillStopTimes(filePath: str,
                  trips: dict[int, TripsEntry],
                  stops: dict[int, StopsEntry],
                  unusedTripIDs: Iterable[int]) -> dict[TripsEntry, list[StopTimesEntry]]:
    """
    fillStopTimes retrieves the stoptime information from a GTFS repository.

    @return A map of TripsEntry to a list of stop entries plus the start and end times
    """
    stopTimes: dict[TripsEntry, list[StopTimesEntry]] = {}
    unusedTripIDs = set(unusedTripIDs)
    
    filename = os.path.join(filePath, "stop_times.txt")
    with open(filename, mode='r', newline='') as inFile:
        csvReader = csv.DictReader(inFile)
        
        # Go through the lines of the file:
        badTrips: set[int] = set()
        for fileLine in csvReader:
            tripID = int(fileLine['trip_id'])
            if tripID not in unusedTripIDs:
                if not tripID in trips:
                    badTrips.add(tripID)
                    continue                
                arrivalTime = parseGTFSTime(fileLine['arrival_time'])
                departureTime = parseGTFSTime(fileLine['departure_time'])
                stopID = int(fileLine['stop_id'])
                if not stopID in stops:
                    logging.warning(f"GTFS Stop Times file expects undefined stop ID {stopID}")
                    continue
                newEntry = StopTimesEntry(trip=trips[tripID],
                                          stop=stops[stopID],
                                          stopSeq=int(fileLine['stop_sequence']),
                                          arrivalTime=arrivalTime,
                                          departureTime=departureTime)
                if newEntry.trip not in stopTimes:
                    stopTimes[newEntry.trip] = []
                stopTimes[newEntry.trip].append(newEntry)

        # Output error message:
        if badTrips:
            strOut = ""
            for tripID in sorted(badTrips):
                if strOut:
                    strOut += ", "
                strOut += str(tripID)
            logging.warning(f"GTFS Stop Times file expects undefined trip ID(s) {strOut}")

    # Sort the stop times by stop sequence:
    stopTimesList: list[StopTimesEntry]
    for stopTimesList in stopTimes.values():
        stopTimesList.sort(key=operator.attrgetter('stopSeq'))

    # Return the stop times file contents:
    return stopTimes

class GTFSSet:
    """
    Represents the entire contents of a GTFS set of files
    """
    shapes: dict[int, list[ShapesEntry]]
    routes: dict[int, RoutesEntry]
    trips: dict[int, TripsEntry]
    stops: dict[int, StopsEntry]
    stopTimes: dict[TripsEntry, list[StopTimesEntry]]

    def __init__(self, filepath: str):
        """
        Loads in the contents of a GTFS set
        
        @param filepath: Diretory in which GTFS set sits
        """
        self.shapes = fillShapes(filepath)
        self.routes = fillRoutes(filepath)
        self.trips, unusedTripIDs = fillTrips(filepath, self.shapes, self.routes)
        self.stops = fillStops(filepath)
        self.stopTimes = fillStopTimes(filepath, self.trips, self.stops, unusedTripIDs)
