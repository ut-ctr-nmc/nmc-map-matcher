"""
gtfs.py: Definitions for GTFS entities including shapes and routes
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
from __future__ import print_function
from nmc_mm_lib import shape, shape_data_provider
import os, operator, sys
from datetime import datetime, timedelta

class GTFS(shape_data_provider.ShapeDataProvider):
    """
    Represents a GTFS file set that represents shape data and associated transit data.
    
    @ivar filePath: str that is the path to the GTFS directory
    @ivar shapeFilename: str that holds the name of the shapefile; useful for something other than "shapes.txt"
    """
    def __init__(self):
        # These lists are here in case we get to the point that we need to maintain two different
        # shape sets and each needs its own list of exclusions or inclusions.
        self.filePath = None
        self.shapeFilename = None
        
        super(shape_data_provider.ShapeDataProvider, self).__init__()
    
    def shapeIDConstructor(self):
        """
        Returns a constructor that signifies the data type for a shape ID. For GTFS, this is an int.
        """
        return lambda x: int(x)
    
    def provideCmdLineOpts(self, options):
        """
        Stuffs new command-line options into the given argparse object. This performs the
        shapefile-specific command line specification setup. 
        @type options: argparse
        """
        raise NotImplementedError()        

    def checkCmdLineOpts(self, args):
        """
        Checks the given command-line options and determines if there are sufficient options
        for successfully loading in data for this class.
        @type args: argparse.Namespace
        @return True if sufficient command-line options are provided for this class
        @rtype boolean
        """
        raise NotImplementedError()
    
    def readCmdLineOpts(self, args):
        """
        Reads in all relevant command-line options and returns false if there are insufficient
        options provided.
        @type args: argparse.Namespace
        @return False if there are insufficient options provided, or True if successful.
        @rtype boolean
        """
        raise NotImplementedError()
    
    def readAllShapes(self):
        """
        Using information from the command-line options, reads in all shape data from the geographic
        data source that's handled by this class.
        @return A collection of shape ID mapped to Shape; that is, lists of ShapesEntry objects.
        @rtype shape.Shapes
        """
        ret = shape.Shapes()
        "@type ret: shape.Shape"
        
        for shapeEntry in readShapesEntries(suppressOutOfOrder=True):
            if shapeEntry.shapeID not in ret:
                ret[shapeEntry.shapeID] = shape.Shape(shapeEntry.shapeID)
            ret[shapeEntry.shapeID][shapeEntry.shapeSeq] = shapeEntry
            
        # In case there were out of order entries:
        for shapeID in ret:
            ret[shapeID].sort(key=lambda x: x.shapeSeq)
            
        # Return all of the shapes:
        return ret    
        
    def readShapesEntries(self, suppressOutOfOrder=False):
        """
        readShapesEntries returns a generator that gives out shape entries one by one.
        At each call, the next valid constructed shape.ShapesEntry is yielded. You must check the
        shapeID to verify whether you are in the same shape as the previous call.
        @param suppressOutOfOrder: Set this to False to check that all elements come back in
            sequence order, grouped together per shape ID, otherwise that isn't guaranteed.
        """
        raise NotImplementedError()
    
    def _checkShape(self, shapeEntry, suppressOutOfOrder=False):
        """
        Checks for out-of-order shape entries, as well as whether they are included or excluded.
        @param shapeEntry: A filled-out shape.ShapesEntry. Make sure ShapesEntry.shapeID is cast to the
            supported type; use shapeIDConstructor().
        @type shapeEntry: shape.ShapesEntry
        @return True if the worked properly, or false if there was an out of order entry or excluded.
        """
        # Process exclusions and inclusions:
        if self.exclusionList is not None and shapeEntry.shapeID in self.exclusionList:
            return False
        if self.inclusionList is not None and shapeEntry.shapeID not in self.inclusionList:
            return False
        
        if not suppressOutOfOrder:
            if self._currentShapeID is not None:
                if self._currentShapeID != shapeEntry.shapeID: 
                    if shapeEntry.shapeID in self.encounteredIDs:
                        # There's a shape ID that has unexpectedly appeared that had been encountered before!
                        print("WARNING: Shape ID " + str(shapeEntry.shapeID) + " has appeared unexpectedly in the shape data input.",
                            file=sys.stderr)
                        return False
                else:
                    if shapeEntry.shapeSeq <= self._prevSeq:
                        # Shape sequence is out of order!
                        print("WARNING: Shape ID " + str(shapeEntry.shapeID) + ", seq. " + str(shapeEntry.shapeSeq) + " is out of order.",
                            file=sys.stderr)
                        return False 
        self._currentShapeID = shapeEntry.shapeID
        if shapeEntry.shapeID not in self.encounteredIDs:
            self.encounteredIDs.add(shapeEntry.shapeID)
        self._prevSeq = shapeEntry.shapeSeq
        return True
    
    def readRoutes(self, shapes):
        """
        readRoutes generates a series of routes that are to correspond with the dataset this object points to.
        This default implementation creates a Route for each Shape
        @type shapes: shape.Shapes
        @return A map from routeID to a shape.Route object, or empty dict if there is none.
        @rtype shape.Routes
        """
        routes = shape.Routes()

        ctr = 0
        for shapeID in shapes:
            routes[ctr] = shapes.Route(ctr, shapeID, "")
            ctr += 1
        
        return routes
    
    def readTrips(self, shapes, routes, unusedShapeIDs=set(), restrictService=set()):
        """
        readTrips retrieves the trip information that corresponds with this dataset. The default implementation
        creates one trip per shape.
        @type shapes: dict<int, list<ShapesEntry>>
        @type routes: dict<int, RoutesEntry>
        @type unusedShapeIDs: set<?>
        @type restrictService: set<string>
        @return A map of tripID to shape.Trip records, as well as a list of unused trip IDs
        @rtype (shape.Trips, set<int>)
        """
        trips = shape.Trips()
        unusedTripIDs = set()
        "@type unusedTripIDs: set<?>"

        ctr = 0
        for routeID in routes:
            trips[ctr] = shapes.Trip(ctr, routes[routeID], routes[routeID].shortName, shapes[routes[routeID].shortName])
                
        return (ret, unusedTripIDs)
    
    def readStops(self):
        """
        readStops retrieves the stop information from this dataset. The default implementation returns an empty mapping.
        @return A mapping of stopID to a StopsEntry.
        @rtype shape.Stops
        """
        return shape.Stops()
        
    def readStopTimes(self, trips, stops, unusedTripIDs):
        """
        readStopTimes retrieves stop times for this dataset. The default implementation will create an empty list for
        each Trip ID.
        @type trips: Trips
        @type stops: Stops
        @type unusedTripIDs: set<?>
        @return A map of TripsEntry to a list of stop entries plus the start and end times
        @rtype StopTimes
        """
        stopTimes = shape.StopTimes()

        for tripID in trips:
            stopTimes[tripID] = list() # Fake the system by having no stops defined.
            
        return stopTimes





def fillShapes(filePath, gps):
    """
    fillShapes retrieves the shape information from a shape file and returns a list of shape entries.
    @type filePath: str
    @type gps: gps.GPS
    @return A map of shape_id to a list of shape entries
    @rtype dict<int, Shape>
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
                newEntry = shape.ShapesEntry(int(lineElems[0]), int(lineElems[3]), float(lineElems[1]),
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

def fillTrips(filePath, shapes, routes, unusedShapeIDs = set(), restrictService = set()):
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
