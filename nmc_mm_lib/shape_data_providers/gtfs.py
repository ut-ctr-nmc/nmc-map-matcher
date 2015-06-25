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
import os, sys, csv
from datetime import datetime, timedelta

class GTFS(shape_data_provider.ShapeDataProvider):
    """
    Represents a GTFS file set that represents shape data and associated transit data.
    
    @ivar filePath: str that is the path to the GTFS directory
    @ivar shapeFilename: str that holds the name of the shapefile; useful for something other than "shapes.txt"
    """
    def __init__(self):
        self.filePath = None
        self.shapeFilename = None
        
        super(shape_data_provider.ShapeDataProvider, self).__init__()
    
    def provideCmdLineOpts(self, options):
        """
        Stuffs new command-line options into the given argparse object. This performs the
        shapefile-specific command line specification setup. 
        @type options: argparse
        """
        options.add_argument("-gp", nargs=1, type=str, metavar="GTFSPATH", dest="gtfsPath",
            help="If using GTFS for shape/transit input, specifies the path to the GTFS set")
        options.add_argument("-gs", nargs=1, type=str, metavar="GTFSSHAPES", dest="gtfsShapes",
            help="If using GTFS for shape/transit input, specifies an alternate GTFS shapes filename")

    def checkCmdLineOpts(self, args):
        """
        Checks the given command-line options and determines if there are sufficient options
        for successfully loading in data for this class.
        @type args: argparse.Namespace
        @return True if sufficient command-line options are provided for this class
        @rtype boolean
        """
        return args.gtfsPath is not None or args.gtfsShapes is not None
    
    def readCmdLineOpts(self, args):
        """
        Reads in all relevant command-line options and returns false if there are insufficient
        options provided.
        @type args: argparse.Namespace
        @return False if there are insufficient options provided, or True if successful.
        @rtype boolean
        """
        self.filePath = args.gtfsPath if args.gtfsPath is not None else os.getcwd()
        self.shapeFilename = args.gtfsShapes if args.gtfsShapes is not None else "shapes.txt"
        
        if not os.path.isabs(self.shapeFilename):
            self.shapeFilename = os.path.join(self.filePath, self.shapeFilename)
        return True
        
    def readShapesEntries(self, suppressOutOfOrder=False):
        """
        readShapesEntries returns a generator that gives out shape entries one by one.
        At each call, the next valid constructed shape.ShapesEntry is yielded. You must check the
        shapeID to verify whether you are in the same shape as the previous call.
        @param suppressOutOfOrder: Set this to False to check that all elements come back in
            sequence order, grouped together per shape ID, otherwise that isn't guaranteed.
        """
        loadedShapes = {}
        "@type loadedShapes: dict<str, shape.Shape>"
        with open(self.shapeFilename, 'r') as inFile:
            csvReader = csv.DictReader(inFile)
            
            firstRun = True
            for fileLine in csvReader:
                if firstRun:
                    # Sanity check:
                    if not (all(x in fileLine for x in ["shape_id", "shape_pt_lat", "shape_pt_lon", "shape_pt_sequence"])):
                        print("ERROR: The GTFS shapefile %s doesn't have the expected header." % self.shapeFilename, file=sys.stderr)
                        return None
                    firstRun = False
            
                if fileLine["shape_id"] not in loadedShapes:
                    loadedShapes[fileLine["shape_id"]] = shape.Shape(fileLine["shape_id"])
                newEntry = shape.ShapesEntry(loadedShapes[fileLine["shape_id"]], int(fileLine["shape_pt_sequence"]),
                            float(fileLine["shape_pt_lat"]), float(fileLine["shape_pt_lon"]))
                (newEntry.pointX, newEntry.pointY) = self.gps.gps2feet(newEntry.lat, newEntry.lng)
                if self._checkShape(newEntry, suppressOutOfOrder):
                    yield newEntry
    
    def readRoutes(self, shapes):
        """
        readRoutes generates a series of routes that are to correspond with the dataset this object points to.
        This default implementation creates a Route for each Shape
        @type shapes: shape.Shapes
        @rtype shape.Routes
        """            
        ret = shape.Routes()
        filename = os.path.join(self.filePath, "routes.txt") 
        with open(filename, 'r') as inFile:
            csvReader = csv.DictReader(inFile)
            
            firstRun = True
            for fileLine in csvReader:
                if firstRun:
                    # Sanity check:
                    if not (all(x in fileLine for x in ["route_id", "agency_id", "route_short_name", "route_long_name"])):
                        print("ERROR: The routes.txt file doesn't have the expected header.", file=sys.stderr)
                        return None
                    firstRun = False
            
                newEntry = shape.Route(fileLine["route_id"], fileLine["route_short_name"], fileLine["route_long_name"])
                ret[newEntry.routeID] = newEntry
        
        # Return the stops file contents:
        return ret
    
    def readTrips(self, shapes, routes, unusedShapeIDs=set(), restrictService=set()):
        """
        readTrips retrieves the trip information that corresponds with this dataset.
        @type shapes: shape.Shapes
        @type routes: shape.Routes
        @type unusedShapeIDs: set<str>
        @type restrictService: set<str>
        @return A dict of tripID to shape.Trip records, as well as a list of unused trip IDs
        @rtype (shape.Trips, set<str>)
        """
        ret = shape.Trips()
        unusedTripIDs = set()
        "@type unusedTripIDs: set<str>"
        filename = os.path.join(self.filePath, "trips.txt") 
        with open(filename, 'r') as inFile:
            csvReader = csv.DictReader(inFile)

            firstRun = True
            for fileLine in csvReader:
                if firstRun:
                    # Sanity check:
                    if not (all(x in fileLine for x in ["route_id", "service_id", "trip_id", "trip_headsign", "trip_short_name",
                                                        "direction_id", "block_id", "shape_id"])):
                        print("ERROR: The trips.txt file doesn't have the expected header.", file=sys.stderr)
                        return None
                    firstRun = False

                shapeID = fileLine["shape_id"]
                routeID = fileLine["route_id"]
                serviceID = fileLine["service_id"]
                tripID = fileLine["trip_id"]
                if (shapeID in unusedShapeIDs) or ((len(restrictService) > 0) and (serviceID not in restrictService)):
                    unusedTripIDs.add(tripID)
                else:
                    if shapeID not in shapes:
                        print("WARNING: GTFS Trips file expects undefined shape ID " + shapeID, file=sys.stderr)
                        unusedTripIDs.add(tripID)                        
                    else:
                        if routeID not in routes:
                            print("WARNING: GTFS Trips file expects undefined route ID " + routeID, file=sys.stderr)
                            unusedTripIDs.add(tripID)                        
                        else:
                            newEntry = shape.Trip(tripID, routes[routeID], fileLine["trip_headsign"], shapes[shapeID])
                            ret[newEntry.tripID] = newEntry
                                
        # Return the trips file contents:
        return (ret, unusedTripIDs)

    def readStops(self):
        """
        readStops retrieves the stop information from this dataset. The default implementation returns an empty mapping.
        @return A mapping of stopID to a StopsEntry.
        @rtype shape.Stops
        """
        ret = shape.Stops()
        filename = os.path.join(self.filePath, "stops.txt") 
        with open(filename, 'r') as inFile:
            csvReader = csv.DictReader(inFile)

            firstRun = True
            for fileLine in csvReader:
                if firstRun:
                    # Sanity check:
                    if not (all(x in fileLine for x in ["stop_id", "stop_code", "stop_name", "stop_desc", "stop_lat", "stop_lon"])):
                        print("ERROR: The stops.txt file doesn't have the expected header.", file=sys.stderr)
                        return None
                    firstRun = False
            
                newEntry = shape.Stop(fileLine["stop_id"], fileLine["stop_name"], float(fileLine["stop_lat"]), float(fileLine["stop_lon"]))
                (newEntry.pointX, newEntry.pointY) = self.gps.gps2feet(newEntry.gpsLat, newEntry.gpsLng)
                ret[newEntry.stopID] = newEntry
                        
        # Return the stops file contents:
        return ret
        
    def readStopTimes(self, trips, stops, unusedTripIDs):
        """
        readStopTimes retrieves stop times for this dataset. The default implementation will create an empty list for
        each Trip ID.
        @type trips: Trips
        @type stops: Stops
        @type unusedTripIDs: set<str>
        @return A map of TripsEntry to a list of stop entries plus the start and end times
        @rtype StopTimes
        """
        stopTimes = shape.StopTimes()
        
        filename = os.path.join(self.filePath, "stop_times.txt")
        with open(filename, 'r') as inFile:
            csvReader = csv.DictReader(inFile)

            firstRun = True
            for fileLine in csvReader:
                if firstRun:
                    # Sanity check:
                    if not (all(x in fileLine for x in ["trip_id", "arrival_time", "departure_time", "stop_id", "stop_sequence"])):
                        print("ERROR: The stops_times.txt file doesn't have the expected header.", file=sys.stderr)
                        return None
                    firstRun = False
            
                if fileLine["tripID"] not in unusedTripIDs:
                    if not fileLine["tripID"] in trips:
                        print("WARNING: GTFS Stop Times file expects undefined trip ID " + fileLine["tripID"], file=sys.stderr)
                        continue
                    
                    # Split apart time string because hours can roll around.
                    timeElems = fileLine["arrival_time"].split(':')
                    timeHour = int(timeElems[0])
                    timeDays = timeHour / 24
                    timeHour = timeHour % 24
                    arrivalTime = datetime(1900, 1, 1, timeHour, int(timeElems[1]), int(timeElems[2]))
                    arrivalTime += timedelta(days = timeDays)

                    timeElems = fileLine["departure_time"].split(':')
                    timeHour = int(timeElems[0])
                    timeDays = timeHour / 24
                    timeHour = timeHour % 24
                    departureTime = datetime(1900, 1, 1, timeHour, int(timeElems[1]), int(timeElems[2]))
                    departureTime += timedelta(days = timeDays)
                    
                    if not fileLine["stop_id"] in stops:
                        print("WARNING: GTFS Stop Times file expects undefined stop ID " + fileLine["stop_id"], file=sys.stderr)
                        continue
                    newEntry = shape.StopTime(fileLine["tripID"], fileLine["stop_id"], int(fileLine["stop_sequence"]))
                    newEntry.arrivalTime = arrivalTime
                    newEntry.departureTime = departureTime
                    if newEntry.trip not in stopTimes:
                        stopTimes[newEntry.trip] = []
                    stopTimes[newEntry.trip].append(newEntry)
    
        # Return the stop times file contents:
        return stopTimes
