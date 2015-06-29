"""
shape_data_provider.py: Abstract class for functionality that has to do with reading in
the shape data that's to be mapped to underlying topology.
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
import shape, nmc_mm_lib.shape_data_providers, sys, inspect

importers = []

class ShapeDataProvider:
    """
    The scheme for ShapeDataProvider is to keep the shape import abstract (as we always need
    shape data for anything we're doing) but to create fake routes and trips such that each
    trip corresponds with one shape, and each route corresponds with one trip.
    @todo This is a good opportunity to detect similarities in shapes and consolidate.
    
    @ivar exlusionList: set<str> of shape IDs to exclude; set when getShapeImporter() is called.
    @ivar inclusionList: set<str> of shape IDs to exclude; set when getShapeImporter() is called.
    @ivar _encounteredIDs: set<str> keeps track of the IDs that have been encountered.
    @ivar _currentShapeID: str that is the current Shape ID when keeping track of order.
    @ivar _currentShape: shape.Shape that is the current shape.
    @ivar _prevSeq: int that is used to keep track of out-of-order entries in readShape().
    @ivar gps: The resident GPS object used with this importer
    """
    def __init__(self):
        # These lists are here in case we get to the point that we need to maintain two different
        # shape sets and each needs its own list of exclusions or inclusions.
        self.exclusionList = set()
        self.inclusionList = set()
        self._encounteredIDs = None
        self._currentShapeID = None
        self._currentShape = None
        self._prevSeq = -1
        self.gps = None
    
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
        
        for shapeEntry in self.readShapesEntries(suppressOutOfOrder=True):
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
        @return True if the worked properly, or False if there was an out of order entry or excluded.
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
        @type shapes: shape.Shapes
        @type routes: shape.Routes
        @type unusedShapeIDs: set<str>
        @type restrictService: set<string>
        @return A map of tripID to shape.Trip records, as well as a list of unused trip IDs
        @rtype (shape.Trips, set<str>)
        """
        trips = shape.Trips()
        unusedTripIDs = set()
        "@type unusedTripIDs: set<str>"

        ctr = 0
        for routeID in routes:
            trips[str(ctr)] = shapes.Trip(ctr, routes[routeID], routes[routeID].shortName, "")
                
        return (trips, unusedTripIDs)
    
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
        @type unusedTripIDs: set<str>
        @return A map of TripsEntry to a list of stop entries plus the start and end times
        @rtype StopTimes
        """
        stopTimes = shape.StopTimes()

        for tripID in trips:
            stopTimes[tripID] = list() # Fake the system by having no stops defined.
            
        return stopTimes
    
def loadShapeProviders(options):
    """
    Dynamically load all of the importers that are available in the geo_data_providers directory
    and inform the command line processor.
    @type options: argparse 
    """
    global importers
    modules = []
    for module in nmc_mm_lib.shape_data_providers.__all__:
        __import__(module)
        modules.append(sys.modules[module])
        
    for module in modules:
        for name, obj in inspect.getmembers(module):
            if inspect.isclass(obj) and issubclass(obj, GeoDataProvider):
                newImporter = obj()
                importers.append(newImporter)
                newImporter.provideCmdLineOpts(options)
                
    # Finally, provide a couple of command-line options that are supported among all:
    options.add_argument("-sx", nargs=1, type=str, metavar="EXCLUDEFILENAME", dest="excludeFilename",
        help="Specifies a file that lists all shape IDs to be excluded")
    options.add_argument("-si", nargs=1, type=str, metavar="INCLUDEFILENAME", dest="includeFilename",
        help="Specifies a file that lists all shape IDs to be included")

def getValidShapeProvider(args, gps):
    """
    Returns an object that can be used to retrieve shape data (and routes, etc). Preloads that object with
    the exclusion and inclusion lists if specified.
    @type args: argparse.Namespace
    @type gps: gps.GPS
    @return The ShapeDataProvider instance that is implied by the command-line parameters that had been set.
    @rtype ShapeDataProvider
    """
    global importers
    
    validCount = 0
    validImporter = None
    for importer in importers:
        if importer.checkCmdLineOpts(args):
            validImporter = importer
            validCount += 1
    if validCount == 0:
        print("ERROR: Insufficient command line parameters have been given to load shape data.", file=sys.stderr)
        return False 
    if validCount > 1:
        print("ERROR: Too many command line parameters were given such that there are conflicts among supported shape"
              " data importers. Only one import style may be used.", file=sys.stderr)
        return False
    validImporter.readCmdLineOpts(args)
    
    if args.excludeFilename is not None and args.includeFilename is not None:
        print("WARNING: It is ideal not to use both exclusion lists and inclusion lists at the same time.", file=sys.stderr)
    if args.excludeFilename is not None:
        with open(args.excludeFilename, 'r') as inFile:
            for fileLine in inFile:
                id = fileLine.strip()
                if len(id) > 0:
                    validImporter.exclusionList.add(id)
    if args.includeFilename is not None:
        with open(args.includeFilename, 'r') as inFile:
            for fileLine in inFile:
                id = fileLine.strip()
                if len(id) > 0:
                    validImporter.inclusionList.add(id)

    validImporter._encounteredIDs = set()        
    validImporter._currentShapeID = None
    validImporter._prevSeq = -1
    validImporter.gps = gps
    return validImporter
