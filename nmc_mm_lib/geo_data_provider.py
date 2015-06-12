"""
geo_data_provider.py: Abstract class for functionality that has to do with reading in
underlying geographic topology.
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
import nmc_mm_lib.geo_data_providers, sys, inspect

importers = []
validImporter = None

class GeoDataProvider:
    def __init__(self):
        pass
    
    def provideCmdLineOpts(self, options):
        """
        Stuffs new command-line options into the given argparse object.
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
    
    def readData(self):
        """
        Using information from the command-line options, reads in data from the geographic data
        source that's handled by this class.
        @return A graph representing the geographic data
        @rtype graph.GraphLib
        """
        raise NotImplementedError()
    
def loadGeoDataProviders(options):
    """
    Dynamically load all of the importers that are available in the geo_data_providers directory
    and inform the command line processor.
    @type options: argparse 
    """
    global importers
    modules = []
    for module in nmc_mm_lib.geo_data_providers.__all__:
        __import__(module)
        modules.append(sys.modules[module])
        
    for module in modules:
        for name, obj in inspect.getmembers(module):
            if inspect.isclass(obj) and issubclass(obj, GeoDataProvider):
                newImporter = obj()
                importers.append(newImporter)
                newImporter.provideCmdLineOpts(options)

def readData(args):
    """
    Reads in geographic input data using the importer that was found in checkCmdLineOpts().
    @type args: argparse.Namespace
    @return A graph representing the geographic data, or None if there is a problem.
    @rtype graph.GraphLib
    """
    global importers, validImporter
    if not _readCmdLineOpts(args):
        return None
    if validImporter is None:
        print("ERROR: No valid importer for underlying geographic topology was found.", file=sys.stderr)
        return None
    graphData = validImporter.readData()
    return graphData

def _readCmdLineOpts(args):
    """
    Checks the given command-line options and determines if there are sufficient options
    for successfully loading in data for any supported subclasses. Call loadGeoDataProviders first.
    @type args: argparse.Namespace
    @return True if sufficient command-line options are provided for this class
    @rtype boolean
    """
    global importers, validImporter
    validCount = 0
    validImporter = None
    for importer in importers:
        if importer.checkCmdLineOpts(args):
            validImporter = importer
            validCount += 1
    if validCount == 0:
        print("ERROR: Insufficient command line parameters have been given to load underlying geographic topology.", file=sys.stderr)
        return False 
    if validCount > 1:
        print("ERROR: Too many command line parameters were given such that there are conflicts among supported geographic"
              " topology importers. Only one import style may be used.", file=sys.stderr)
        return False
    validImporter.readCmdLineOpts(args)
    return True
