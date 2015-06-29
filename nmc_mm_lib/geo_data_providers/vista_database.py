"""
vista_database.py is an underlying node-link data representation of a
    roadway network coded within a database per the VISTA traffic simulator.
    This requires psycopg to read a PostgreSQL DB of the format defined in
    file vista-def.sql.
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
import nmc_mm_lib.graph, nmc_mm_lib.geo_data_provider, sys
# Also requires psycopg2, but the import is down in readData() to avoid errors
# if psycopg2 is not installed. 

class VISTADatabase(nmc_mm_lib.geo_data_provider.GeoDataProvider):
    """
    @ivar connection: psycopg2.connection object
    @ivar serverName: str representing the database server name
    @ivar userName: str representing the database user name
    @ivar password: str for the database password
    @ivar dbName: str for the name of the database
    """
    def __init__(self):
        "@type self.connection: psycopg2.connection"
        self.connection = None
        self.serverName = ""
        self.userName = ""
        self.password = ""
        self.dbName = ""
        
        super(nmc_mm_lib.geo_data_provider.GeoDataProvider, self).__init__()
    
    def provideCmdLineOpts(self, options):
        """
        Stuffs new command-line options into the given argparse object.
        @type options: argparse
        """
        options.add_argument("-ds", nargs=1, type=str, metavar="DBSERVER", dest="serverName", default="localhost",
            help="If using a VISTA database for underlying topology input, specifies the database server")
        options.add_argument("-du", nargs=1, type=str, metavar="DBUSERNAME", dest="userName",
            help="If using a VISTA database for underlying topology input, specifies a user name")
        options.add_argument("-dp", nargs=1, type=str, metavar="DBPASSWORD", dest="password", default="",
            help="If using a VISTA database for underlying topology input, specifies a database password")
        options.add_argument("-dd", nargs=1, type=str, metavar="DBNAME", dest="databaseName",
            help="If using a VISTA database for underlying topology input, specifies a database name directly. Cannot be used with -dn.")
        options.add_argument("-dn", nargs=1, type=str, metavar="NETWORKNAME", dest="networkName",
            help='If using a VISTA database for underlying topology input, specifies a network name. The resulting database name '
            'will then be "DBUSERNAME_NETWORKNAME". Cannot be used with -dd.')

    def checkCmdLineOpts(self, args):
        """
        Checks the given command-line options and determines if there are sufficient options
        for successfully loading in data for this class.
        @type args: argparse.Namespace
        @return True if sufficient command-line options are provided for this class
        @rtype boolean
        """
        return args.databaseName is not None or args.networkName is not None
    
    def readCmdLineOpts(self, args):
        """
        Reads in all relevant command-line options and returns false if there are insufficient
        options provided.
        @type args: argparse.Namespace
        @return False if there are insufficient options provided, or True if successful.
        @rtype boolean
        """
        self.serverName = args.serverName
        self.userName = args.userName
        if self.userName is None:
            print("ERROR: A database user name parameter -du is required if -dd or -dn were specified.", file=sys.stderr)
            return False; 
        if args.databaseName is not None:
            self.dbName = args.databaseName
        if args.networkName is not None:
            if self.dbName is not None:
                print("ERROR: Both the -dd and -dn parameters cannot be used together.", file=sys.stderr)
                return False; 
            self.dbName = args.userName + "_" + args.networkName
        if len(self.dbName.strip()) == 0:
            print("ERROR: No database name was specified either through the -dd or -dn parameters.", file=sys.stderr)
            return False; 
        self.password = args.password
        return True;
    
    def readData(self):
        """
        Using information from the command-line options, reads in data from the geographic data
        source that's handled by this class.
        @return A graph representing the geographic data
        @rtype graph.GraphLib
        """
        print("INFO: Reading geographic input from DB server '%s', database '%s'..." % (self.serverName, self.dbName), file=sys.stderr)
        import psycopg2
        self.connection = psycopg2.connect(host = self.serverName, user = self.userName, password = self.password,
            database = self.dbName)
        return self._fillGraph()
    
    def _fillGraph(self):
        """
        fillGraph fills up the Graph structure from the VISTA database.
        @return A Graph representing the VISTA network model
        @rtype graph.GraphLib
        """
        cursor = self.connection.cursor()
        
        # Step 1: Figure out the geographic center of the network and create the Graph:  
        cursor.execute('SELECT AVG(x), AVG(y) FROM nodes WHERE type = 1')
        row = cursor.fetchone()
        graphLib = nmc_mm_lib.graph.GraphLib(row[1], row[0])
        
        # Step 2: Fill out the nodes:
        cursor.execute('SELECT id, x, y FROM nodes WHERE type = 1')
        for row in cursor:
            node = nmc_mm_lib.graph.GraphNode(row[0], row[2], row[1])
            graphLib.addNode(node)
        
        # Step 3: Fill out the links:
        cursor.execute("SELECT id, source, destination, length FROM linkdetails WHERE type = 1")
        for row in cursor:
            if (row[1] not in graphLib.nodeMap) or (row[2] not in graphLib.nodeMap):
                print("WARNING: Link %d has bad Node IDs %d and/or %d" % (row[0], row[1], row[2]), file = sys.stderr)
                continue 
            link = nmc_mm_lib.graph.GraphLink(row[0], graphLib.nodeMap[row[1]], graphLib.nodeMap[row[2]])
            #link.distance = row[3] # Use reported distance rather than measured distance.
            graphLib.addLink(link)
            
        # There we are.
        return graphLib
