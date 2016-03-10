"""
vista_network.py is an underlying node-link data representation of a
    roadway network. This requires psycopg to read a PostgreSQL DB
    of the format defined in file vista-def.sql.
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
from nmc_mm_lib import graph
import sys, psycopg2

def connect(dbServer, userName, password, networkName):
    """
    Connects to the VISTA database.
    @type dbServer: str
    @type userName: str
    @type password: str
    @type networkName: str
    @rtype psycopg2.connection
    """
    dbName = userName + "_" + networkName
    database = psycopg2.connect(host = dbServer, user = userName, password = password, database = dbName)
    return database

def fillGraph(database, useDirectDist=True):
    """
    fillGraph fills up the Graph structure from the VISTA database.
    @type database: psycopg2.connection
    @return A Graph representing the VISTA network model
    @rtype graph.GraphLib
    """
    cursor = database.cursor()
    
    # Step 1: Figure out the geographic center of the network and create the Graph:  
    cursor.execute('SELECT AVG(x), AVG(y) FROM nodes WHERE type = 1')
    row = cursor.fetchone()
    graphLib = graph.GraphLib(row[1], row[0], useDirectDist)
    
    # Step 2: Fill out the nodes:
    cursor.execute('SELECT id, x, y FROM nodes WHERE type = 1')
    for row in cursor:
        node = graph.GraphNode(row[0], row[2], row[1])
        graphLib.addNode(node)
    
    # Step 3: Fill out the links:
    cursor.execute("SELECT id, source, destination, length FROM linkdetails WHERE type = 1")
    for row in cursor:
        if (row[1] not in graphLib.nodeMap) or (row[2] not in graphLib.nodeMap):
            print("WARNING: Link %d has bad Node IDs %d and/or %d" % (row[0], row[1], row[2]), file = sys.stderr)
            continue 
        link = graph.GraphLink(row[0], graphLib.nodeMap[row[1]], graphLib.nodeMap[row[2]])
        link.distance = row[3] # This is the reported distance, but may be replaced if GraphLib.useDirectDist is true.
        graphLib.addLink(link)
        
    # There we are.
    return graphLib
