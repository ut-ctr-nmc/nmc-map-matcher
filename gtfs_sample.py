"""
gtfs_sample.py resolves a GTFS shapefile to an MPO downtown grid network
    series of links and outputs a CSV format of data showing the GTFS
    tracks with respect to the grid network.
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
import logging
from typing import Final, Any, Hashable
from nmc_mm_lib import graph, gtfs, path_engine
import csv, os

MPO_PATH: Final[str] = os.path.join("samples", "mpo")
GTFS_PATH: Final[str] = os.path.join("samples", "gtfs")

# Grab MPO model: we'll derive the topography from node locations (nodes.csv)
# and connectivity (cnx.csv). The links.csv file is there to help with
# visualization in a GIS program, but it isn't necessary here.
# TODO: This is much more compact with Pandas!
nodes: dict[Hashable, dict[str, Any]] = {}
filename = os.path.join(MPO_PATH, "small_atx_nodes.csv")
with open(filename, mode='r', newline='') as nodesFile:
    csvReader = csv.DictReader(nodesFile)
    for fileLine in csvReader:
        nodes[fileLine["id"]] = {"id": fileLine["id"],
                                 "lon": float(fileLine["lon"]),
                                 "lat": float(fileLine["lat"])}
cnxs: dict[Hashable, dict[str, Any]] = {}
filename = os.path.join(MPO_PATH, "small_atx_cnx.csv")
with open(filename, mode='r', newline='') as cnxsFile:
    csvReader = csv.DictReader(cnxsFile)
    for fileLine in csvReader:
        cnxs[fileLine["id"]] = {"id": fileLine["id"],
                                "source": nodes[fileLine["source"]],
                                "dest": nodes[fileLine["dest"]]}

# Create map of it:
map = graph.Map() # Use default GPS to Web Mercator scheme
for nodeID, node in nodes.items():
    # Define each node from lon/lat, ID, w/o optional metadata dict
    map.addNode(nodeID, node["lon"], node["lat"])
for linkID, cnx in cnxs.items():
    # Define each link by using node IDs. By saying that we hadn't specified
    # endpoints, GPS endpoints for each link are grabbed from the nodes:
    map.addLink(cnx["source"]["id"], cnx["dest"]["id"], linkID=linkID,
                hasEndpoints=False)
# Commit our geometry:
map.completeMap()

# Grab GTFS:
gtfsSet = gtfs.GTFSSet(GTFS_PATH)

# Express trackpoints derived from GTFS shapes in terms of map:
gtfsShapesTracks: dict[Hashable, tuple[graph.Map.Trackpoint, ...]] = {}
for shapeID, shapeEntries in gtfsSet.shapes.items():
    gtfsShapesTracks[shapeID] = tuple(map.makeTrackpoint(shapeEntry.lng,
        shapeEntry.lat, shapeID, shapeEntry.shapeSeq)
            for shapeEntry in shapeEntries)

# Run path match for each GTFS route:
matchedPaths: dict[Hashable, list[path_engine.PathEnd]] = {}
pathEngine = path_engine.PathEngine() # Use default match parameters
for shapeID, gtfsTrack in gtfsShapesTracks.items():
    logging.info(f"GTFS Shape ID {shapeID}:")
    path: list[path_engine.PathEnd] | None \
        = pathEngine.constructPath(gtfsTrack, map)
    if path is not None:
        matchedPaths[shapeID] = path

# Output matched results:
for shapeID, matchedPath in matchedPaths.items():
    with open(f"gtfs_matched_{shapeID}.csv", mode='wt') as outputFile:
        path_engine.dumpStandardInfo(matchedPath, outputFile)
