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

from nmc_mm_lib import graph, path_engine, dump_io, reporter
from support import mpo_read, gtfs
from typing import Final, Hashable
import os
import sys
import logging

REPORT_STEP_SIZE: Final[float] = 20.0

# Configure logging to use stdout
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    stream=sys.stdout,
)

MPO_PATH: Final[str] = os.path.join("samples", "mpo")
GTFS_PATH: Final[str] = os.path.join("samples", "gtfs", "small_atx")

# Grab tiny MPO map:
map: graph.Map = mpo_read.readMPOModel(
    os.path.join(MPO_PATH, "small_atx_nodes.csv"),
    os.path.join(MPO_PATH, "small_atx_links.csv"),
    os.path.join(MPO_PATH, "small_atx_cnx.csv"),
)

# Grab GTFS:
gtfsSet = gtfs.GTFSSet(GTFS_PATH)

# Express trackpoints derived from GTFS shapes in terms of map:
gtfsShapesTracks: dict[Hashable, tuple[graph.Trackpoint, ...]] = {}
for shapeID, shapeEntries in gtfsSet.shapes.items():
    gtfsShapesTracks[shapeID] = tuple(
        map.makeTrackpoint(shapeEntry.lng, shapeEntry.lat, shapeID, shapeEntry.shapeSeq)
        for shapeEntry in shapeEntries
    )

# Run path match for each GTFS route:
matchedPaths: dict[Hashable, list[path_engine.PathEnd]] = {}
pathEngine = path_engine.PathEngine()  # Use default match parameters
for shapeID, gtfsTrack in gtfsShapesTracks.items():
    logging.info(f"GTFS Shape ID {shapeID}:")
    path: list[path_engine.PathEnd] | None = pathEngine.constructPath(gtfsTrack, map)
    if path is not None:
        matchedPaths[shapeID] = path

# Get a list of output trackpoints:
trackpointLists: dict[Hashable, list[reporter.OutputTrackpoint]] = {}
for shapeID, treeNodes in matchedPaths.items():
    trackpointLists[shapeID] = reporter.prepareTrackpath(
        map, treeNodes, intermediary=True, increment=REPORT_STEP_SIZE
    )

# Output matched results:
with open("gtfs_small_matched.csv", mode="wt") as outputFile:
    dump_io.dumpStandardInfo(trackpointLists, outputFile, includeHeader=True)

# Explain series of streets for each Shape ID:
for shapeID in matchedPaths.keys():
    logging.info(f"GTFS Shape ID: {shapeID}")
    streetName = ("", "")
    pathPoint: path_engine.PathEnd
    for pathPoint in matchedPaths[shapeID]:
        link: graph.Map.LinkRecord
        for link in [pathPoint.pointOnLink.link] + pathPoint.routeInfo:
            newStreetName = (link.data["name"], link.data["dir"])
            if newStreetName != streetName:
                logging.info(f'  {link.data["name"]} going {link.data["dir"]}')
                streetName = newStreetName
