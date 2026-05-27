"""
avl_osm_sample.py resolves AVL (Automatic Vehicle Location) tracks collected
from CapMetro (Austin, TX area transit agency)to an OpenStreetMap network.
Map data comes from the Overpass API. The sample captures two rides through
Austin, TX.
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

from support import osm_overpass
from shapely import from_wkt
from nmc_mm_lib import graph, path_engine, dump_io, reporter
from typing import Final, Any, Hashable, Generator
from datetime import datetime
import csv
import os
import sys
import logging

STEP_SIZE: Final[float] = 20.0

# Configure logging to use stdout
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
    stream=sys.stdout,
)

AVL_PATH: Final[str] = os.path.join("samples", "avl", "capmetro_rapid_20241101.csv")
OVERPASS_API: Final[str] = "https://overpass-api.de/api/interpreter"
OVERPASS_BOUNDS: Final[osm_overpass.OSMReader.OverpassBounds] = (
    osm_overpass.OSMReader.OverpassBounds(
        # Depicts the GPS bounding box for the greater Austin, TX metro area:
        minLat=29.582,
        maxLat=30.672,
        minLon=-98.050,
        maxLon=-97.513,
        stepsNS=3,  # How many north-south chunks to request
        stepsEW=3,  # How many east-west chunks to request
        overlap=0.05,  # Degrees of overlaps in rectangular chunks
    )
)

# Create OSM base map:
osmReader = osm_overpass.OSMReader(OVERPASS_API, OVERPASS_BOUNDS)
osmReader.geoRead()
map = graph.Map(workingCRS="EPSG:3081")  # Use Texas system in meters
osmReader.addToMap(map)
map.completeMap()


# Now, load in AVL bus tracks:
def avlRead(filename: str) -> Generator[dict[str, Any]]:
    """
    Reads AVL CSV file and yields line by line
    """
    with open(filename, mode="r", newline="") as fileHandle:
        csvReader = csv.DictReader(fileHandle)
        for fileLine in csvReader:
            yield fileLine


AVLCollection = dict[str, dict[datetime, Any]]

avl: AVLCollection = {}
for fileLine in avlRead(AVL_PATH):
    if fileLine["trip_id"] not in avl:
        avl[fileLine["trip_id"]] = {}
    datetimeKey = datetime.fromisoformat(fileLine["avl_timestamp"])
    avl[fileLine["trip_id"]][datetimeKey] = {
        "lat": float(fileLine["latitude"]),
        "lon": float(fileLine["longitude"]),
        "timestamp": datetimeKey,
        "stop_id": fileLine["stop_id"],
        "current_status": fileLine["current_status"],
        "current_stop_seq": fileLine["current_stop_seq"],
        "speed": float(fileLine["speed"]),
        "bearing": float(fileLine["bearing"]) if fileLine["bearing"] else None,
    }

# Express trackpoints derived from AVL data in terms of map:
avlTracks: dict[Hashable, tuple[graph.Trackpoint, ...]] = {}
for tripID, avlEntries in avl.items():
    avlTracks[tripID] = tuple(
        map.makeTrackpoint(
            lonHoriz=avlEntry["lon"], latVert=avlEntry["lat"], ident=tripID, seq=index
        )
        for index, avlEntry in enumerate(avlEntries.values())
    )

# Run path match for each AVL track:
matchedPaths: dict[Hashable, list[path_engine.PathEnd]] = {}
pathEngine = path_engine.PathEngine(
    path_engine.PathEngine.Params(
        maxHops=16
    )
)
for tripID, avlTrack in avlTracks.items():
    '''DEBUGGING'''
    if tripID == "2811287_12440":
        continue
    ''''''
    logging.info(f"AVL Trip ID {tripID} with {len(avlTrack)} trackpoints")
    path: list[path_engine.PathEnd] | None = pathEngine.constructPath(avlTrack, map)
    if path is not None:
        matchedPaths[tripID] = path

# Get a list of output trackpoints:
trackpointLists: dict[Hashable, list[reporter.OutputTrackpoint]] = {}
for tripID, treeNodes in matchedPaths.items():
    trackpointLists[tripID] = reporter.prepareTrackpath(
        map, treeNodes, intermediary=True, increment=STEP_SIZE
    )

# Output matched results:
with open("avl_matched.csv", mode="wt") as outputFile:
    dump_io.dumpStandardInfo(trackpointLists, outputFile, includeHeader=True)

# Explain series of streets for each AVL Trip ID:
for tripID in matchedPaths.keys():
    logging.info(f"AVL Trip ID: {tripID}")
    streetName = ""
    pathPoint: path_engine.PathEnd
    for pathPoint in matchedPaths[tripID]:
        link: graph.Map.LinkRecord
        for link in [pathPoint.pointOnLink.link] + pathPoint.routeInfo:
            newStreetName = link.data["name"].strip().lower()
            if newStreetName != streetName:
                logging.info(f'  {link.data["name"].strip()}')
                streetName = newStreetName
