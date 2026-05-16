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

from collections import namedtuple

from shapely import from_wkt
import urllib.parse
import requests
from nmc_mm_lib import graph, path_engine, dump_io, reporter
from support import gtfs
from typing import Final, Any, Hashable, Generator, TypedDict, NamedTuple
from datetime import datetime
import csv
import os
import sys
import logging

STEP_SIZE: Final[float] = 20.0

# Configure logging to use stdout
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    stream=sys.stdout,
)

AVL_PATH: Final[str] = os.path.join("samples", "avl", "samples/avl/capmetro_rapid_20241101.csv")
OVERPASS_API: Final[str] = "https://overpass-api.de/api/interpreter"
OVERPASS_BOUNDS: Final[OSMReader.OverpassBounds] = {
    # Depicts the GPS bounding box for the greater Austin, TX metro area:
    "minLat": 29.582,
    "maxLat": 30.672,
    "minLon": -98.050,
    "maxLon": -97.513
}
STEPS_NS: Final[int] = 3 # How many north-south chunks to request
STEPS_EW: Final[int] = 3 # How many east-west chunks to request
OVERLAP: Final[float] = 0.05 # Degrees of overlaps in rectangular chunks

def avlRead(filename: str) -> Generator[dict[str, Any]]:
    """
    Reads AVL CSV file and yields line by line
    """
    filename = os.path.join(AVL_PATH, filename)
    with open(filename, mode="r", newline="") as fileHandle:
        csvReader = csv.DictReader(fileHandle)
        for fileLine in csvReader:
            yield fileLine



class OSMReader:
    """
    OSMReader is a class that handles reading OSM data from the Overpass API
    and converting it into a graph.Map object.
    """
    # Determine roadway types we're interested in:
    HIGHWAY_CLAUSE: str = '["highway"~"^(motorway|trunk|primary|secondary' \
        '|tertiary|motorway_link|trunk_link|primary_link|unclassified' \
        '|residential|living_street)$"]'
    nodeCache: dict[Hashable, Intersection] = {}
    waySets: dict[Hashable, dict[Way, bool]] = {}

    class OverpassBounds(TypedDict):
        """
        Definition for Overpass bounding box to limit queries
        """
        minLat: float
        maxLat: float
        minLon: float
        maxLon: float

    def __init__(self, endpoint: str, bounds: OverpassBounds) -> None:
        self.endpoint = endpoint
        self.bounds = bounds

    class Intersection(NamedTuple):
        lat: float
        lon: float
        signal: bool
        junction: bool
        midblock_sig: bool

    class Way(NamedTuple):
        type: str
        name: str

    def geoRead(self) -> graph.Map:
        """
        Queries Overpass API for OSM data within the specified bounding box
        """



        # Commit our geometry:
        map.completeMap()
        return map

    def getChunk(self, lowCoords, highCoords):
        queryStr = f'[out:json];way({lowCoords[0]},{lowCoords[1]},{highCoords[0]},{highCoords[1]}){self.HIGHWAY_CLAUSE};(._;>;);out meta;'
        print(f"Fetching from Overpass API ({lowCoords[0]:.3f}, {lowCoords[1]:.3f})-({highCoords[0]:.3f}, {highCoords[1]:.3f})")
        queryStr = urllib.parse.quote(queryStr)
        response = requests.get(self.endpoint + "?data=" + queryStr)
        response.raise_for_status()
        result = response.json()
        
        nodeCount = 0
        for element in result["elements"]:
            if "type" in element and element["type"] == "node" and element["id"] not in nodeCache:
                sigFlag = False
                junctFlag = False
                if "tags" in element and "highway" in element["tags"]:
                    sigFlag = element["tags"]["highway"] == "traffic_signals"
                    junctFlag = element["tags"]["highway"] == "motorway_junction"
                nodeCache[element["id"]] = OSMReader.Intersection(lat=element["lat"],
                                                            lon=element["lon"],
                                                            signal=sigFlag,
                                                            junction=junctFlag,
                                                            midblock_sig=None)
                waySets[element["id"]] = {} # That's way -> True if endpoint
                nodeCount += 1
        for element in result["elements"]:
            if "type" in element and element["type"] == "way":
                if "tags" in element and "highway" in element["tags"]:
                    ourName = element["tags"]["name"].strip().upper() if "name" in element["tags"] else "none"
                    way = OSMReader.Way(type=element["tags"]["highway"], name=ourName)
                    index = 0
                    numNodes = len(element["nodes"])
                    for nodeID in element["nodes"]:
                        if nodeID in waySets:
                            waySets[nodeID][way] = nodeID == 0 or nodeID == numNodes - 1
                        index += 1
        print("New nodes: %d." % nodeCount)

'''



def process():
    nodeCache = {}
    waySets = {}
    
    nodeCount = 0
    blockWidth = (CORNER_HIGH[1] - CORNER_LOW[1]) / STEPS_EW
    blockHeight = (CORNER_HIGH[0] - CORNER_LOW[0]) / STEPS_NS
    for vStep in range(STEPS_NS):
        for hStep in range(STEPS_EW):
            lowCoords = (CORNER_LOW[0] + blockHeight * vStep - blockHeight * OVERLAP, CORNER_LOW[1] + blockWidth * hStep - blockWidth * OVERLAP)
            highCoords = (CORNER_LOW[0] + blockHeight * (vStep + 1) + blockHeight * OVERLAP, CORNER_LOW[1] + blockWidth * (hStep + 1) + blockWidth * OVERLAP)
            print("Getting (%.4f,%.4f)-(%.4f,%.4f)..." % (lowCoords[0], lowCoords[1], highCoords[0], highCoords[1]))
            nodeCount = getChunk(nodeCache, waySets, lowCoords, highCoords)

    print("Sorting through final geometry...")
    intList = []
    for nodeID, node in nodeCache.items():
        nonMotorwayCnt = 0
        motorwayCnt = 0
        endCnt = 0
        for way, endFlag in waySets[nodeID].items():
            if way.type.startswith("motorway"):
                motorwayCnt += 1
            else:
                nonMotorwayCnt += 1
            if endFlag:
                endCnt += 1
        if node.signal or (motorwayCnt + nonMotorwayCnt > 1 and not (endCnt == 2 and motorwayCnt + nonMotorwayCnt == 2)):
            motorwayFlag = node.junction or nonMotorwayCnt == 0
            intList.append(Intersection(lat=node.lat, lon=node.lon, signal=node.signal, junction=motorwayFlag,
                                        midblock_sig=node.signal and motorwayCnt + nonMotorwayCnt < 2)) 
    print("Number of intersections: %d" % len(intList))
    
    print("Outputting CSV '%s'..." % OUTFILE)
    outHandle = open(OUTFILE, "w")
    csvWriter = csv.writer(outHandle)
    csvWriter.writerow(['lat', 'lon', 'signal', 'junction', 'midblock_sig'])
    for node in intList:
        csvWriter.writerow([node.lat, node.lon, tf(node.signal), tf(node.junction), tf(node.midblock_sig)])
    outHandle.close()
    print("Done.")
'''

AVLCollection = dict[str, dict[datetime, Any]]

avl: AVLCollection = {}
for fileLine in avlRead(AVL_PATH):
    if fileLine["trip_id"] not in avl:
        avl[fileLine["trip_id"]] = {}
    datetimeKey = datetime.fromisoformat(fileLine["timestamp"])
    avl[fileLine["trip_id"]][datetimeKey] = {
        "id": fileLine["id"],
        "lat": float(fileLine["lat"]),
        "lon": float(fileLine["lon"]),
        "timestamp": fileLine["timestamp"],
        "route_id": fileLine["route_id"],
        "shape_id": fileLine["shape_id"]
    }


'''
# Create map of it:
map = graph.Map()  # Use default GPS to Web Mercator scheme
for nodeID, node in nodes.items():
    # Define each node from lon/lat, ID, w/o optional metadata dict. Keep in
    # mind that we actually don't need these locations since we're using the
    # link geometry.
    map.addNode(nodeID, node["lon"], node["lat"])
for linkID, cnx in cnxs.items():
    # Define each link by using node IDs, but use link geometry. Also add in
    # the street name metadata:
    map.addLink(
        cnx["source"]["id"],
        cnx["dest"]["id"],
        controlPoints=links[linkID]["geog"],
        linkID=linkID,
        metadata=links[linkID],
    )
'''






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
        map, treeNodes, intermediary=True, increment=STEP_SIZE
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
