"""
osm_overpass.py: Code for using Overpass to read OpenStreetMap data
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
import urllib.parse
import requests
from nmc_mm_lib import graph
from typing import Hashable, TypedDict, NamedTuple
import logging

class OSMReader:
    """
    OSMReader is a class that handles reading OSM data from the Overpass API
    and converting it into a graph.Map object.
    """
    class OverpassBounds(TypedDict):
        """
        Definition for Overpass bounding box to limit queries
        """
        minLat: float
        maxLat: float
        minLon: float
        maxLon: float

    class Intersection(NamedTuple):
        """
        Internal record-keeping for OSMReader
        """
        lat: float
        lon: float
        signal: bool
        junction: bool

    class Way(NamedTuple):
        """
        Internal record-keeping for OSMReader
        """
        type: str
        name: str

    # Determine roadway types we're interested in:
    HIGHWAY_CLAUSE: str = '["highway"~"^(motorway|trunk|primary|secondary' \
        '|tertiary|motorway_link|trunk_link|primary_link|unclassified' \
        '|residential|living_street)$"]'
    nodeCache: dict[Hashable, Intersection] = {}
    waySets: dict[Hashable, dict[Way, bool]] = {}

    def __init__(self, endpoint: str, bounds: OverpassBounds) -> None:
        self.endpoint = endpoint
        self.bounds = bounds

    def geoRead(self) -> graph.Map:
        """
        Queries Overpass API for OSM data within the specified bounding box
        """



        # Commit our geometry:
        map.completeMap()
        return map

    def getChunk(self, lowCoords, highCoords):
        queryStr = f'[out:json];way({lowCoords[0]},{lowCoords[1]},{highCoords[0]},{highCoords[1]}){self.HIGHWAY_CLAUSE};(._;>;);out meta;'
        logging.info(f"Fetching from Overpass API ({lowCoords[0]:.3f}, {lowCoords[1]:.3f})-({highCoords[0]:.3f}, {highCoords[1]:.3f})")
        queryStr = urllib.parse.quote(queryStr)
        response = requests.get(self.endpoint + "?data=" + queryStr)
        response.raise_for_status()
        result = response.json()
        
        nodeCount = 0
        for element in result["elements"]:
            if "type" in element and element["type"] == "node" and element["id"] not in self.nodeCache:
                sigFlag = False
                junctFlag = False
                if "tags" in element and "highway" in element["tags"]:
                    sigFlag = element["tags"]["highway"] == "traffic_signals"
                    junctFlag = element["tags"]["highway"] == "motorway_junction"
                self.nodeCache[element["id"]] = OSMReader.Intersection(lat=element["lat"],
                                                            lon=element["lon"],
                                                            signal=sigFlag,
                                                            junction=junctFlag)
                self.waySets[element["id"]] = {} # That's way -> True if endpoint
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
