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
from typing import Hashable, NamedTuple
import logging

class OSMReader:
    """
    OSMReader is a class that handles reading OSM data from the Overpass API
    and converting it into a graph.Map object.
    """
    class OverpassBounds(NamedTuple):
        """
        Definition for Overpass bounding box to limit queries
        """
        minLat: float
        maxLat: float
        minLon: float
        maxLon: float
        stepsNS: int = 1
        stepsEW: int = 1
        overlap: float = 0.05

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
    nodeCache: dict[Hashable, Intersection]
    waySets: dict[Hashable, dict[Way, bool]]

    def __init__(self, endpoint: str, bounds: OverpassBounds) -> None:
        self.endpoint = endpoint
        self.bounds = bounds

    def geoRead(self) -> graph.Map:
        """
        Queries Overpass API for OSM data within the specified bounding box
        """
        self.nodeCache = {}
        self.waySets = {}
        
        blockWidth = (self.bounds.maxLon - self.bounds.minLon) / self.bounds.stepsEW
        blockHeight = (self.bounds.maxLat - self.bounds.minLat) / self.bounds.stepsNS
        for vStep in range(self.bounds.stepsNS):
            for hStep in range(self.bounds.stepsEW):
                lowCoords = (self.bounds.minLat + blockHeight * vStep - blockHeight * self.bounds.overlap, self.bounds.minLon + blockWidth * hStep - blockWidth * self.bounds.overlap)
                highCoords = (self.bounds.minLat + blockHeight * (vStep + 1) + blockHeight * self.bounds.overlap, self.bounds.minLon + blockWidth * (hStep + 1) + blockWidth * self.bounds.overlap)
                print("Getting (%.4f,%.4f)-(%.4f,%.4f)..." % (lowCoords[0], lowCoords[1], highCoords[0], highCoords[1]))
                self.getChunk(lowCoords, highCoords)

        print("Sorting through final geometry...")
        intList = []
        for nodeID, node in self.nodeCache.items():
            nonMotorwayCnt = 0
            motorwayCnt = 0
            endCnt = 0
            for way, endFlag in self.waySets[nodeID].items():
                if way.type.startswith("motorway"):
                    motorwayCnt += 1
                else:
                    nonMotorwayCnt += 1
                if endFlag:
                    endCnt += 1
            if node.signal or (motorwayCnt + nonMotorwayCnt > 1 and not (endCnt == 2 and motorwayCnt + nonMotorwayCnt == 2)):
                motorwayFlag = node.junction or nonMotorwayCnt == 0
                intList.append(OSMReader.Intersection(lat=node.lat, lon=node.lon, signal=node.signal, junction=motorwayFlag)) 
        print("Number of intersections: %d" % len(intList))
            
        # Commit our geometry:
        logging.info("Committing geometry.")
        map.completeMap()
        return map

    def getChunk(self, lowCoords, highCoords):
        queryStr = f'[out:json];way({lowCoords[0]},{lowCoords[1]},{highCoords[0]},{highCoords[1]}){self.HIGHWAY_CLAUSE};(._;>;);out meta;'
        logging.info(f"Fetching from Overpass API ({lowCoords[0]:.3f}, {lowCoords[1]:.3f})-({highCoords[0]:.3f}, {highCoords[1]:.3f})")
        queryStr = urllib.parse.quote(queryStr)
        response = requests.get(self.endpoint + "?data=" + queryStr)
        response.raise_for_status()
        result = response.json()
        
        # Pass #1: Nodes:
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
        
        # Pass #2: Ways:
        wayCount = 0
        for element in result["elements"]:
            if "type" in element and element["type"] == "way":
                if "tags" in element and "highway" in element["tags"]:
                    ourName = element["tags"]["name"].strip().upper() if "name" in element["tags"] else "none"
                    way = OSMReader.Way(type=element["tags"]["highway"], name=ourName)
                    wayCount += 1
                    numNodes = len(element["nodes"])
                    for nodeID in element["nodes"]:
                        if nodeID in self.waySets:
                            self.waySets[nodeID][way] = nodeID == 0 or nodeID == numNodes - 1
        logging.info(f"New nodes: {nodeCount}; New ways: {wayCount}.")
        return nodeCount
    