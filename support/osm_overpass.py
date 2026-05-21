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

    class Node(NamedTuple):
        """
        Internal record-keeping for OSMReader
        """

        id: Hashable
        lat: float
        lon: float
        signal: bool
        junction: bool

    class Way(NamedTuple):
        """
        Internal record-keeping for OSMReader
        """

        id: Hashable
        type: str
        name: str
        oneWay: bool
        motorway: bool
        nodes: tuple["OSMReader.Node", ...]

    # Determine roadway types we're interested in:
    HIGHWAY_CLAUSE: str = (
        '["highway"~"^(motorway|trunk|primary|secondary'
        "|tertiary|motorway_link|trunk_link|primary_link|unclassified"
        '|residential|living_street)$"]'
    )
    nodeCache: dict[Hashable, Node]
    waySet: set[Way]
    wayNodeLkp: dict[Hashable, set[Way]]

    def __init__(self, endpoint: str, bounds: OverpassBounds) -> None:
        self.endpoint = endpoint
        self.bounds = bounds

    def geoRead(self) -> None:
        """
        Queries Overpass API for OSM data within the specified bounding box
        """
        self.nodeCache = {}
        self.waySet = {*()}
        self.wayNodeLkp = {}

        logging.info(f"Fetching {self.bounds.stepsEW}x{self.bounds.stepsNS} OSM data")
        blockWidth = (self.bounds.maxLon - self.bounds.minLon) / self.bounds.stepsEW
        blockHeight = (self.bounds.maxLat - self.bounds.minLat) / self.bounds.stepsNS
        for vStep in range(self.bounds.stepsNS):
            for hStep in range(self.bounds.stepsEW):
                lowCoords = (
                    self.bounds.minLat
                    + blockHeight * vStep
                    - blockHeight * self.bounds.overlap,
                    self.bounds.minLon
                    + blockWidth * hStep
                    - blockWidth * self.bounds.overlap,
                )
                highCoords = (
                    self.bounds.minLat
                    + blockHeight * (vStep + 1)
                    + blockHeight * self.bounds.overlap,
                    self.bounds.minLon
                    + blockWidth * (hStep + 1)
                    + blockWidth * self.bounds.overlap,
                )
                self.getChunk(lowCoords, highCoords)

    def addToMap(self, map: graph.Map) -> None:
        """
        Utility for loading the OSM Nodes and Ways to a None/Link representation in given graph.Map
        """
        logging.info("Loading into map")
        addedNodes: set[OSMReader.Node] = {*()}

        # Traverse through nodes in each Way, breaking apart sections between
        # nodes into links:
        linkCount = 0
        reversedLinkCount = 0
        for way in (w for w in self.waySet if len(w.nodes) > 1):
            startIdx = 0
            for index, node in enumerate(way.nodes):
                endFlag = False
                if index == len(way.nodes) - 1 or len(self.wayNodeLkp[node.id]) > 1:
                    endFlag = True
                if index == 0 or endFlag:
                    if node not in addedNodes:
                        map.addNode(
                            nodeID=node.id,
                            lonHoriz=node.lon,
                            latVert=node.lat,
                            metadata=node._asdict(),
                        )
                        addedNodes |= {node}
                    if endFlag and index > 0:
                        controlPoints = [
                            (n.lon, n.lat) for n in way.nodes[startIdx : index + 1]
                        ]
                        map.addLink(
                            origNodeID=way.nodes[startIdx].id,
                            destNodeID=node.id,
                            controlPoints=controlPoints,
                            metadata=way._asdict(),
                        )
                        linkCount += 1
                        if not way.oneWay:
                            # We need our map to be unidirectional, so create reverse link:
                            map.addLink(
                                origNodeID=node.id,
                                destNodeID=way.nodes[startIdx].id,
                                controlPoints=reversed(controlPoints),
                                metadata=way._asdict(),
                            )
                            reversedLinkCount += 1
                        startIdx = index
        logging.info(
            f"Number of nodes: {len(addedNodes)}; Number of links: {linkCount}; Reversed: {reversedLinkCount}"
        )

    def getChunk(self, lowCoords, highCoords):
        """
        Internal function to make a single Overpass API query for the specified bounding box
        """
        queryStr = f"[out:json];way({lowCoords[0]},{lowCoords[1]},{highCoords[0]},{highCoords[1]}){self.HIGHWAY_CLAUSE};(._;>;);out meta;"
        logging.info(
            f"Fetching from Overpass API ({lowCoords[0]:.3f}, {lowCoords[1]:.3f})-({highCoords[0]:.3f}, {highCoords[1]:.3f})"
        )
        queryStr = urllib.parse.quote(queryStr)
        response = requests.get(self.endpoint + "?data=" + queryStr)
        response.raise_for_status()
        result = response.json()

        # Pass #1: Nodes:
        nodeCount = 0
        for element in result["elements"]:
            if (
                "type" in element
                and element["type"] == "node"
                and element["id"] not in self.nodeCache
            ):
                sigFlag = False
                junctFlag = False
                if "tags" in element and "highway" in element["tags"]:
                    sigFlag = element["tags"]["highway"] == "traffic_signals"
                    junctFlag = element["tags"]["highway"] == "motorway_junction"
                self.nodeCache[element["id"]] = OSMReader.Node(
                    id=element["id"],
                    lat=element["lat"],
                    lon=element["lon"],
                    signal=sigFlag,
                    junction=junctFlag,
                )
                self.wayNodeLkp[element["id"]] = {*()}
                nodeCount += 1

        # Pass #2: Ways:
        wayCount = 0
        for element in result["elements"]:
            if "type" in element and element["type"] == "way":
                if "tags" in element and "highway" in element["tags"]:
                    ourName = (
                        element["tags"]["name"].strip()
                        if "name" in element["tags"]
                        else "none"
                    )
                    motorway = element["tags"]["highway"] == "motorway"
                    oneWay = motorway or (
                        "oneway" in element["tags"]
                        and (
                            element["tags"]["oneway"] == "yes"
                            or element["tags"]["oneway"] == "1"
                        )
                    )
                    assert (
                        not "oneway" in element["tags"]
                        or str(element["tags"]["oneway"]) != "-1"
                    ), f"Way {element['id']}: reverse one-way not supported"
                    nodes = tuple(self.nodeCache[nodeID] for nodeID in element["nodes"])
                    way = OSMReader.Way(
                        id=element["id"],
                        type=element["tags"]["highway"],
                        name=ourName,
                        oneWay=oneWay,
                        motorway=motorway,
                        nodes=nodes,
                    )
                    self.waySet |= {way}
                    for node in nodes:
                        self.wayNodeLkp[node.id] |= {way}
                    wayCount += 1
        logging.info(f"New nodes: {nodeCount}; New ways: {wayCount}.")
        return nodeCount
