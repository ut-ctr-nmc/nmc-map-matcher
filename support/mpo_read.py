"""
mpo_read.py contains functions for reading the sample MPO CSV files
and building up a node-link representation of a roadway network.
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

from shapely import from_wkt
from nmc_mm_lib import graph
from typing import Any, Hashable, Generator
import csv


def mpoRead(filename: str) -> Generator[dict[str, Any]]:
    """
    Reads MPO CSV file and yields line by line
    """
    with open(filename, mode="r", newline="") as fileHandle:
        csvReader = csv.DictReader(fileHandle)
        for fileLine in csvReader:
            yield fileLine


MPOCollection = dict[Hashable, dict[str, Any]]

def readMPOModel(nodesFile: str, linksFile: str, cnxFile: str) -> graph.Map:
    """
    Reads MPO CSV files and builds up a graph.Map representation of the
    roadway network.
    """
    # Grab MPO model: we'll derive the topography from node ids (nodes.csv),
    # link segments (links.csv), and connectivity (cnx.csv). Nodes and links
    # have extra naming information that we'll put in as metadata.
    # TODO: This is much more compact with Pandas!
    nodes: MPOCollection = {}
    for fileLine in mpoRead(nodesFile):
        nodes[fileLine["id"]] = {
            "id": fileLine["id"],
            "lon": float(fileLine["lon"]),
            "lat": float(fileLine["lat"]),
            "name": fileLine["name"],
        }
    links: MPOCollection = {}
    for fileLine in mpoRead(linksFile):
        links[fileLine["id"]] = {
            "geog": from_wkt(fileLine["wkt"]),
            "name": fileLine["name"],
            "dir": fileLine["dir"],
        }
    cnxs: MPOCollection = {}
    for fileLine in mpoRead(cnxFile):
        cnxs[fileLine["id"]] = {
            "source": nodes[fileLine["source"]],
            "dest": nodes[fileLine["dest"]],
        }

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
    # Commit our geometry:
    map.completeMap()
    return map
