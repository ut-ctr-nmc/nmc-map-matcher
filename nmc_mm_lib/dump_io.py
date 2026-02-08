"""
dump_io.py contains code for inputting and outputting trackpoints that
    are matched to base map paths.
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

from nmc_mm_lib import graph, path_engine
from collections.abc import Hashable, Iterable
from typing import IO, Mapping, TypedDict
from numbers import Number
import csv, json, sys
import logging


class StdFieldNames(TypedDict):
    trackID: Hashable
    trackSeq: int
    linkID: Hashable
    linkDist: float
    totalDist: float
    lon: float
    lat: float
    numLinksTrav: int
    linksTrav: str


def dumpStandardInfo(
    treeNodesLists: Mapping[Hashable, Iterable[path_engine.PathEnd]],
    outFile: IO = sys.stdout,
    includeHeader: bool = True,
) -> None:
    """
    Outputs the body of a CSV format of track path information.
    """
    writer = csv.DictWriter(outFile, fieldnames=StdFieldNames.__annotations__.keys())
    if includeHeader:
        writer.writeheader()
    treeNodes: Iterable[path_engine.PathEnd]
    for treeNodes in treeNodesLists.values():
        treeNode: path_engine.PathEnd
        for treeNode in treeNodes:
            outData: StdFieldNames
            outData = {
                "trackID": treeNode.refPoint.id,
                "trackSeq": (
                    treeNode.refPoint.seq if treeNode.refPoint.seq is not None else -1
                ),
                "linkID": treeNode.pointOnLink.link.id,
                "linkDist": treeNode.pointOnLink.getDistanceAlong(),
                "totalDist": treeNode.totalDist,
                "lon": treeNode.pointOnLink.point.x,
                "lat": treeNode.pointOnLink.point.y,
                "numLinksTrav": len(treeNode.routeInfo) if not treeNode.restart else -1,
                "linksTrav": str(
                    [routeTraverse.id for routeTraverse in treeNode.routeInfo]
                    if not treeNode.restart
                    else []
                ),
            }
            # A links traversed length of -1 shall be a special indication
            # saying that we are restarting, and the link list doesn't exist.
            writer.writerow(outData)


def stringToList(string: str) -> list[Number]:
    """
    stringToList is a helper function to safely convert a string representation
    of a list of numbers back into a list.
    """
    try:
        myList = json.loads(string)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid list format: {e}")
    if not isinstance(myList, list) or not all(isinstance(x, Number) for x in myList):
        raise ValueError("Input string must be a list of numbers.")
    return myList


def readStandardDump(
    baseMap: graph.Map, inFile: IO, includeHeader: bool = True
) -> dict[Hashable, list[path_engine.PathEnd]]:
    """
    readStandardDump reconstructs the tree entries that PathEngine had created.

    @return A dictionary of trackID to a list of PathEnds
    """
    ret: dict[Hashable, list[path_engine.PathEnd]] = {}
    params = {}
    if not includeHeader:
        params["fieldnames"] = list(StdFieldNames.__annotations__.keys())
    reader = csv.DictReader(inFile, **params)

    for inData in reader:
        trackPoint = baseMap.makeTrackpoint(
            float(inData["lon"]),
            float(inData["lat"]),
            ident=inData["trackID"],
            seq=int(inData["trackSeq"]),
        )
        link = baseMap.getLinkByID(inData["linkID"])
        if not link:
            logging.warning(
                "The path match file refers to a nonexistent link"
                + f" ID {inData['linkID']}."
            )
            continue
        refDist, percentAlong, isPerpendicular, pointAlong = baseMap.pointDist(
            trackPoint, link
        )
        matchPoint = graph.Map.PointOnLink(
            link, percentAlong, not isPerpendicular, refDist, pointAlong
        )
        newEntry = path_engine.PathEnd(trackPoint, matchPoint)
        newEntry.totalDist = float(inData["totalDist"])
        if int(inData["numLinksTrav"]) >= 0:
            newEntry.routeInfo = [
                baseMap.getLinkByID(linkID)
                for linkID in stringToList(inData["linksTrav"])
            ]
        newEntry.restart = int(inData["numLinksTrav"]) == -1
        # TODO: totalCost and totalLinkCount aren't being stored in the dump,
        # so they're not set.

        if inData["trackID"] not in ret:
            ret[inData["trackID"]] = []
        ret[inData["trackID"]].append(newEntry)

    # Now we need to make sure sequences are sorted:
    for treeNodes in ret.values():
        treeNodes.sort(
            key=lambda node: node.refPoint.seq if node.refPoint.seq is not None else -1
        )
    return ret
