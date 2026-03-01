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

from nmc_mm_lib import graph, path_engine, reporter
from collections.abc import Hashable, Iterable
from typing import IO, Mapping, TypedDict
from numbers import Number
import csv
import json
import sys
import logging


class StdFieldNames(TypedDict):
    trackID: Hashable
    trackSeq: int | float
    linkID: Hashable
    linkDist: float
    totalDist: float
    lon: float
    lat: float
    numLinksTrav: int | None
    linksTrav: str | None


def dumpStandardInfo(
    trackpointLists: Mapping[Hashable, Iterable[reporter.OutputTrackpoint]],
    outFile: IO = sys.stdout,
    includeHeader: bool = True,
) -> None:
    """
    Outputs the body of a CSV format of track path information.
    """
    writer = csv.DictWriter(outFile, fieldnames=StdFieldNames.__annotations__.keys())
    if includeHeader:
        writer.writeheader()
    trackpointList: Iterable[reporter.OutputTrackpoint]
    for trackpointList in trackpointLists.values():
        linkIDInt = True
        trackpoint: reporter.OutputTrackpoint
        for trackpoint in trackpointList:
            if trackpoint.linksTrav is not None:
                for linkID in trackpoint.linksTrav:
                    try:
                        int(linkID)
                    except ValueError:
                        linkIDInt = False
                        break
            if not linkIDInt:
                break
        for trackpoint in trackpointList:
            linkTravStr = None
            if trackpoint.linksTrav is not None:
                if linkIDInt:
                    linkTravStr = str([int(linkID) for linkID in trackpoint.linksTrav])
                else:
                    linkTravStr = str(trackpoint.linksTrav)
            # TODO: Facilitate rounding:
            outData: StdFieldNames = {
                "trackID": trackpoint.id,
                "trackSeq": trackpoint.seq if trackpoint.seq is not None else -1,
                "linkID": trackpoint.linkID,
                "linkDist": trackpoint.linkDist,
                "totalDist": trackpoint.totalDist,
                "lon": trackpoint.lonHoriz,
                "lat": trackpoint.latVert,
                "numLinksTrav": (
                    len(trackpoint.linksTrav)
                    if trackpoint.linksTrav is not None
                    else -1
                ),
                "linksTrav": linkTravStr,
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
    # TODO: Will currently not read dumps that have intermediary points
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
            linkList: list[graph.Map.LinkRecord] = []
            for linkID in stringToList(inData["linksTrav"]):
                linkRecord = baseMap.getLinkByID(linkID)
                if not linkRecord:
                    logging.warning(
                        "The path match file refers to a nonexistent link"
                        + f" ID {linkID} in the linksTrav list."
                    )
                    continue
                linkList.append(linkRecord)
            newEntry.routeInfo = linkList
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
