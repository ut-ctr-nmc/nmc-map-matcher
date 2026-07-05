"""
reporter.py contains utility functions for assembling together tracks
    from match results.
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

from typing import Hashable, NamedTuple, Iterable
from dataclasses import dataclass

import shapely
from nmc_mm_lib import graph, path_engine


"Identifies minimum length of a span of intermediate points"
MINIMUM_SPAN: float = 1e-4


@dataclass(frozen=True)
class OutputTrackpoint(graph.Trackpoint):
    """
    A Trackpoint with additional fields for outputting track information.
    """

    linkID: Hashable  # The link ID of the corresponding point
    linkDist: (
        float  # Distance along the corresponding link that this point corresponds to
    )
    totalDist: float  # Cumulative distance along corresponding track
    distanceAway: float  # Distance of match away from topology. Measure of quality
    linksTrav: (
        list[Hashable] | None
    )  # A list of link IDs traversed to get to this point, None if restarted
    subseqFlag: bool  # Whether this point is part of a subsequence


def prepareTrackpath(
    map: graph.Map,
    treeNodes: Iterable[path_engine.PathEnd],
    intermediary: bool = False,
    increment: float | None = None,
    fractionalSeq: bool = True,
) -> list[OutputTrackpoint]:
    """
    Takes the results from the path engine and prepares a list of Trackpoint objects
    expressed in untransformed coordinates

    @param map: The map used in the path engine, needed to convert points back to
    untransformed coordinates
    @param treeNodes: The path engine results to be processed
    @param intermediary: Also outputs lon/lat of link starts; uses sub-sequences
    @param increment: Put a point once every given meters; uses sub-sequences
    @param fractionalSeq: Whether to use fractional sequence numbers
    """
    # TODO: Add abilities to space points equally between tree nodes, and to put
    # points at underlying geometry bends

    class StartsRecord(NamedTuple):
        link: graph.Map.LinkRecord
        dist: float

    ret: list[OutputTrackpoint] = []
    linksTrav: list[Hashable] = []
    treeNode: path_engine.PathEnd
    priorTreeNode: path_engine.PathEnd | None = None
    for treeNode in treeNodes:
        # Special points considerations: periodic increments and/or points at link starts:
        if (
            (intermediary or increment and increment > 0)
            and priorTreeNode
            and not treeNode.restart
            and treeNode.refPoint.seq is not None
        ):
            dist = treeNode.totalDist
            starts: list[StartsRecord] = [
                StartsRecord(link=treeNode.pointOnLink.link, dist=dist)
            ]
            dist = (
                dist
                - treeNode.pointOnLink.getDistanceAlong()
                + treeNode.pointOnLink.link.getLength()
            )
            for routeTraverse in reversed(treeNode.routeInfo):
                dist -= routeTraverse.getLength()
                starts.append(StartsRecord(link=routeTraverse, dist=dist))
            starts.append(
                StartsRecord(
                    link=priorTreeNode.pointOnLink.link,
                    dist=priorTreeNode.totalDist,
                )
            )
            startDist = starts[-1].dist
            span = starts[0].dist - startDist
            offset = priorTreeNode.pointOnLink.getDistanceAlong()

            if span >= MINIMUM_SPAN:
                curStart: StartsRecord = starts.pop()
                dist = curStart.dist
                while len(starts) >= 1:
                    popFlag = False
                    if increment:
                        dist += increment
                        offset += increment
                        if dist >= starts[-1].dist:
                            if not intermediary:
                                # We carry on in the path of the next link, but
                                # not necessarily starting at that link:
                                offset = dist - starts[-1].dist
                                # Retract one step so we end up in same spot
                                # when popped:
                                dist -= increment
                                offset -= increment
                                curStart = starts.pop()
                                linksTrav.append(curStart.link.id)
                                continue
                            # Otherwise, We ensure we mark the beginning of each link:
                            popFlag = True
                    else:
                        popFlag = True
                    if popFlag:
                        dist = starts[-1].dist
                        offset = 0
                        curStart = starts.pop()
                        if not starts:
                            break
                        linksTrav.append(curStart.link.id)

                    lon, lat = map.revertPoint(
                        *curStart.link.getPointAlong(offset, normalize=False)
                    )
                    if fractionalSeq:
                        seqValue = treeNode.refPoint.seq - 1 + (dist - startDist) / span
                    else:
                        seqValue = len(ret)
                    newRec = OutputTrackpoint(
                        lonHoriz=lon,
                        latVert=lat,
                        point=shapely.geometry.Point(lon, lat),
                        id=treeNode.refPoint.id,
                        seq=seqValue,
                        linkID=curStart.link.id,
                        linkDist=offset,
                        totalDist=dist,
                        distanceAway=0.0,
                        linksTrav=linksTrav.copy(),
                        subseqFlag=True,
                    )
                    ret.append(newRec)
                    linksTrav.clear()
        else:
            if not treeNode.restart:
                linksTrav.extend(
                    routeTraverse.id for routeTraverse in treeNode.routeInfo
                )

        lon, lat = map.revertPointOnLink(treeNode.pointOnLink)
        seqValue = len(ret) if not fractionalSeq else treeNode.refPoint.seq
        newRec = OutputTrackpoint(
            lonHoriz=lon,
            latVert=lat,
            point=shapely.geometry.Point(lon, lat),
            id=treeNode.refPoint.id,
            seq=seqValue,
            linkID=treeNode.pointOnLink.link.id,
            linkDist=treeNode.pointOnLink.getDistanceAlong(),
            totalDist=treeNode.totalDist,
            distanceAway=treeNode.pointOnLink.refDist,
            linksTrav=linksTrav.copy() if not treeNode.restart else None,
            subseqFlag=False,
        )
        ret.append(newRec)
        linksTrav.clear()
        priorTreeNode = treeNode
    return ret
