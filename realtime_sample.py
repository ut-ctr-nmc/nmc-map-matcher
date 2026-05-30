"""
realtime_sample.py demonstrates the capability of map-matching a stream of
    trackpoints as they arrive in real time through a small simulation.
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
from support import mpo_read
from typing import Final, Iterator
import csv, time, os, sys, logging

# Configure logging to use stdout
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    stream=sys.stdout,
)

TRACK_FILE: Final[str] = os.path.join("samples", "small_tracks.csv")
MPO_PATH: Final[str] = os.path.join("samples", "mpo")

# Grab tiny MPO map:
map: graph.Map = mpo_read.readMPOModel(
    os.path.join(MPO_PATH, "small_atx_nodes.csv"),
    os.path.join(MPO_PATH, "small_atx_links.csv"),
    os.path.join(MPO_PATH, "small_atx_cnx.csv"),
)

# We'll need to keep track of this in the generator, and afterwards.
stableIndex: int = -1
currentIndex: int = -1


def trackpointStreamer(
    filename: str, pathEngine: path_engine.PathEngine
) -> Iterator[graph.Trackpoint]:
    """
    Simulates a stream of trackpoints by reading them from a file one by one.
    Assumes the file is formatted as sequence,lat,lon (with no header).
    """
    global stableIndex, currentIndex
    streetName = ("", "")
    logging.info("Begin streaming trackpoints.")
    with open(filename, "r") as inFile:
        csvReader = csv.reader(inFile)
        for line in csvReader:
            # Simulate the passage of time for the trackpoint to arrive:
            time.sleep(1)
            # TODO: Change over to asyncio for more modern streaming support
            # It is also possible to use stream blocking, say, on STDIN or a socket.

            # Bring this new trackpoint into where we are processing:
            logging.info(f"Received trackpoint: {line}")
            trackpoint: graph.Trackpoint = map.makeTrackpoint(
                lonHoriz=float(line[2]), latVert=float(line[1]), seq=int(line[0])
            )
            yield trackpoint

            # Now that we fed in a trackpoint, let's see if we can make sense of the track
            # so far. First, drill down through the hypotheses until we find a parent layer
            # that has just one path. Assumption then is that all others have been pruned.
            # Alternatively, it is possible to look at the score for each path (.totalCost)
            # and make decisions based on that.
            currentIndex += 1
            index = currentIndex
            parentsList: list[set[path_engine.PathEnd | None]] = [
                {*pathEngine.pathPointsPrev}
            ]
            while index > stableIndex:
                parentsList.append(
                    {
                        parent.prevTreeNode
                        for parent in parentsList[-1]
                        if parent is not None
                    }
                )
                index -= 1
            parentsList.reverse()  # We want to start with the oldest layer and work forward.
            while parentsList:
                if len(parentsList[0]) == 1:
                    # A single parent means we have a single lowest-cost path, and we are
                    # early enough in the tree that the algorithm has made this sole
                    # conclusion.
                    parent = next(iter(parentsList[0]))
                    if parent is not None and not parent.restart:
                        newStreetName = (
                            parent.pointOnLink.link.data["name"],
                            parent.pointOnLink.link.data["dir"],
                        )
                        if newStreetName != streetName:
                            logging.info(
                                f"+ Conclusion (@ {parent.totalDist:.1f} m): {newStreetName[0]} going {newStreetName[1]}"
                            )
                            streetName = newStreetName
                    del parentsList[0]
                    stableIndex += 1
                else:
                    break
    logging.info("Finished streaming trackpoints.")


# Run the map matcher on the trackpoint generator:
pathEngine = path_engine.PathEngine(
    # Use parameters that maintain fewer hypotheses so we can get conclusive results
    # sooner for this demo:
    path_engine.PathEngine.Params(limitClosestPoints=4, limitSimulPaths=3, maxHops=3)
)
finalList: list[path_engine.PathEnd] | None = pathEngine.constructPath(
    trackpointStreamer(TRACK_FILE, pathEngine), map
)

# The final list given by constuctPath() is the lowest-cost path through the tree. We'll
# report just on that, starting at the point we reported earlier while streaming.
index = stableIndex
streetName = ("", "")
while index <= currentIndex and finalList is not None:
    if finalList[index] is not None and not finalList[index].restart:
        newStreetName = (
            finalList[index].pointOnLink.link.data["name"],
            finalList[index].pointOnLink.link.data["dir"],
        )
        if newStreetName != streetName:
            logging.info(f"+ Most likely (@ {finalList[index].totalDist:.1f} m): {newStreetName[0]} going {newStreetName[1]}")
            streetName = newStreetName
    index += 1