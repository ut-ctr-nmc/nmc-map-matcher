# Path Match <!-- omit in toc -->

- [PathEngine](#pathengine)
  - [Restarts and Errors](#restarts-and-errors)
- [Reporting](#reporting)

## PathEngine

Path matching functionality is implemented in `nmc_mm_lib.path_engine`, class `PathEngine`. In short, given an already completed base map and list of trackpoints, to start:

```python
engine = nmc_mm_lib.path_engine.PathEngine(params: PathEngine.Params | None)
results: list[path_engine.PathEnd] = engine.constructPath(trackpoints: Iterable[graph.Trackpoint], baseMap: graph.Map)
```

In this, underlying map topology is found that allows one to traverse from one trackpoint to the next.

When instantiated, `PathEngine` takes a `PathEngine.Params` NamedTuple that has the following members. Note that all units here are expressed in meters *only if* the *working CRS* is set to a projection that provides projections in meters, and that the projections are nearly Euclidian. The default `EPSG:3857` working CRS supposedly does this on a worldwide scale, although accuracy should improve if a more localized CRS is used. For example, `EPSG:3082` is what had been used for related projects in Texas.

Parameters for instantiating `PathEngine`:

* **searchRadius:** Radius (m) to search from each trackpoint to perpendicular underlying map links. (Default: 100 m)
* **radiusPrimary:** Radius (m) to search from each trackpoint to underlying map link endpoints. (Default: 100 m). These first two settings correspond mostly to the expected accuracy of trackpoints with respect to the underlying map topology. For example, if there is significant GPS distortion caused by "urban canyons", these values would need to be large enough to cover that distortion. But, values too large can cause the algorithm to consider areas on the underlying map that are clearly irrelevant, resulting in longer processing time.
* **radiusSecondary:** Radius (m) to search from each trackpoint to a previous underlying map point. This is used in cases where a trackpoint may be far away, but the previously found underlying map point is nearby. This can happen if a bus travels into, say, a parking lot, where no corresponding links exist in the underlying map network. (Default: 50 m)
* **limitPathDist:** Maximum path distance (m) allowed in finding paths from a previous underlying map point to a newly proposed underlying map point.  (Default: 500 m). This affects how extensively the pathfinding algorithm can search in terms of traversed topology, from one trackpoint to the next.
* **limitDirectDist:** Maximum radius (m) from a previous underlying map point to a newly proposed map point allowed in pathfinding. (Default: 500 ft). This affects the range of traversed topology "as far as the crow flies" from one trackpoint to the next. A gap would definitely appear if the distance from one trackpoint to another exceeds this setting.
* **limitDirectDistRev:** Radius (m) to allow backtracking on an existing link during pathfinding. This can happen if several trackpoints "pile up" on a underlying map link when a bus travels into areas not represented in the underlying map network, such as a parking lot. (Default: 160 m)
* **distFactor:** Cost multiplier for linear path distance in terms of traversed topology (Default: 1.0). This assesses uniform penalty per meter traveled.
* **driftFactor:** Cost multiplier for distance from a trackpoint to its corresponding underlying map link. This is used to favor proposed underlying map points that are closer to trackpoints. (Default: 2.0)
* **nonPerpPenalty:** Penalty multiplier for trackpoints that aren't perpendicular to underlying map links. While endpoints of underlying map links may be considered as candidates for path-matching, underlying map links that are perpendicular to the trackpoints are favored. (Default: 1.5). The paper goes into detail about the significance of this.
* **limitClosestPoints:** Number of close-proximity underlying map points that are considered for each trackpoint. Those that are perpendicularly closest to existing underlying map points are prioritized first. An extremely densely-featured underlying map may need a higher value than a sparsely-detailed map. This directly affects the speed of pathfinding. (Default: 12 points)
* **limitSimulPaths:** Number of proposed paths to maintain during pathfinding stage. This is the maximum number of "hypotheses" or "candidate paths" that are maintained at any one time. This is prioritized according to the calculated match cost. This also directly affects the speed of pathfinding.  (Default: 8 paths)
* **maxHops:** Maximum number of underlying map links to traverse in a path-finding operation. A densely-featured map with many short links between nodes needs a higher number, so that all of the necessary links can be traversed from one trackpoint to the next. (Default: 12 hops). Speed of pathfinding is also affected by this value.
* **tossRatio:** Allows for the invalidation of short paths that are missing successful matches on one or both of their ends. This especially affects the validity of trackpoint paths that cover more area than the extentes of an underlying map. A path is called "invalid" if the number of unmatched trackpoints divided by the total number of trackpoints exceeds this value. (Default: 1.0, which disables this feature)

A couple other settings:

* In `nmc-mm-lib.graph.WalkPathProcessor.Params`: **uTurnInterPenalty:** Penalty to add to U-turns in intersections, or None to disable U-turns in intersections. (Default: None)
* In `nmc-mm-lib.graph.WalkPathProcessor.Params`: **uTurnDeadEndPenalty:** Penalty to add to U-turns at dead-ends, or None for `uTurnInterPenalty`. (Default: 50)
* In sample codes: **REPORT_STEP_SIZE:** Distance from one reported map-matched trackpoint to the next, generated equidistantly from "critical trackpoints" that match with the originally submitted trackpoints and intersections.

**NOTE** that if a CRS expressing distances in units other than *meters* (e.g. U.S. American *feet*) would require all default parameters assuming *meters* to be explicitly provided in terms of the alternative units. Again, the map matching algorithm assumes a nearly Euclidian ("approximately flat Earth") measurement scheme.

### Restarts and Errors

A "restart" (marked as `None` in a `path_engine.PathEnd` describing a matched series of points) marks a point where the map matcher was unsuccessful in finding a continuous path from the previous path sequence to the current one. This happens if:

* The trajectory goes off of the underlying topology enough for the map matcher to lose track of its location, and must pick up again or "restart" at another location when proximity to the underlying topology is regained. This happens, for example, when a trajectory goes through a parking lot that isn't represented in the underlying topology. There may even be the case where several matched points may "bunch up" on one end of the topology (with an increasing offset distance), and then suddenly appear at a different location when a "restart" happens.
* The vehicle made a U-turn along a roadway. This can happen when the map matcher is not configured to allow U-turns. Instead, the trajectory may proceed one direction, and then "restart" going back the other direction.
* The trajectory exits and re-enters the edge of the available underlying topology.

Note that when a "restart" happens, the total distance traveled is likely not cumulated accurately because the map matcher didn't find trajectory to run the route through. In this case, the direct distance is recorded across the problematic gap.

These are a couple key error messages that you may see:

* `WARNING: No map paths were found for path ##, sequence ##.` This means that no topology in the underlying map was found to coincide with the reported trackpoint. This may happen if the trackpoint exists well outside of a map's bounding rectangle, or if trackpoints are located in spaces not represented by the underlying map. For example, if the underlying map only represents major state highways, but the trackpoint exists in the middle of a large parking lot, and parameters aren't loose enough for a nearby highway to be found, this error will appear. Parameters to check are: `searchRadius`, `radiusPrimary`, and `radiusSecondary`.
* `WARNING: No closest links found for trackpoint ##, seq. ##.` This means that no paths were found through the underlying map that would allow connection from the previous trackpoint to the reported one. That doesn't necessarily mean that the underlying map *doesn't* have a connection; however, it may mean that the matching parameters are set too conservative. Parameters to check include: `limitPathDist`, `limitDirectDist`, `limitDirectDistRev`, and `maxHops`. Others beyond these may exhaust further searching by penalizing viable matches too aggressively.

## Reporting

Results are output from `path_engine.PathEngine.constructPath()` throgh a series of `path_engine.PathEnd` objects that represent the last remaining, winning series of points tracing back to the path's origin. The `reporter.prepareTrackpath()` method will go through such a list and create a series of `reporter.OutputTrackpoint` (extends `graph.TrackPoint`) that represent the path in several more useful ways. These include:

* One point snapped to underlying topology per input trackpoint, along with `.distanceAway` a measure of physical distance between them
* Optionally throwing in points representing the starts of each link that is traversed in the underlying topology (called "intermediary points"). This greatly helps in creating a list of driving directions, for example. For these, `.subseqFlag` is `True` and `.distanceAway` is always `0`.
* Optionally adding in points between matches at a predefined distance interval, which can be useful for later sampling of underlying topology at a finer resolution than the original input trackpoints

To help with outputting, a utility method for writing lists of output trackpoints to a CSV file is found in `dump_io.dumpStandardInfo()`.
