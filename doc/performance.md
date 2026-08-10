# Performance <!-- omit in toc -->

- [Where the Time Goes](#where-the-time-goes)
- [Caching Invariant Map Properties](#caching-invariant-map-properties)
  - [Reverse Links](#reverse-links)
  - [Link Length and Origin](#link-length-and-origin)
  - [Outgoing Link Adjacency](#outgoing-link-adjacency)
  - [Bulk Extraction from Shapely](#bulk-extraction-from-shapely)
- [Priority Queues](#priority-queues)
- [Target Link Accumulation](#target-link-accumulation)
- [Per-Instance State](#per-instance-state)
- [Results](#results)
- [Verifying Changes: Run-to-Run Nondeterminism](#verifying-changes-run-to-run-nondeterminism)
- [Remaining Opportunities](#remaining-opportunities)

This document records the profiling of the map matcher that was performed against `avl_osm_sample.py`, the optimizations that followed, and the measurement methodology needed to confirm that those optimizations left results unchanged.

The workload throughout is the AVL sample: 283 trackpoints from one CapMetro bus trip matched against an OpenStreetMap network of the greater Austin, TX metro area, comprising 35,603 nodes and 77,967 directed links. Timings come from a run under Python 3.14 with shapely 2.1.2, networkx 3.6.1, and pyproj 3.7.2.

## Where the Time Goes

Two facts emerged from profiling that shaped everything else.

First, cost is concentrated in a small minority of trackpoints. Sorting trackpoints by the time `PathEngine._findShortestPaths()` spends on each:

| Trackpoints | Share of match time |
| --- | --- |
| Worst 10 (3.5%) | 22% |
| Worst 28 (10%) | 43% |

The worst single trackpoint required 7,728 link expansions in `WalkPathProcessor._walkPath()` and took 0.97 s on its own. These spikes occur in dense downtown grids, where a high branching factor combined with `maxHops=16` lets the search enumerate a very large number of candidate paths.

Second, and more actionable, the bulk of the time was not spent on the search itself but on **repeatedly recomputing properties of a base map that never changes**. Ranked by cumulative time within the 22.9 s of matching:

| Call site | Calls | Cumulative |
| --- | --- | --- |
| `Map.isReverseLink()` | 284,964 | 11.1 s |
| `createNext()` → `geometry.length` | 210,636 | 4.7 s |
| `PointOnLink.findDistanceFrom()` | 219,931 | 4.0 s |
| `LinkRecord.getOrigin()` (`shapely.get_point`) | 172,785 | 3.4 s |
| `Map.outgoingLinks()` (NetworkX) | 194,300 | 2.6 s |
| `LinkRecord.getLength()` | 168,144 | 2.1 s |
| `queue.PriorityQueue` `put`/`get` | 452,792 | 3.1 s |

`isReverseLink()` alone accounted for roughly 28% of match time. Every one of these entries is a per-call cost paid on a value that could have been computed once when the map was built.

This also explains an earlier negative result: introducing priority queues to order the pathfinding work did not speed things up. The ordering was sound, but each queue *element* carried microseconds of GEOS and NetworkX round-trips, so reordering the work could not reduce the work.

## Caching Invariant Map Properties

`Map.completeMap()` already existed as the point where a base map is finalized and its spatial index built. It now also precomputes the per-link quantities the pathfinder needs. Nothing about the search algorithm changed.

### Reverse Links

`isReverseLink()` determines whether one link retraces another backwards, which `_walkPath()` needs in order to assess U-turn penalties. Establishing this requires comparing rounded lengths and then walking both coordinate sequences — but for a fixed base map the answer never changes.

`completeMap()` now resolves the relation once into `Map.reverseLinkLookup`, a mapping of link ID to the frozenset of link IDs that reverse it. Only links joining the same pair of nodes in opposite directions can possibly qualify, so the comparison runs about once per opposing pair rather than hundreds of thousands of times during matching. The innermost loop calls the new `Map.isReverseLinkCached()`, which is a set membership test.

The comparison rule itself was factored into `Map._reverseMatches()` so that `completeMap()` and the public `isReverseLink()` share one implementation; `completeMap()` passes it coordinates it has already extracted.

### Link Length and Origin

`LinkRecord.getLength()` and `LinkRecord.getOrigin()` previously called into GEOS on every invocation. Both now read values cached into the edge data (`length`, `origin`, `originXY`) by `completeMap()`.

In `_walkPath()`, the length of the origin link is invariant across an entire `walkPath()` call, so it is hoisted into `WalkPathProcessor.origLinkLength` once per search.

The distance computation was deliberately left as `shapely` `Point.distance()`, reached through `PointOnLink.findDistanceFrom()` and the cached origin `Point`. Substituting planar arithmetic would have been faster still, but keeping every distance behind one call site preserves the code's dependence on the underlying projection scheme, which may not be planar. Note that shapely's `distance` is itself Cartesian; the value of the abstraction is that a future geodesic implementation would need to change only one place.

### Outgoing Link Adjacency

`Map.outgoingLinks()` queried `networkx.MultiDiGraph.edges()` on every expansion. `completeMap()` now builds `Map.outLinkLookup`, a plain dictionary from node ID to a tuple of outgoing `LinkRecord`s, and `outgoingLinks()` is a dictionary lookup. Its return type changed from a generator to a tuple, which also lets `_walkPath()` test for a dead end with `len(outgoing) == 1` instead of consuming an iterator.

One subtlety: `_walkPath()` may replace the list of links to consider with a single shortcut from `GlobalPathCache`. The dead-end test must consult the node's true out-degree, not the possibly-substituted list, so `outgoing` and `myList` are tracked separately.

### Bulk Extraction from Shapely

The first implementation of the above made map building about 9 s *slower* than it saved, because reading `.length`, `.coords`, and `get_point()` costs on the order of 34 µs per geometry and there are 77,967 of them.

`completeMap()` therefore extracts these in bulk using shapely's vectorized entry points — `shapely.length()`, `shapely.get_point()`, `shapely.get_coordinates()`, and `shapely.get_num_coordinates()` — each called once over the whole list of geometries. Per-link coordinate sequences are then sliced out of the flat coordinate array using the per-geometry counts. This reduced `completeMap()` from 11.4 s to roughly 2 s.

The net effect on map building is about +1.3 s over the original, paid once, in exchange for 17.8 s saved during matching.

## Priority Queues

Both `WalkPathProcessor.walkPath()` and `PathEngine._findShortestPaths()` used `queue.PriorityQueue`, which acquires a lock on every `put()`, `get()`, `empty()`, and `qsize()` call. The matcher is single-threaded, so this synchronization is pure overhead. Both now use `heapq` over a plain list.

The tie-breaking counters were preserved exactly: `_findShortestPaths()` previously used `pairQueue.qsize()` as its monotonic insertion counter, which is equivalent to `len(pairQueue)` before a push, since every push happens before any pop.

## Target Link Accumulation

`WalkPathProcessor.Params.allTargetLinkIDs` gates whether `_walkPath()` writes a discovered path into `GlobalPathCache`. `PathEngine.constructPath()` only ever called `update()` on it, never clearing it, so the set grew monotonically across an entire track. The cache-writing branch therefore fired on an ever-larger set of links as a run progressed — links that had already been discarded by the `GlobalPathCache.trimExcept()` call immediately preceding, so the work was wasted and the cache polluted.

The set is now replaced per trackpoint, matching the trim that runs alongside it.

## Per-Instance State

Several mutable containers were declared as class attributes and so were shared by every instance in the process:

- `Map.graph`, `Map.linkIDLookup`
- `GlobalPathCache.pathCache`
- `PathEngine.globalPathCache`, `PathEngine.PrevCosts.costs`
- the `set()` default of `WalkPathProcessor.Params.allTargetLinkIDs`, which is a NamedTuple field default and therefore also shared

Constructing two `Map` or two `PathEngine` objects in one process would have silently entangled their state. All are now initialized per instance. This has no effect on the single-instance sample codes, but it removes a hazard for anything that matches against more than one base map.

## Results

Measured with `PYTHONHASHSEED=0` (see below), on the AVL sample:

| | Before | After | Speedup |
| --- | --- | --- | --- |
| Matching (`constructPath`) | 22.9 s | 5.1 s | 4.5× |
| Worst single trackpoint | 0.97 s | 0.17 s | 5.7× |
| Whole sample, including map build | 33.1 s | 16.4 s | 2.0× |
| Link expansions in `_walkPath()` | 209,422 | 209,422 | unchanged |

The expansion count being identical is the confirmation that matters: the search explores exactly the same territory as before, each step simply costs less. Output was byte-identical for both `avl_osm_sample.py` and `gtfs_small_sample.py`.

## Verifying Changes: Run-to-Run Nondeterminism

**The matcher produces different results on different runs of identical code.** Unmodified code was observed to yield three distinct outputs across four consecutive runs of `avl_osm_sample.py`.

The cause is in `support/osm_overpass.py`. `OSMReader.addToMap()` iterates `self.waySet`, a `set` whose members are `OSMReader.Way` NamedTuples containing a string `name` field. Python randomizes string hashing per process, so the set's iteration order differs every run. That determines the order links are added to the graph, which determines `treeIndex`, which determines the order `shapely.strtree.STRtree` returns query results, which decides ties in the `foundLinks.sort(key=refDist)` at the end of `Map.findPointsOnLinks()`. A different candidate ordering can produce a different matched path.

The practical consequence is that **comparing matcher output before and after a change is meaningless unless the hash seed is pinned.** During this work an optimization was initially recorded as byte-identical when the match was coincidence; the same code produced a different file on the next run.

To compare two revisions:

```
PYTHONHASHSEED=0 python3 avl_osm_sample.py
```

on each side, then compare `avl_matched.csv`. Sorting `waySet` by `way.id` before iterating would make runs reproducible without a pinned seed, if that is wanted.

## Remaining Opportunities

The optimizations above removed per-step overhead but did not change the shape of the search, so the concentration of cost remains: after the changes, 10% of trackpoints still account for 38% of match time.

That tail comes from `WalkPathProcessor.Next.backtrackSet`, which forbids revisiting a link within a path and so makes the search enumerate *simple paths* rather than compute *shortest paths*. In a dense grid at `maxHops=16` that count grows combinatorially. Three changes target it:

1. **Label each link once.** All costs are non-negative and the U-turn penalty depends only on the (previous link, link) pair, so a link-labeled search — effectively Dijkstra on the line graph — is exact and bounds work by the number of links within `limitPathDist` rather than the number of paths within `maxHops`.
2. **True A\* with early exit.** `PathElement` currently sorts by `destDistance` before `cost`, and `destDistance` holds the *parent's* crow-flies distance to the goal, making the search greedy best-first; `walkPath()` then drains the entire queue looking for something better. Keying on `cost + distFactor × crowDistance` is admissible, which permits stopping at the first pop of the destination link.
3. **One multi-target search per origin.** `_findShortestPaths()` runs up to `limitSimulPaths × limitClosestPoints` (8 × 12 = 96) separate searches per trackpoint over nearly identical territory. Since the candidates are clustered around a single trackpoint, one search per origin that terminates when all candidates are settled would cut this by roughly an order of magnitude. This is the intent of the existing `TODO` in that method.

These change which path wins in ambiguous cases, so they should be evaluated against match quality and not only against wall-clock time.
