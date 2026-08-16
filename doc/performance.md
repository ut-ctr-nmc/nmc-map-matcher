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
- [Reproducibility](#reproducibility)
  - [What Was Nondeterministic](#what-was-nondeterministic)
  - [Ordering the OSM Reader](#ordering-the-osm-reader)
  - [Verification](#verification)
- [Where the Time Goes Now](#where-the-time-goes-now)
- [Remaining Opportunities](#remaining-opportunities)
  - [1. Redundant Searching Within a Trackpoint](#1-redundant-searching-within-a-trackpoint)
  - [2. Simple-Path Enumeration](#2-simple-path-enumeration)
  - [3. Crow-Flies Distance](#3-crow-flies-distance)
  - [Leftovers from the Caching Pass](#leftovers-from-the-caching-pass)
  - [Concurrency](#concurrency)
  - [Evaluating These](#evaluating-these)
- [Measurement Notes](#measurement-notes)

This document records the profiling of the map matcher that was performed against `avl_osm_sample.py`, the optimizations that followed, the work that reproducible measurement depends on, and what is left to do.

`avl_osm_sample.py` is the most involved end-to-end exercise in the repository and is the workload throughout. It matches CapMetro bus trips against an OpenStreetMap network of the greater Austin, TX metro area, comprising 35,603 nodes and 77,967 directed links, fetched once from the Overpass API and thereafter read from `osm_cache.json`. The original profiling described below covered a single trip of 283 trackpoints; the sample now runs two trips totaling 674 trackpoints, which is what the later measurements reflect. Timings come from runs under Python 3.14 with shapely 2.1.2, networkx 3.6.1, and pyproj 3.7.2.

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

Measured with `PYTHONHASHSEED=0`, on the single-trip AVL workload, under `cProfile` (see [Measurement Notes](#measurement-notes)):

| | Before | After | Speedup |
| --- | --- | --- | --- |
| Matching (`constructPath`) | 22.9 s | 5.1 s | 4.5× |
| Worst single trackpoint | 0.97 s | 0.17 s | 5.7× |
| Whole sample, including map build | 33.1 s | 16.4 s | 2.0× |
| Link expansions in `_walkPath()` | 209,422 | 209,422 | unchanged |

The expansion count being identical is the confirmation that matters: the search explores exactly the same territory as before, each step simply costs less. Output was byte-identical for both `avl_osm_sample.py` and `gtfs_small_sample.py`.

## Reproducibility

Everything above depends on being able to tell whether a change altered the match. That was not previously possible, and fixing it came before any further optimization.

### What Was Nondeterministic

**The matcher used to produce different results on different runs of identical code.** Three consecutive runs of `avl_osm_sample.py`, unmodified, produced three distinct `avl_matched.csv` files.

The cause was in `support/osm_overpass.py`. `OSMReader.addToMap()` iterated `self.waySet`, a `set` whose members were `OSMReader.Way` NamedTuples containing a string `name` field. Python randomizes string hashing per process, so the set's iteration order differed every run. That determined the order links were added to the graph, which determined `treeIndex`, which determined the order `shapely.strtree.STRtree` returns query results, which decided ties in the `foundLinks.sort(key=refDist)` at the end of `Map.findPointsOnLinks()`. A different candidate ordering can produce a different matched path.

The practical consequence was that comparing matcher output before and after a change was meaningless unless the hash seed was pinned. During the optimization work an improvement was initially recorded as byte-identical when the match was coincidence; the same code produced a different file on the next run.

Two further orderings were unreliable even with a pinned seed, because they were decided outside the process:

- The order the Overpass API lists elements in a response. Nothing in the API guarantees it, and it governed which way records existed and in what order they were stored.
- The order of tag keys within a way. `Way.tags` was built by iterating `element["tags"]` as parsed from JSON, so two responses describing the same way with its tags in a different order produced two unequal `Way` records, both of which were kept.

### Ordering the OSM Reader

`OSMReader` now keys ways by their OSM ID instead of collecting them in a set:

- `waySet: set[Way]` became `wayLkp: dict[Hashable, Way]`, and `wayNodeLkp` maps a node ID to the set of *way IDs* passing through it rather than to `Way` records. Deduplication is now by identity rather than by whole-record equality, so a way that arrives twice from overlapping bounding-box chunks is stored once even if the two copies differ in some incidental respect. Previously both copies were kept, which inflated `wayNodeLkp` counts and could split a way at a node that is not really a junction.
- `addWay()` is the single place ways are recorded, shared by the cache reader and the Overpass reader.
- `sortedWays()` and `sortedNodes()` return records in ID order. `addToMap()` walks `sortedWays()`, so link insertion order — and therefore `treeIndex`, and therefore tie-breaking — is fixed by OSM IDs rather than by arrival order or hash order.
- `Way.tags` is built with `tuple(sorted(...))` in both the fetch path and the cache path, so a way is described identically no matter how the API ordered its tags.
- `cacheWrite()` writes nodes and ways through those same sorted accessors, so `osm_cache.json` is byte-identical for a given query. The file is diffable, and a run from cache matches the fetch it came from.

Note that the cache file's own order no longer influences anything, since ordering is imposed when the records are consumed rather than when they are read. An older, unsorted cache file therefore still yields the correct, deterministic result.

### Verification

- Three unpinned runs of `avl_osm_sample.py` plus one with `PYTHONHASHSEED=0` produced four byte-identical `avl_matched.csv` files. Before the change, three unpinned runs produced three different files.
- `cacheWrite()` run in three processes under different hash seeds produced byte-identical 29,290,100-byte cache files.
- Matching from a reordered cache file produced output identical to matching from the original, confirming that cache order is no longer an input to the result.
- A synthetic Overpass response was reconstructed from the cached data and replayed four times with elements and tag keys shuffled differently each time, with no network access. All four produced the same 28,108 way records and the same link insertion order.
- Link expansions moved from 401,070 to 402,050, a 0.2% difference that reflects ties now breaking a different way. The change is performance-neutral.

Because link ordering changed, output from this revision does not match output from before it. That is expected and unavoidable: any fixed order differs from whatever order a particular earlier run happened to draw.

## Where the Time Goes Now

Current wall-clock for the whole of `avl_osm_sample.py`, both trips, reading the OSM data from cache: **5.5 s**.

| Stage | Time |
| --- | --- |
| `geoRead()` — parse 29 MB cache | 0.8 s |
| `addToMap()` — build nodes and links | 1.9 s |
| `completeMap()` — precompute and index | 1.1 s |
| Matching, 674 trackpoints | 2.0 s |
| `prepareTrackpath()` and CSV output | 0.1 s |

Map building is now the larger half of the run and is dominated by constructing 78,041 shapely `LineString`s. Within matching, 402,050 link expansions are performed across 34,761 `walkPath()` calls, and the cost concentration the earlier work noted persists: the worst 10% of trackpoints account for 32% of match time, and the worst single trackpoint still drives 7,728 expansions.

## Remaining Opportunities

The optimizations above removed per-step overhead but did not change the shape of the search. Instrumenting `walkPath()` over both AVL trips shows where the remaining work actually goes:

| | |
| --- | --- |
| `walkPath()` calls | 34,761 (51.6 per trackpoint; the worst hit the full 8 × 12 = 96) |
| Calls that return a path | 10,178 |
| Calls rejected cheaply, before searching, by the `limitDirectDist` gate | 4,508 |
| **Calls that search the graph and return nothing** | **20,075** |
| Expansions spent on those fruitless calls | 323,015 of 402,050 (**80%**) |
| Time spent on those fruitless calls | 1.49 s of 2.08 s (**72%**) |

The dominant cost is not finding the winning path. It is re-exploring the same neighborhood up to 96 times per trackpoint and coming up empty in two thirds of those attempts. The three items below are ordered by what that measurement implies, which is not the order they were first written down in.

### 1. Redundant Searching Within a Trackpoint

`_findShortestPaths()` runs up to `limitSimulPaths × limitClosestPoints` separate searches per trackpoint over nearly identical territory. Since the candidates are clustered around a single trackpoint, one search per origin that terminates when all candidates are settled would replace roughly 96 searches with 8 — and it attacks the 80% directly, because a search that fails to reach any candidate currently pays the full cost of exploring the reachable region once per candidate. This is the intent of the existing `TODO` in that method.

### 2. Simple-Path Enumeration

`WalkPathProcessor.Next.backtrackSet` forbids revisiting a link within a path, which makes the search enumerate *simple paths* rather than compute *shortest paths*. In a dense grid at `maxHops=16` that count grows combinatorially; expansions (402,050) and queue pushes (367,289) being nearly equal confirms the queue is essentially drained rather than pruned.

All costs are non-negative and the U-turn penalty depends only on the (previous link, link) pair, so a link-labeled search — effectively Dijkstra on the line graph — is exact and bounds work by the number of links within `limitPathDist` rather than the number of paths within `maxHops`.

A related change is to make the queue ordering admissible. `PathElement` sorts by `destDistance` before `cost`, and `destDistance` holds the *parent's* crow-flies distance to the goal, so the search is greedy best-first rather than A\*; `walkPath()` then drains the queue looking for something better. Keying on `cost + distFactor × crowDistance` would permit stopping at the first pop of the destination link. Note that early exit *by itself* is worth less than it appears: only 4% of expansions currently happen after a winner has been recorded. The value is in the ordering doing the pruning, not in the stopping.

### 3. Crow-Flies Distance

`shapely`'s `distance` is now the largest single library cost in matching: 439,311 calls, 1.51 s of 4.56 s profiled (33%), plus the numpy scalar unwrapping behind each one.

The earlier work deliberately kept this behind `PointOnLink.findDistanceFrom()` and the cached origin `Point` so that the code's dependence on the projection scheme stays at one call site, in case a future implementation is geodesic rather than planar. That reasoning still holds, but the price is now measurable, and shapely's `distance` is Cartesian regardless. Caching each point's `(x, y)` and using `math.hypot()` inside that one method would recover most of it without spreading the flat-earth assumption any further than it already goes.

### Leftovers from the Caching Pass

Small items that the caching work did not reach:

- `Map.findPointsOnLinks()` still calls `self.graph.out_degree()` — the last NetworkX call on a matching hot path — when `outLinkLookup` already holds the answer.
- `Map.PointOnLink.getDistanceAlong()` still reads `geometry.length` from GEOS instead of the cached `data["length"]`.
- `graph.hasMoreThan()` and `graph.hasExactly()` are now dead. They existed for the generator-based `outgoingLinks()`, which returns a tuple.
- `WalkPathProcessor.Params.limitDirectDistRev` is plumbed through from `PathEngine.Params` and documented in [Path Match](path_match.md), but nothing reads it. Either the backtracking limit it describes was lost, or the parameter should be removed.
- `GlobalPathCache` stores no cost alongside a cached next link, so `_walkPath()` substitutes a cached shortcut for the full outgoing list without being able to check whether it is still the cheaper option. Adding the "expense" third element noted in the `TODO` there is what would make that substitution safe to trust.

### Concurrency

`README.md` advertises that the algorithm "may be parallelized *(To be implemented)*", and accelerating with concurrency is item 2 of the roadmap in [Home](home.md). Nothing runs concurrently today; there is no use of `threading`, `multiprocessing`, or `concurrent.futures` anywhere in the project.

The structure anticipates it — `_findShortestPaths()` bundles a `WalkPathProcessor` into each `pairQueue` element, and `path_engine.py` has a branch commented as handling the case where "during multithreading better scores were logged" — but the work is sequential. This overlaps with item 1: collapsing the per-trackpoint searches first would shrink what is left to parallelize, so it is worth doing in that order.

### Evaluating These

Items 1 and 2 change which path wins in ambiguous cases, so they should be evaluated against match quality and not only against wall-clock time. With the ordering now deterministic, a diff of `avl_matched.csv` across a change is meaningful evidence, and an unchanged expansion count is the strongest signal that a change was purely mechanical.

## Measurement Notes

The figures in [Results](#results) were taken under `cProfile`, which inflates absolute times by roughly 2.2× for this workload — matching measures 4.56 s profiled against 2.08 s of wall clock. The before/after ratios there are sound; the absolute numbers are not comparable to a plain run, which is why the current figures are reported separately in [Where the Time Goes Now](#where-the-time-goes-now).

To compare two revisions, run

```
python3 avl_osm_sample.py
```

on each side and compare `avl_matched.csv`. Pinning `PYTHONHASHSEED` is no longer necessary, though it remains harmless. If a comparison ever shows spurious differences again, suspect that a new unordered container has entered the map-building path.
