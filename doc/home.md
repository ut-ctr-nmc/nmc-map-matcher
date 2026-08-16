# NMC Map Matcher

The NMC Map Matcher is designed to take a set of georeference points, such as GTFS, GPS tracks, etc., and project them to an underlying node-link representation of a transportation network. The algorithm maintains multiple candidate paths, which turns out to be very effective in dealing with a variety of ambiguities in the mapping.

Navigate to the following:

* [Paper Abstract](abstract.md)
* [Introductory Presentation](intro_pres.md)
* [Modules](modules.md)
* [Path Match](path_match.md)
* [Performance](performance.md)

This documentation continues to be updated. Also, a number of comments are provided from within the source code that may be of help.

## Usage Notes
This code had been written to support research activities within the Network Modeling Center, and continues to be an active resource for NMC activities. Version 2.0 improvements have brought about the "shapely" and "pyproj" back-ends that allow for customized georeferencing formats, and "networkx" that provides for efficient handling of large road networks. Sample codes include a free-standing "fake" set of small bus routes, and the mapping of real bus trajectories on a regional OpenStreetMap road network fetched from Overpass API.

## Installation Notes
This project appears in the TRB Annual Meeting 2015 program as a "practice-ready" paper. While this project contains functioning code, be advised that it takes some coding to apply for a particular purpose. The quickest way to get the code running with the samples on your computer is to do the following:

1. Download and install [Python 3.13 or later](https://www.python.org/downloads/).
1. Install the shapely, pyproj, and networkx Python packages.
1. Run sample codes `gtfs_small_sample.py` and `avl_osm_sample.py`. There is also `realtime_sample.py`.
1. Visualize CSV inputs and results using Google My Maps, QGIS, or other comparable software.

## Development Roadmap
This section contains a description of the branches and steps that are planned for future development. The branches include:
* **master:** Contains baseline functioning code, now Version 2.0.
* **dev:** The working area for current development efforts. These may occasionally not successfully run.

The development roadmap is to:

1. Add more free-standing samples that use more input formats.
1. Accelerate map-matching with optimizations and concurrency. A first round of optimization has been done and is written up in [Performance](performance.md), which also lists what remains; concurrency has not been started.
1. Improve documentation on how to use class methods and attributes.
1. How about creating a QGIS plug-in?