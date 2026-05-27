# NMC Map Matcher
*Version 2.0*

Network Modeling Center<br>
Center for Transportation Research<br>
Cockrell School of Engineering<br>
The University of Texas at Austin

Repository: https://github.com/ut-ctr-nmc/nmc-map-matcher

Kenneth Perrine<br>
kperrine@utexas.edu

Copyright (C) 2014-2026, The University of Texas at Austin

## Overview

The NMC Map Matcher is designed to take a set of georeference points, such as GPS tracks, geographic shapes, etc., and project them to an underlying node-link representation of a geographic network, such as a map of a transportation system. As it runs, the algorithm maintains multiple "hypotheses", or candidate paths, which turns out to be very effective in dealing with a variety of ambiguities in the mapping. Features of this approach:

* Starts and ends of trajectories may occur at midpoints along links in the underlying map.
* The algorithm will often survive and highlight gaps or errors in the underlying mapping.
* Most any optimization on a Euclidian-mapped graph representing linear relations on non-negative numbers can be solved using this algorithm.
* The algorithm is structured to map-match and maintain several "best guesses" while new trajectory points arrive; the algorithm does not expect the entire trajectory to be known when it starts. *(An example is to be implemented)*
* The algorithm may be parallelized *(To be implemented)*
* Code is provided to return GPS coordinates at critical points "snapped to" the underlying map, or GPS coordinates along the underlying map at a desired spacing.
* "Out-of-the-box" sample codes include the mapping of GTFS shape subsets on a tiny network, and the mapping of bus trajectories on an OpenStreetMap representation of a major metro area.

This code accompanies this paper:

>   Perrine, Kenneth A., Alireza Khani, and Natalia Ruiz-Juri. A
   map-matching algorithm for applications in multimodal transportation
   network modeling. Transportation Research Board 94th Annual Meeting,
   Jan. 2015, Washington DC.

Entry points for this Python 3 project consists of a couple of samples that are documented at the project website on GitHub.

For documentation, please access the wiki on the GitHub site:
   https://github.com/ut-ctr-nmc/nmc-map-matcher/wiki
   
This project appears in the TRB Annual Meeting 2015 program as a
"practice-ready" paper. While this project contains functioning code, be advised that it isn't quite finished and until some further adaptations are made for your use case. The quickest way to get the code running on your computer is to do the following:

1. Download and install Python 3.13 or newer
2. Install the dependencies shapely, pyproj, and networkx
3. Follow the examples in the wiki pages on "Theory of Operation" and so forth.

## Licensing

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
