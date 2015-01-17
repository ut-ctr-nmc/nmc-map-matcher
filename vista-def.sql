-- PostgreSQL script for defining the minimal VISTA database specification for
-- using with the NMC Map Matcher. This is provided to allow general use of
-- the map-matcher without requiring the use of VISTA at all.  This database
-- contains the node-link representation of underlying maps.

-- In VISTA, it is expected that all databases follow this naming convention:
-- "user_network".

-- To use, first define an empty database, and then run the statements within
-- this script on it to create the correctly-defined tables. Then you'll need
-- to supply the node-link information that defines your underlying map.

-- Kenneth Perrine, kperrine@utexas.edu
-- Network Modeling Center, Center for Transportation Research,
--   Cockrell School of Engineering, The University of Texas at Austin 
-- Version 1.0

-- Copyright (C) 2014, The University of Texas at Austin

-- This program is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.

-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.

-- You should have received a copy of the GNU General Public License
-- along with this program.  If not, see <http://www.gnu.org/licenses/>.

BEGIN;

-- nodes contains a geographic point (x = GPS lat, y = GPS lng). All nodes
-- where type is not 1 will be ignored.
CREATE TABLE nodes(
  id INTEGER UNIQUE NOT NULL CHECK (id > 0),
  type INTEGER NOT NULL,
  x DOUBLE PRECISION NOT NULL,
  y DOUBLE PRECISION NOT NULL,
  PRIMARY KEY (id)
);

-- linkdetails is a unidirectional link that contains references to an origin
-- node and destination node. All links where type is not 1 will be ignored.
CREATE TABLE linkdetails(
  id INTEGER NOT NULL CHECK (id > 0),
  type INTEGER NOT NULL,
  source INTEGER NOT NULL,
  destination INTEGER NOT NULL,
  length REAL,
  PRIMARY KEY (id)
);

COMMIT;
