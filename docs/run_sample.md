# Running the Sample

This document describes how to run the NMC Map Matcher with the sample that's stored in the repository. The roadway sample is a snippet of downtown Austin, TX that was extracted from a mesoscopic model provided by the local municipal planning organization. The GTFS sample is provided by CapMetro, the local transit agency.

To run the sample, one must upload the mesoscopic simulator map subset to a PostgreSQL database. That is done like this:

```bash
# You may need to -h hostname and -U username to what you need.
createdb -h localhost -U postgres test_map
cat samples/vista/vista_small_atx.sql | psql -U postgres -d test_map -h localhost
```

For evaluation purposes, a shapefile representations of these is located in the `samples/vista` directory that can be viewed in your favorite GIS software. To create the shapefiles, I did this:

In the `test_map` database:

```sql
CREATE EXTENSION postgis;

CREATE VIEW nodes_gis AS
  SELECT id, "type", ST_MakePoint(x, y)::geography AS geog
  FROM nodes;

CREATE VIEW links_gis AS
  SELECT ld.id, ld."type", source, destination, length, speed, capacity, lanes,
    ST_MakeLine(ST_MakePoint(n1.x, n1.y), ST_MakePoint(n2.x, n2.y))::geography AS geog
  FROM linkdetails ld, nodes n1, nodes n2
  WHERE n1.id = source
    AND n2.id = destination; 
```

Then, in the shell:

```bash
mkdir samples/vista/shp
pgsql2shp -f samples/vista/shp/nodes_gis -h localhost -U postgres -g geog test_map "SELECT * FROM nodes_gis;"
pgsql2shp -f samples/vista/shp/links_gis -h localhost -U postgres -g geog test_map "SELECT * FROM links_gis;"
```

Then, unfortunately, for the code to work, we also need a "links" table. (The sample was incomplete.) In the database:

```sql
CREATE TABLE links (
  id serial,
  points path
);
INSERT INTO links
  SELECT ld.id, ('[(' || n1.x || ',' || n1.y || '),(' || n2.x || ',' || n2.y || ')]')::path AS points
  FROM linkdetails ld, nodes n1, nodes n2
  WHERE n1.id = source
    AND n2.id = destination; 
```

Now, prepare to run on that database. This is currently done by running commands as listed in [the wiki](https://github.com/ut-ctr-nmc/nmc-map-matcher/wiki/Theory-of-Operation). Specifically:

```bash
# Perform initial path match:
python path_match.py localhost test_map postgres **** samples/gtfs/small_atx > path_match.csv

# Generate a problem report:
python problem_report.py localhost test_map postgres **** samples/gtfs/small_atx path_match.csv > report.csv

# Problem Code 0 = Great!
```
