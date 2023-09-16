-- Toy network of 8 blocks of downtown Austin, TX --

CREATE TABLE nodes (
    id integer NOT NULL PRIMARY KEY,
    "type" integer NOT NULL,
    x double precision NOT NULL,
    y double precision NOT NULL
);

CREATE TABLE linkdetails (
    id integer NOT NULL PRIMARY KEY,
    "type" integer NOT NULL,
    source integer NOT NULL REFERENCES nodes(id),
    destination integer NOT NULL REFERENCES nodes(id),
    length real,
    speed real,
    capacity real,
    lanes integer
);

CREATE TABLE links (
    id integer NOT NULL PRIMARY KEY REFERENCES linkdetails(id),
    points path
);

COPY nodes (id, type, x, y) FROM stdin;
1	1	-97.7447999999999979	30.2824999999999989
2	1	-97.7437999999999931	30.2822999999999993
3	1	-97.7426576614379883	30.2819538698366735
4	1	-97.7419281005859375	30.281750039774888
5	1	-97.745211124420166	30.281546209713099
6	1	-97.7440953254699707	30.2811941396063773
7	1	-97.7430438995361328	30.2808976595165049
8	1	-97.7419924736022949	30.2805826494210137
9	1	-97.7455759048461914	30.2806752994491006
10	1	-97.744438648223877	30.2803602893536095
11	1	-97.7433872222900391	30.2800452792581183
12	1	-97.742314338684082	30.2797487991682459
13	1	-97.7404689788818359	30.2812867896344606
14	1	-97.7408766746520996	30.2802861693311414
15	1	-97.7411985397338867	30.2793967290615242
\.

COPY linkdetails (id, type, source, destination, length, speed, capacity, lanes) FROM stdin;
1	1	9	5	0.000944218948	0.583333313	2000	1
2	1	5	9	0.000944218948	0.583333313	2000	1
3	1	9	10	0.0011800779	0.583333313	2000	1
4	1	10	9	0.0011800779	0.583333313	2000	1
5	1	10	6	0.000901763153	0.583333313	2000	1
6	1	6	10	0.000901763153	0.583333313	2000	1
11	1	7	11	0.000918924692	0.583333313	2000	1
12	1	11	7	0.000918924692	0.583333313	2000	1
21	1	10	11	0.00109760091	0.583333313	2000	1
22	1	11	10	0.00109760091	0.583333313	2000	1
23	1	8	12	0.000893813965	0.583333313	2000	1
24	1	15	14	0.000945886422	0.583333313	2000	1
25	1	14	13	0.00108048914	0.583333313	2000	1
26	1	4	8	0.00116916385	0.583333313	2000	1
31	1	2	3	0.0011936262	0.583333313	2000	1
32	1	3	2	0.0011936262	0.583333313	2000	1
37	1	3	4	0.000757499656	0.583333313	2000	1
38	1	4	3	0.000757499656	0.583333313	2000	1
39	1	4	13	0.00153089408	0.583333313	2000	1
40	1	13	4	0.00153089408	0.583333313	2000	1
41	1	7	8	0.00109760091	0.583333313	2000	1
42	1	8	14	0.00115451624	0.583333313	2000	1
43	1	15	12	0.00117002591	0.583333313	2000	1
44	1	12	11	0.00111309462	0.583333313	2000	1
29	1	1	2	340	0.5	2000	1
30	1	2	1	340	0.5	2000	1
27	1	5	1	400	0.583333313	2000	1
28	1	1	5	400	0.583333313	2000	1
7	1	5	6	330	0.5	2000	1
8	1	6	5	330	0.5	2000	1
35	1	2	6	382	0.5	2000	1
36	1	6	2	382	0.5	2000	1
33	1	3	7	382	0.5	2000	1
34	1	7	3	382	0.5	2000	1
9	1	6	7	358	0.5	2000	1
10	1	7	6	358	0.5	2000	1
\.

COPY links (id, points) FROM stdin;
27	[(-97.7452111244201944,30.281546209713099),(-97.7447999999999979,30.2824999999999989)]
30	[(-97.7437999999999931,30.2822999999999993),(-97.7447999999999979,30.2824999999999989)]
29	[(-97.7447999999999979,30.2824999999999989),(-97.7437999999999931,30.2822999999999993)]
36	[(-97.7440953254699991,30.2811941396063986),(-97.7437999999999931,30.2822999999999993)]
32	[(-97.7426576614380025,30.2819538698366983),(-97.7437999999999931,30.2822999999999993)]
31	[(-97.7437999999999931,30.2822999999999993),(-97.7426576614380025,30.2819538698366983)]
34	[(-97.7430438995361044,30.2808976595165014),(-97.7426576614380025,30.2819538698366983)]
38	[(-97.7419281005858949,30.2817500397748987),(-97.7426576614380025,30.2819538698366983)]
37	[(-97.7426576614380025,30.2819538698366983),(-97.7419281005858949,30.2817500397748987)]
40	[(-97.7404689788817933,30.2812867896344997),(-97.7419281005858949,30.2817500397748987)]
1	[(-97.7455759048462056,30.2806752994491006),(-97.7452111244201944,30.281546209713099)]
28	[(-97.7447999999999979,30.2824999999999989),(-97.7452111244201944,30.281546209713099)]
8	[(-97.7440953254699991,30.2811941396063986),(-97.7452111244201944,30.281546209713099)]
10	[(-97.7430438995361044,30.2808976595165014),(-97.7440953254699991,30.2811941396063986)]
35	[(-97.7437999999999931,30.2822999999999993),(-97.7440953254699991,30.2811941396063986)]
7	[(-97.7452111244201944,30.281546209713099),(-97.7440953254699991,30.2811941396063986)]
5	[(-97.7444386482239054,30.2803602893535988),(-97.7440953254699991,30.2811941396063986)]
33	[(-97.7426576614380025,30.2819538698366983),(-97.7430438995361044,30.2808976595165014)]
12	[(-97.7433872222899964,30.2800452792581005),(-97.7430438995361044,30.2808976595165014)]
9	[(-97.7440953254699991,30.2811941396063986),(-97.7430438995361044,30.2808976595165014)]
41	[(-97.7430438995361044,30.2808976595165014),(-97.7419924736022949,30.2805826494209995)]
26	[(-97.7419281005858949,30.2817500397748987),(-97.7419924736022949,30.2805826494209995)]
2	[(-97.7452111244201944,30.281546209713099),(-97.7455759048462056,30.2806752994491006)]
4	[(-97.7444386482239054,30.2803602893535988),(-97.7455759048462056,30.2806752994491006)]
6	[(-97.7440953254699991,30.2811941396063986),(-97.7444386482239054,30.2803602893535988)]
3	[(-97.7455759048462056,30.2806752994491006),(-97.7444386482239054,30.2803602893535988)]
22	[(-97.7433872222899964,30.2800452792581005),(-97.7444386482239054,30.2803602893535988)]
21	[(-97.7444386482239054,30.2803602893535988),(-97.7433872222899964,30.2800452792581005)]
11	[(-97.7430438995361044,30.2808976595165014),(-97.7433872222899964,30.2800452792581005)]
44	[(-97.7423143386840962,30.2797487991681997),(-97.7433872222899964,30.2800452792581005)]
23	[(-97.7419924736022949,30.2805826494209995),(-97.7423143386840962,30.2797487991681997)]
43	[(-97.7411985397339009,30.2793967290614994),(-97.7423143386840962,30.2797487991681997)]
25	[(-97.7408766746520996,30.2802861693310987),(-97.7404689788817933,30.2812867896344997)]
39	[(-97.7419281005858949,30.2817500397748987),(-97.7404689788817933,30.2812867896344997)]
42	[(-97.7419924736022949,30.2805826494209995),(-97.7408766746520996,30.2802861693310987)]
24	[(-97.7411985397339009,30.2793967290614994),(-97.7408766746520996,30.2802861693310987)]
\.
