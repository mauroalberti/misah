#!/usr/bin/env python3
"""
Regenerate export_synthetic_v5.gpkg.

The fixture beside this script is not written by hand: it comes out of qgSurf's
own exporter, through its public insert functions, so that what these tests read
is what the plugin actually writes rather than a transcription of its schema.
That is the whole point of it -- the two projects agree on a file format that no
test spans, and a hand-built copy of the DDL would drift from the original
without anything noticing.

It needs a qgSurf checkout at the path below, at or after the commit that
introduced schema v5 ("Record a line intersection as the span it is, not as one
distance"). qgis and osgeo are stubbed, since export.py imports them at module
level and neither is available outside QGIS; nothing the exporter does for these
tables touches either.

    python3 make_export_synthetic_v5.py

One profile, and every one of the thirteen tables filled -- including the case
schema v5 exists for: a line crossing the section at a point (300..300) beside
one running along it over a stretch (600..750).
"""
import sys, types, os

for name in ("qgis", "qgis.core", "qgis.PyQt", "osgeo"):
    sys.modules.setdefault(name, types.ModuleType(name))
class _Lvl: Warning=1; Critical=2; Info=0
class _Qgis: MessageLevel=_Lvl
class _Log:
    CRITICAL=2
    @staticmethod
    def logMessage(*a, **k): pass
sys.modules["qgis.core"].Qgis=_Qgis; sys.modules["qgis.core"].QgsMessageLog=_Log
sys.modules["osgeo"].ogr = types.SimpleNamespace()
QGSURF_PARENT = "/home/mauro/Documenti/projects"  # the directory holding the qgSurf checkout
sys.path.insert(0, QGSURF_PARENT)

from qgSurf.utils.geoprofiler import export as E

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "export_synthetic_v5.gpkg")
if os.path.exists(OUT): os.remove(OUT)

CRS = 'PROJCS["WGS 84 / UTM zone 33N",GEOGCS["WGS 84",DATUM["WGS_1984"]],AUTHORITY["EPSG","32633"]]'

rs = E.upsert_result_set(OUT, "synthetic v5", "one profile, every table filled", overwrite=True)
E.insert_source(OUT, rs, "profiles_layer", "sezioni", "memory://sezioni", CRS, {"note": "synthetic"})
E.insert_source(OUT, rs, "polygons_layer", "geologia", "memory://geologia", CRS, None)

pid = E.insert_profile(OUT, rs, 1, "sezione A-A'", CRS, 1000.0, 200.0, 850.0)

# topography: eleven samples over a rise, s in metres
E.insert_profile_samples(OUT, pid, [
    (i, i * 100.0, 200.0 + 650.0 * (1 - abs(i - 5) / 5.0), 600000.0 + i * 80.0, 4440000.0 + i * 60.0)
    for i in range(11)
])

E.insert_profile_vertices(OUT, pid, [
    {"vertex_ndx": 0, "kind": "start", "s": 0.0,    "x": 600000.0, "y": 4440000.0, "lon": 16.0,   "lat": 40.10, "z": 200.0},
    {"vertex_ndx": 1, "kind": "break", "s": 500.0,  "x": 600400.0, "y": 4440300.0, "lon": None,   "lat": None,  "z": None},
    {"vertex_ndx": 2, "kind": "end",   "s": 1000.0, "x": 600800.0, "y": 4440600.0, "lon": 16.009, "lat": 40.105,"z": 200.0},
])

E.insert_projected_points(OUT, pid, [
    {"label": "campione 1", "s": 250.0, "z": 520.0, "dist_to_profile": 35.0,
     "src_x": 600200.0, "src_y": 4440150.0, "src_z": 520.0, "src_fid": 11},
])

E.insert_projected_attitudes(OUT, pid, [
    {"label": "S0", "s": 400.0, "z": 690.0, "slope_degr": 32.0, "down_sense": "right",
     "src_dip_dir": 135.0, "src_dip_ang": 35.0, "dist_to_profile": 12.0,
     "src_x": 600320.0, "src_y": 4440240.0, "src_z": 690.0, "src_fid": 21},
    {"label": "S1", "s": 720.0, "z": 490.0, "slope_degr": 18.0, "down_sense": "left",
     "src_dip_dir": 315.0, "src_dip_ang": 20.0, "dist_to_profile": 48.0,
     "src_x": 600576.0, "src_y": 4440432.0, "src_z": 490.0, "src_fid": 22},
])

E.insert_projected_focal_mechanisms(OUT, pid, [
    {"label": "ML 4.3 2026-03-11", "s": 500.0, "z": -8000.0,
     "strike": 135.0, "dip": 60.0, "rake": -90.0, "profile_azimuth": 45.0,
     "dist_to_profile": 120.0, "src_x": 600400.0, "src_y": 4440300.0,
     "src_z": -8000.0, "src_fid": 31},
])

# what v5 exists for: a crossing, then a stretch run along, then another crossing
E.insert_intersected_lines(OUT, pid, [
    {"feature_id": "faglia",   "s_from": 300.0, "s_to": 300.0},
    {"feature_id": "faglia",   "s_from": 600.0, "s_to": 750.0},
    {"feature_id": "contatto", "s_from": 880.0, "s_to": 880.0},
])

E.insert_intersected_polygons(OUT, pid, [
    {"unit_name": "Flysch di Albidona", "s_from": 0.0,   "s_to": 450.0},
    {"unit_name": "Formazione di Timpa", "s_from": 450.0, "s_to": 1000.0},
])

E.write_graphical_params(OUT, rs, {"vertical_exaggeration": 2.0, "grid": True})
E.write_source_categories(OUT, rs, "line_intersections", [("faglia", "#d62728"), ("contatto", None)])
E.write_source_categories(OUT, rs, "polygon_intersections",
                          [("Flysch di Albidona", "#1f77b4"), ("Formazione di Timpa", "#2ca02c")])

print("scritta:", OUT, os.path.getsize(OUT), "byte")

import sqlite3
c = sqlite3.connect(OUT)
print("schema_version:", c.execute("SELECT value FROM gp_meta WHERE key='schema_version'").fetchone()[0])
for t in sorted(r[0] for r in c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'gp_%'")):
    print(f"  {t}: {c.execute('SELECT count(*) FROM ' + t).fetchone()[0]} righe")
print("colonne gp_intersected_lines:", [r[1] for r in c.execute("PRAGMA table_info(gp_intersected_lines)")])
