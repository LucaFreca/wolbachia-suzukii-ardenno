#!/usr/bin/env python3
"""
Download archived observations of a Rete CML station (default: Ardenno, SO)
from the Centro Meteo Lombardo "Cartina temporeale" snapshot archive.

Endpoint (found in /realtime/js/framework.js):
    http://www.centrometeolombardo.com/Moduli/refx.php?t=arc&d=YYYYMMDDHHMM
Each snapshot (~185 KB) contains ALL network stations at that time:
  - top.coords        -> station list (id, name, prov, x, y, ..., elevation, page id)
  - top.datostazione  -> one row per station, same order as coords
Archive starts 2013-01-01, 5-minute granularity, times in Italian LOCAL time.

Usage:
    pip install requests pandas
    python cml_ardenno_download.py --start 2025-06-01 --end 2026-11-02 --step 60
    python cml_ardenno_download.py --step 1440 --hhmm 2355     # 1 snapshot/day (daily extremes)

Be polite: default 1 request / 1.5 s. Data are (c) CML - private use only,
cite "Dati Rete CML - www.centrometeolombardo.com". Ask CML before redistributing.
"""
import argparse, csv, os, re, time, datetime as dt
import requests

URL = "http://www.centrometeolombardo.com/Moduli/refx.php?t=arc&d={}"

COLS = {
    0: "status",            
    1: "delay",             
    2: "date", 3: "time",   
    4: "temp_C", 5: "temp_max_C", 6: "temp_max_time", 7: "temp_min_C", 8: "temp_min_time",
    9: "rh_pct", 10: "rh_max_pct", 11: "rh_max_time", 12: "rh_min_pct", 13: "rh_min_time",
    14: "dewpt_C", 15: "dewpt_max_C", 16: "dewpt_max_time", 17: "dewpt_min_C", 18: "dewpt_min_time",
    19: "windchill_C", 20: "windchill_min_C", 21: "windchill_min_time",
    22: "heatindex_C", 23: "heatindex_max_C", 24: "heatindex_max_time",
    25: "wind_kmh", 26: "gust_max_kmh", 27: "gust_max_time", 28: "wind_10min_kmh",
    29: "wind_dir_deg", 30: "wind_dir",
    31: "press_hPa", 32: "press_trend", 33: "press_max_hPa", 34: "press_max_time",
    35: "press_min_hPa", 36: "press_min_time",
    37: "rain_rate_mmh",
    38: "rain_col38", 39: "rain_col39", 
    40: "rain_year_mm",
    41: "rain_col41", 42: "rain_rate_max_mmh",
    43: "col43", 44: "col44",
}
HEADER = ["snapshot"] + [COLS[i] for i in range(45)]


def parse_snapshot(html, station):
    """Return the row (list of str) for `station` or None."""
    coords = re.search(r"top\.coords\s*=\s*(\[.*?\]);", html, re.S)
    data = re.search(r"top\.datostazione\s*=\s*(\[.*?\]);", html, re.S)
    if not (coords and data):
        return None
    row_re = re.compile(r"\[('(?:[^'\\]|\\.)*'(?:\s*,\s*'(?:[^'\\]|\\.)*')*)\]")
    split = lambda s: [x.replace("\\'", "'") for x in re.findall(r"'((?:[^'\\]|\\.)*)'", s)]
    clist = [split(m.group(1)) for m in row_re.finditer(coords.group(1))]
    dlist = [split(m.group(1)) for m in row_re.finditer(data.group(1))]
    ids = [c[0] for c in clist]
    if station not in ids:
        return None
    i = ids.index(station)
    return dlist[i] if i < len(dlist) else None


def daterange(start, end, step_min, hhmm=None):
    t = dt.datetime.combine(start, dt.time(int(hhmm[:2]), int(hhmm[2:]))) if hhmm \
        else dt.datetime.combine(start, dt.time(0, 0))
    stop = dt.datetime.combine(end, dt.time(23, 59))
    while t <= stop:
        yield t
        t += dt.timedelta(minutes=step_min)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--station", default="ardenno")
    ap.add_argument("--start", default="2025-06-01")
    ap.add_argument("--end", default=dt.date.today().isoformat())
    ap.add_argument("--step", type=int, default=60, help="minutes (multiple of 5)")
    ap.add_argument("--hhmm", default=None, help="fixed first time, e.g. 2355 for daily")
    ap.add_argument("--sleep", type=float, default=1.5)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    start, end = dt.date.fromisoformat(a.start), dt.date.fromisoformat(a.end)
    out = a.out or f"cml_{a.station}_{a.start}_{a.end}_{a.step}min.csv"
    done = set()
    if os.path.exists(out):                       # resume
        with open(out) as f:
            done = {r["snapshot"] for r in csv.DictReader(f)}
    new = not os.path.exists(out)
    s = requests.Session()
    s.headers["User-Agent"] = "Mozilla/5.0 (personal research download)"

    with open(out, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(HEADER)
        for t in daterange(start, end, a.step, a.hhmm):
            key = t.strftime("%Y%m%d%H%M")
            if key in done or t > dt.datetime.now():
                continue
            row = None
            for back in (0, 5, 10):               # fallback if snapshot missing
                k2 = (t - dt.timedelta(minutes=back)).strftime("%Y%m%d%H%M")
                try:
                    r = s.get(URL.format(k2), timeout=30)
                    r.encoding = "windows-1252"
                    row = parse_snapshot(r.text, a.station)
                except requests.RequestException as e:
                    print(key, "ERR", e)
                time.sleep(a.sleep)
                if row:
                    break
            if row:
                w.writerow([key] + (row + [""] * 45)[:45])
                f.flush()
                print(key, row[2], row[3], "T", row[4], "RH", row[9], "P", row[31])
            else:
                print(key, "missing")


if __name__ == "__main__":
    main()
