"""Sparse original-TIFF audit of the putative Yan et al. 95-degree Q1 pore burst.

The deposited filename says 90_Q1, whereas its embedded ND2 Location has
Theta_95/Q1. Reads only TIFF header, image metadata, frame 0, and frames
4800..4847 (~28 MB total), via HTTP Range (206) with integrity checks.
Do not infer joint camera/pressure clock from matching experimental labels.
Source: doi:10.18419/DARUS-4683 (CC BY 4.0); no source TIFF redistributed.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timedelta, timezone
import json
import re
import struct
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from scipy.ndimage import gaussian_filter

SOURCE_ID = 351440
SOURCE_BYTES = 2_544_624_129
SOURCE_MD5_ARCHIVE = "cc6afebdb4acba25cb94c0f71d15f97b"
PRESSURE_ID = 351431
PRESSURE_BYTES = 725_225
URL = "https://darus.uni-stuttgart.de/api/access/datafile/"
PITCH_UM = 2.4859
PORE_WIDTH_UM = 800.0
DEPTH_UM_ASSUMED = 120.0
PUMP_UL_S = .01 * 1000 / 3600
PORE_VOLUME_UL = 800 * 800 * DEPTH_UM_ASSUMED / 1e9
EXPECTED_FRAMES = 4848
EXPECTED_DT = .033328503370285034
FIRST = 4800
LAST = 4847
BYTES_PER_FRAME = 512 * 512 * 2
SETTINGS = ((100,.4),(180,.4),(180,.5),(180,.6),(260,.6))


def get_range(url: str, start: int, length: int, expected_total: int):
    """Strict bounded 206 response; keep signed redirect URL for next ranges."""
    if start < 0 or length <= 0 or length > 30_000_000:
        raise ValueError("Range is outside narrow audit cap")
    req = Request(url, headers={
        "Range": f"bytes={start}-{start+length-1}",
        "User-Agent": "IndependentCAHScientificAudit/sparse-TIFF-1"})
    with urlopen(req, timeout=150) as reply:
        content_range = reply.headers.get("Content-Range", "")
        expected = f"bytes {start}-{start+length-1}/{expected_total}"
        if reply.status != 206 or content_range != expected:
            raise ValueError(f"Expected HTTP 206 {expected}; received "
                             f"{reply.status} {content_range}")
        data = reply.read(length+1)
        if len(data) != length:
            raise ValueError("Truncated or oversized HTTP byte range")
        return data, reply.url


def parse_first_ifd(header: bytes):
    if header[:4] != b"MM\x00*" or len(header) < 4096:
        raise ValueError("Expected Big-Endian classic TIFF header")
    pos = struct.unpack_from(">I", header, 4)[0]
    n = struct.unpack_from(">H", header, pos)[0]
    if n != 15:
        raise ValueError("Unexpected IFD layout")
    tags = {}
    for j in range(n):
        tag, kind, count, value = struct.unpack_from(">HHII", header, pos+2+12*j)
        tags[tag] = (kind, count, value)
    next_ifd = struct.unpack_from(">I", header, pos+2+12*n)[0]
    if tags[256][2] != 512 or tags[257][2] != 512 or tags[279][2] != BYTES_PER_FRAME:
        raise ValueError("Unexpected frame dimensions/stride")
    if tags[258] != (3,1,16 << 16):
        raise ValueError("Expected 16-bit original pixels")
    desc_kind, desc_size, desc_offset = tags[270]
    desc = header[desc_offset:desc_offset+desc_size].decode("ascii", "replace")
    if ("images=4848" not in desc or
        "finterval=0.033328503370285034" not in desc):
        raise ValueError("TIFF acquisition metadata mismatch")
    strip_offset = tags[273][2]
    if next_ifd != strip_offset + EXPECTED_FRAMES*BYTES_PER_FRAME:
        raise ValueError("Frames are not provably contiguous")
    return tags, strip_offset


def parse_nd2_info(text: str):
    """Recover embedded original ND2 folder and Julian-like dTimeAbsolute."""
    fields = {}
    for name in ("Location", "dTimeAbsolute", "dCalibration", "Series 0 Name"):
        match = re.search(r"^\s*" + re.escape(name) + r"\s*=\s*(.+)$",
                          text, re.MULTILINE)
        fields[name] = match.group(1).strip() if match else None
    if fields["Location"] is None or "/Theta_95/Q1/" not in fields["Location"]:
        raise ValueError("File's claimed original Theta_95/Q1 path not verified")
    if fields["dTimeAbsolute"] is None:
        raise ValueError("No time field in source ND2 metadata")
    jd = float(fields["dTimeAbsolute"])
    fields["julian_date_interpretation_utc"] = (
        datetime(1970,1,1,tzinfo=timezone.utc) +
        timedelta(days=jd-2440587.5)).isoformat()
    return fields


def nd2_information(header, tags, signed_url):
    # 50838: ImageJ metadata bytecounts (first metadata record is IJIJ header);
    # 50839: ImageJ metadata payload. Only 95-degree TIFF has an info record.
    _, count, bytecounts_offset = tags[50838]
    kind, payload_size, payload_offset = tags[50839]
    if count < 2 or kind != 1:
        raise ValueError("ImageJ metadata header missing")
    first, second = struct.unpack_from(">II", header, bytecounts_offset)
    if first != 20 or not (100_000 < second < 2_000_000):
        raise ValueError("Expected Nikon info/labl metadata, not labels-only")
    raw, _ = get_range(signed_url, payload_offset, first+second, SOURCE_BYTES)
    if raw[:20] != b"IJIJinfo\x00\x00\x00\x01labl\x00\x00\x12\xf0":
        raise ValueError("Unexpected ImageJ metadata magic")
    return parse_nd2_info(raw[first:first+second].decode("utf-16-be",
                                                         errors="replace"))


def find_pore(first):
    image = gaussian_filter(first, 1.2)
    pixels = round(PORE_WIDTH_UM/PITCH_UM)
    row = np.median(image[:,180:335],axis=1)
    ys = np.arange(20,165)
    y0 = int(ys[np.argmin(row[ys]+row[ys+pixels])])
    col = np.median(image[y0+70:y0+pixels-70,:],axis=0)
    xs = np.arange(60,175)
    x0 = int(xs[np.argmin(col[xs]+col[xs+pixels])])
    if not (110<=x0<=123 and 64<=y0<=74):
        raise ValueError(f"Unexpected measured pore-wall pair: {(x0,y0)}")
    return (x0,y0,x0+pixels,y0+pixels)


def crossing(times, values, target):
    envelope = np.maximum.accumulate(np.asarray(values, dtype=float))
    indexes = np.flatnonzero(envelope >= target)
    if not len(indexes):
        return None
    j = int(indexes[0])
    if j == 0:
        return float(times[0])
    return float(times[j-1]+(target-envelope[j-1])*
                 (times[j]-times[j-1])/(envelope[j]-envelope[j-1]))


def analyze(frames0, final48):
    if final48.shape != (48,512,512):
        raise ValueError("Must use exactly 48 late frames")
    first = np.asarray(frames0, dtype=float)
    last = np.asarray(final48[-1],dtype=float)
    x0,y0,x1,y1 = find_pore(first)
    reference = gaussian_filter(first, 1.)
    contrast = reference-gaussian_filter(last,1.)
    roi = np.zeros((512,512),bool)
    roi[y0+3:y1-3,x0+3:x1-3] = True
    exterior = np.zeros((512,512),bool)
    exterior[:y0-10,180:330] = True
    exterior[y1+10:,180:330] = True
    active = {s:roi & (contrast>s[0]) for s in SETTINGS}
    times = np.arange(FIRST,LAST+1)*EXPECTED_DT
    signals = {s:[] for s in SETTINGS}
    positions = {f:[] for f in (.4,.5,.6)}
    cy0,cy1 = y0+95,y1-95
    cx0,cx1 = x0+6,x1-6
    eligible = contrast[cy0:cy1,cx0:cx1]>180
    for img in final48:
        current = gaussian_filter(np.asarray(img,dtype=float),1.)
        drift = float(np.median((current-reference)[exterior]))
        darkening = reference-current+drift
        for key in SETTINGS:
            signals[key].append(float(np.count_nonzero(
                active[key]&(darkening>key[1]*contrast)))/int(roi.sum()))
        for frac in positions:
            water = ((darkening[cy0:cy1,cx0:cx1]>
                     frac*contrast[cy0:cy1,cx0:cx1]) & eligible)
            seen = np.sum(water,axis=0)/np.maximum(1,eligible.sum(axis=0))>.5
            matches = np.flatnonzero(seen)
            positions[frac].append(None if not len(matches)
                                   else int(cx0+matches[-1]))
    summaries = {}
    for setting,curve in signals.items():
        start,end = curve[0],curve[-1]
        t10 = crossing(times,curve,start+.1*(end-start))
        t90 = crossing(times,curve,start+.9*(end-start))
        if t10 is None or t90 is None:
            raise ValueError("Incomplete observed optical transition")
        summaries[f"{setting[0]}_{setting[1]:.1f}"]={
            "start_mask":round(start,6),"final_mask":round(end,6),
            "t10_clip_s":t10,"t90_clip_s":t90,"duration_s":t90-t10,
            "conditional_0p8_phase_gain_uL":PORE_VOLUME_UL*(end-start)*.8,
            "nominal_pump_delivery_uL":PUMP_UL_S*(t90-t10)}
    j1,j2=4822-FIRST,4824-FIRST
    displacement_px=positions[.5][j2]-positions[.5][j1]
    result={"source_id":SOURCE_ID,"source_filename_angle_deg":90,
            "inferred_angle_from_embedded_ND2_folder_deg":95,
            "pore_xyxy":(x0,y0,x1,y1),"last48_frame_start":FIRST,
            "frame_dt_s":EXPECTED_DT,"threshold_summary":summaries,
            "front_interval_frames":[4822,4824],
            "front_displacement_um":displacement_px*PITCH_UM,
            "front_speed_mm_s":displacement_px*PITCH_UM/
                               (2*EXPECTED_DT)/1000,
            "caveats":["No full-file MD5 calculated from sparse ranges",
                       "2D mask is not calibrated 3D liquid volume",
                       "Camera/pressure time origins cannot be equated"]}
    rows=[]
    for j in range(48):
        rows.append([FIRST+j,round(float(times[j]),6)] +
                    [round(signals[s][j],6) for s in ((180,.4),(180,.5),(180,.6))] +
                    [positions[f][j] for f in (.4,.5,.6)])
    return result,rows


def run():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-ranges",action="store_true",
                        help="Explicitly permit ~28 MB of public source byte-range reads")
    parser.add_argument("--output-dir",type=Path,default=Path(__file__).parent/"outputs")
    args=parser.parse_args()
    if not args.download_ranges:
        parser.error("Use --download-ranges only for manual runs; no CI downloads")
    source=URL+str(SOURCE_ID)
    header,signed=get_range(source,0,4096,SOURCE_BYTES)
    tags,strip=parse_first_ifd(header)
    info=nd2_information(header,tags,signed)
    first,_=get_range(signed,strip,BYTES_PER_FRAME,SOURCE_BYTES)
    batch,_=get_range(signed,strip+FIRST*BYTES_PER_FRAME,
                      48*BYTES_PER_FRAME,SOURCE_BYTES)
    a=np.frombuffer(first,dtype=">u2").reshape(512,512)
    b=np.frombuffer(batch,dtype=">u2").reshape(48,512,512)
    report,rows=analyze(a,b)
    pressure,_=get_range(URL+str(PRESSURE_ID),0,451,PRESSURE_BYTES)
    m=re.search(rb"Date\t([^\r\n]+)",pressure)
    if not m:
        raise ValueError("Pressure acquisition timestamp missing")
    report["provenance"]={
        "archive_md5_unverified_from_partial_download":SOURCE_MD5_ARCHIVE,
        "original_nd2_info":info,
        "pressure_source_file_id":PRESSURE_ID,
        "pressure_sensor_context_date_unverified_timezone":m.group(1).decode(),
        "clock_comparison":"ND2 dTimeAbsolute interpreted as Julian date gives "
          "2021-11-10; pressure context says 2021-11-24. No synchronized "
          "pressure-vs-volume pairing is warranted."}
    args.output_dir.mkdir(parents=True,exist_ok=True)
    import csv
    with (args.output_dir/"third_wettability_front_30hz.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["original_clip_frame","clip_t_s","mask_lo","mask_reference",
                    "mask_hi","front_x_40_px","front_x_50_px","front_x_60_px"])
        w.writerows(rows)
    (args.output_dir/"third_wettability_audit.json").write_text(
        json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    print("Only lightweight derived results saved; no original TIFF published.")


if __name__ == "__main__":
    run()
