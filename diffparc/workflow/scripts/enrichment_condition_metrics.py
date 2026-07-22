#!/usr/bin/env python3
"""Emit a one-row metrics CSV for a single enrichment_sweep condition.

Pure reporting/logging: reads values already produced by MRtrix (SIFT2 mu, the
SIFT2 per-iteration cost CSV, and streamline counts from .tck headers) and writes
one row. No connectivity math, no tractogram modification.

Columns: subject,condition,enrichment_level,roi_streamlines,wb_streamlines,sift2_mu,sift2_cost
"""

import argparse
import csv
import os
import sys


def tck_count(path):
    """Read the streamline count from an MRtrix .tck ASCII header (cheap)."""
    if not path or not os.path.exists(path):
        return ""
    try:
        with open(path, "rb") as f:
            header = b""
            while b"END" not in header and len(header) < 100000:
                chunk = f.read(4096)
                if not chunk:
                    break
                header += chunk
        text = header.split(b"END", 1)[0].decode("utf-8", errors="ignore")
        for line in text.splitlines():
            if line.lower().startswith("count:"):
                return line.split(":", 1)[1].strip()
    except Exception:
        return ""
    return ""


def read_mu(path):
    if not path or not os.path.exists(path):
        return ""
    try:
        with open(path) as f:
            return f.read().strip().split()[0]
    except Exception:
        return ""


def read_final_cost(csv_path):
    """Final SIFT2 cost = last data row of the tcksift2 -csv file.

    Robust to column-name variation across MRtrix versions: prefer a header
    containing both 'total' and 'cost', else any header containing 'cost'.
    """
    if not csv_path or not os.path.exists(csv_path):
        return ""
    try:
        with open(csv_path) as f:
            rows = [r for r in csv.reader(f) if r and any(c.strip() for c in r)]
        if len(rows) < 2:
            return ""
        header = [h.strip().lower() for h in rows[0]]
        col = None
        for i, h in enumerate(header):
            if "cost" in h and "total" in h:
                col = i
                break
        if col is None:
            for i, h in enumerate(header):
                if "cost" in h:
                    col = i
                    break
        if col is None:
            return ""
        last = rows[-1]
        if col < len(last):
            return last[col].strip()
    except Exception:
        return ""
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--subject", required=True)
    ap.add_argument("--condition", required=True, help='"wb" (baseline) or a level number.')
    ap.add_argument("--enrichment-level", required=True,
                    help='Numeric level, or "baseline" for the WB condition.')
    ap.add_argument("--wb-tck", required=True, help="WB tractogram (for constant WB count).")
    ap.add_argument("--roi-tck", default="", help="Subsampled ROI tractogram (empty for baseline).")
    ap.add_argument("--mu", required=True)
    ap.add_argument("--sift2-csv", required=True)
    ap.add_argument("--out-row", required=True)
    args = ap.parse_args()

    row = {
        "subject": args.subject,
        "condition": args.condition,
        "enrichment_level": args.enrichment_level,
        "roi_streamlines": tck_count(args.roi_tck) if args.roi_tck else "0",
        "wb_streamlines": tck_count(args.wb_tck),
        "sift2_mu": read_mu(args.mu),
        "sift2_cost": read_final_cost(args.sift2_csv),
    }

    os.makedirs(os.path.dirname(args.out_row) or ".", exist_ok=True)
    with open(args.out_row, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        w.writeheader()
        w.writerow(row)

    print("enrichment_condition_metrics: " + ",".join(f"{k}={v}" for k, v in row.items()),
          file=sys.stderr)


if __name__ == "__main__":
    main()
