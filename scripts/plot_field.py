#!/usr/bin/env python3
"""Render a femheat node-field CSV (x,y,T) as a PNG heat map -- stdlib only.

Usage:
    plot_field.py <field.csv> <nx> <ny> <out.png> [pixels]

The CSV (written by `femheat --csv`) has a header line ``x,y,T`` followed by one
row per node in structured order (i fastest, then j). No third-party packages
are required: the image is encoded with the standard-library ``zlib`` module, so
this works anywhere Python 3 runs.
"""

import csv
import struct
import sys
import zlib

# Control points sampled from the "viridis" colormap.
_ANCHORS = [
    (0.00, (68, 1, 84)),
    (0.25, (59, 82, 139)),
    (0.50, (33, 145, 140)),
    (0.75, (94, 201, 98)),
    (1.00, (253, 231, 37)),
]


def _colormap(t):
    t = 0.0 if t < 0.0 else 1.0 if t > 1.0 else t
    for i in range(len(_ANCHORS) - 1):
        t0, c0 = _ANCHORS[i]
        t1, c1 = _ANCHORS[i + 1]
        if t <= t1:
            f = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
            return tuple(int(round(a + (b - a) * f)) for a, b in zip(c0, c1))
    return _ANCHORS[-1][1]


def _read_grid(path, nx, ny):
    values = []
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        next(reader)  # skip the header
        for row in reader:
            values.append(float(row[2]))
    if len(values) != nx * ny:
        raise SystemExit(f"expected {nx * ny} rows, got {len(values)}")
    return [values[j * nx:(j + 1) * nx] for j in range(ny)]


def _bilinear(grid, nx, ny, u, v):
    fx, fy = u * (nx - 1), v * (ny - 1)
    i0, j0 = int(fx), int(fy)
    i1, j1 = min(i0 + 1, nx - 1), min(j0 + 1, ny - 1)
    dx, dy = fx - i0, fy - j0
    top = grid[j0][i0] * (1 - dx) + grid[j0][i1] * dx
    bot = grid[j1][i0] * (1 - dx) + grid[j1][i1] * dx
    return top * (1 - dy) + bot * dy


def _write_png(path, pixels):
    height, width = len(pixels), len(pixels[0])
    raw = bytearray()
    for row in pixels:
        raw.append(0)  # per-scanline filter type 0 (none)
        for (r, g, b) in row:
            raw += bytes((r, g, b))

    def chunk(tag, data):
        body = tag + data
        return (struct.pack(">I", len(data)) + body +
                struct.pack(">I", zlib.crc32(body) & 0xFFFFFFFF))

    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    with open(path, "wb") as handle:
        handle.write(b"\x89PNG\r\n\x1a\n")
        handle.write(chunk(b"IHDR", ihdr))
        handle.write(chunk(b"IDAT", zlib.compress(bytes(raw), 9)))
        handle.write(chunk(b"IEND", b""))


def main():
    if len(sys.argv) < 5:
        raise SystemExit(__doc__)
    csv_path, nx, ny, out = (sys.argv[1], int(sys.argv[2]), int(sys.argv[3]),
                             sys.argv[4])
    size = int(sys.argv[5]) if len(sys.argv) > 5 else 600

    grid = _read_grid(csv_path, nx, ny)
    lo = min(min(row) for row in grid)
    hi = max(max(row) for row in grid)
    span = (hi - lo) or 1.0

    image = []
    for py in range(size):
        v = 1.0 - py / (size - 1)  # flip vertically so y points up
        row = []
        for px in range(size):
            u = px / (size - 1)
            t = (_bilinear(grid, nx, ny, u, v) - lo) / span
            row.append(_colormap(t))
        image.append(row)

    _write_png(out, image)
    print(f"wrote {out} ({size}x{size}), T in [{lo:.4f}, {hi:.4f}]")


if __name__ == "__main__":
    main()
