#!/usr/bin/env python
"""Convert a narrowPeak/broadPeak/BED file to the GFF format ROSE expects.
GFF columns: chrom, region_id, '', start, end, '', strand('.'), '', region_id
Usage: rose_peaks_to_gff.py <peaks.bed|narrowPeak> <out.gff>
"""
import sys

def main(inp, outp):
    with open(inp) as fh, open(outp, "w") as out:
        for i, line in enumerate(fh):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, start, end = f[0], f[1], f[2]
            name = f[3] if len(f) > 3 and f[3] not in (".", "") else "peak_%d" % (i + 1)
            # GFF (ROSE flavour): [chrom, name, '', start, end, '', '.', '', name]
            out.write("\t".join([chrom, name, "", start, end, "", ".", "", name]) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: rose_peaks_to_gff.py <peaks> <out.gff>\n")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
