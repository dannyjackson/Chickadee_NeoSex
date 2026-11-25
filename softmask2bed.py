import sys, gzip
fn = sys.argv[1]
fh = gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)
out = sys.stdout
stats = {"masked": 0, "total": 0}
def flush(name, start, end):
    if start is not None:
        out.write(f"{name}\t{start}\t{end}\n")
        stats["masked"] += end - start
name, pos, start = None, 0, None
for line in fh:
    if line.startswith(">"):
        flush(name, start, pos)
        name, pos, start = line[1:].split()[0], 0, None
        continue
    seq = line.rstrip("\n")
    for ch in seq:
        if ch.islower():
            if start is None:
                start = pos
        elif start is not None:
            flush(name, start, pos)
            start = None
        pos += 1
    stats["total"] += len(seq)
flush(name, start, pos)
m, t = stats["masked"], stats["total"]
print(f"soft-masked: {m} of {t} bp ({100*m/max(t,1):.1f}%)", file=sys.stderr)