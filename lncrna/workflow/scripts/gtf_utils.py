"""Minimal GTF helpers shared by the lncRNA scripts (no external deps beyond pandas)."""
import re
from collections import OrderedDict

ATTR_RE = re.compile(r'(\S+) "([^"]*)"')


def parse_attrs(s: str) -> dict:
    return dict(ATTR_RE.findall(s))


def read_gtf(path: str):
    """Yield (fields list, attrs dict) for every non-comment line."""
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            yield p, parse_attrs(p[8])


def transcripts_from_gtf(path: str):
    """Return OrderedDict transcript_id -> dict(chrom, strand, gene_id, exons=[(s,e)...], attrs, lines)."""
    tx = OrderedDict()
    for p, a in read_gtf(path):
        if p[2] not in ("exon", "transcript"):
            continue
        tid = a.get("transcript_id")
        if not tid:
            continue
        t = tx.setdefault(tid, {"chrom": p[0], "strand": p[6], "gene_id": a.get("gene_id", tid),
                                "exons": [], "attrs": a, "lines": []})
        t["lines"].append("\t".join(p))
        if p[2] == "exon":
            t["exons"].append((int(p[3]), int(p[4])))
        else:
            t["attrs"] = {**t["attrs"], **a}
    for t in tx.values():
        t["exons"].sort()
        t["length"] = sum(e - s + 1 for s, e in t["exons"])
        t["start"] = t["exons"][0][0] if t["exons"] else None
        t["end"] = t["exons"][-1][1] if t["exons"] else None
    return tx


def write_fasta_from_exons(tx: dict, genome_fa: str, out_fa: str):
    """Spliced transcript sequences using samtools-style faidx (pure python, reads whole chromosomes)."""
    seqs = {}
    cur = None
    buf = []
    with open(genome_fa) as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = "".join(buf)
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    with open(out_fa, "w") as out:
        for tid, t in tx.items():
            chrom = seqs.get(t["chrom"])
            if chrom is None:
                continue
            s = "".join(chrom[a - 1:b] for a, b in t["exons"])
            if t["strand"] == "-":
                s = s.translate(comp)[::-1]
            out.write(f">{tid}\n")
            for i in range(0, len(s), 60):
                out.write(s[i:i + 60] + "\n")
