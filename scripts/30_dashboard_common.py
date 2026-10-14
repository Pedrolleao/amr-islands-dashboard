"""30 - shared machinery for the dashboard payloads (steps 31-35).

The dashboard is three self-contained HTML pages. Every byte of data they show is spliced into the
file, so the encoding is what decides whether the pages are usable. Three techniques do the work,
and all three are here:

  1. DICTIONARY ENCODING. Every repeated string becomes an integer index into a per-payload table.
     Product names and domain architectures repeat across tens of thousands of clusters, so this is
     worth several megabytes on its own.
  2. POSITIONAL RECORDS. A record is an array whose field order is declared once, never an object
     with the key names repeated per row. At 46,921 islands the key names alone would cost ~1 MB.
  3. BASE64-VARINT PACKING for the two high-cardinality layers, the 350,064 genes of the island maps
     and the 210,971 CA atoms of the structures.

The ids need none of this. `CFR-000001`..`CFR-208757` and `CFRI000001`..`CFRI066005` are dense,
gapless, sorted integer ranges, so a payload stores the bare integer and the page formats it back.

The codec below is deliberately the same one the ArchaeaHQ dashboard uses
(`.../Discovery/scripts/24_atlas_genemaps.py`), so its JavaScript decoder can be ported unchanged.
`js_decoder()` returns that decoder, and step 35 splices it into every page, so the two halves can
never drift apart.
"""
import json
import re

from amr_common import DASHBOARD, CLUSTER_PREFIX

# 5 bits of payload per character, bit 6 marks "another character follows"
B64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_"
B64_INDEX = {c: i for i, c in enumerate(B64)}
ISLAND_PREFIX = CLUSTER_PREFIX + "I"


# --------------------------------------------------------------------------- ids
def cluster_num(cluster_id):
    """CFR-014064 -> 14064. The prefix and the zero padding are page-side concerns."""
    return int(str(cluster_id).split("-")[-1].split(".")[0])


def cluster_module(cluster_id):
    """CFR-014064.m2 -> 2, and 0 for an ordinary cluster id."""
    m = re.search(r"\.m(\d+)$", str(cluster_id))
    return int(m.group(1)) if m else 0


def island_num(island_id):
    """CFRI046929 -> 46929."""
    return int(str(island_id)[len(ISLAND_PREFIX):])


def type_num(type_id):
    """CFRIT000123 -> 123."""
    return int(re.sub(r"^\D+", "", str(type_id)))


# --------------------------------------------------------------------------- varint
def vint(n):
    """Unsigned varint, 5 bits per character, bit 32 continues."""
    if n < 0:
        raise ValueError(f"vint is unsigned, got {n}")
    out = []
    while n >= 32:
        out.append(B64[(n & 31) | 32])
        n >>= 5
    out.append(B64[n])
    return "".join(out)


def zig(n):
    """Zigzag: signed -> unsigned, so small negatives stay one character."""
    return (n << 1) if n >= 0 else ((-n << 1) - 1)


def unzig(n):
    return (n >> 1) if not (n & 1) else -((n + 1) >> 1)


def varints(s):
    """Decode a whole packed string back to its unsigned integers. The mirror of vint()."""
    out, cur, shift = [], 0, 0
    for ch in s:
        v = B64_INDEX[ch]
        cur |= (v & 31) << shift
        if v & 32:
            shift += 5
        else:
            out.append(cur)
            cur, shift = 0, 0
    return out


# --------------------------------------------------------------------------- quality strings
# One printable character per residue for a 0-100 value: a third of what an integer array costs.
QUAL_LO = 35


def qual_string(values):
    return "".join(chr(QUAL_LO + int(round(min(100.0, max(0.0, float(v))) * 0.9))) for v in values)


def qual_values(s):
    return [(ord(c) - QUAL_LO) / 0.9 for c in s]


# --------------------------------------------------------------------------- dictionaries
class Dicts:
    """String -> integer index, one table per named vocabulary, emitted in first-seen order."""

    def __init__(self):
        self._order, self._index = {}, {}

    def enc(self, key, value):
        """Index of `value` in vocabulary `key`, adding it if new. Empty values encode as -1."""
        if value is None or (isinstance(value, float) and value != value):
            return -1
        value = str(value)
        if value == "" or value == "nan":
            return -1
        table = self._index.setdefault(key, {})
        if value not in table:
            table[value] = len(table)
            self._order.setdefault(key, []).append(value)
        return table[value]

    def size(self, key):
        return len(self._order.get(key, []))

    def as_dict(self):
        return self._order


# --------------------------------------------------------------------------- output
def write_payload(obj, name, folder=None):
    """Write one payload as compact JSON and report its size, which is the number that matters."""
    folder = folder or DASHBOARD
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / name
    text = json.dumps(obj, separators=(",", ":"), ensure_ascii=False)
    path.write_text(text, encoding="utf-8")
    mb = len(text.encode("utf-8")) / 1e6
    print(f"wrote {path.parent.name}/{name}  {mb:.2f} MB", flush=True)
    return mb


def js_decoder():
    """The JavaScript half of the codec, spliced into every page by step 35.

    It is generated from this file rather than written into the templates so the two halves cannot
    drift: if the Python codec changes, every page changes with it on the next build.
    """
    return """
const B64="%s";
const B64I=(()=>{const m={};for(let i=0;i<B64.length;i++)m[B64[i]]=i;return m;})();
function varints(s){const out=[];let cur=0,sh=0;
  for(let i=0;i<s.length;i++){const v=B64I[s[i]];cur|=(v&31)<<sh;
    if(v&32){sh+=5;}else{out.push(cur);cur=0;sh=0;}}
  return out;}
const unzig=n=>(n&1)?-((n+1)>>1):(n>>1);
const qualValues=s=>{const o=new Float32Array(s.length);
  for(let i=0;i<s.length;i++)o[i]=(s.charCodeAt(i)-%d)/0.9;return o;};
const clusterId=n=>"%s-"+String(n).padStart(6,"0");
const islandId=n=>"%s"+String(n).padStart(6,"0");
""" % (B64, QUAL_LO, CLUSTER_PREFIX, ISLAND_PREFIX)


if __name__ == "__main__":
    # Self-test: the codec must round-trip, or every payload built on it is silently wrong.
    vals = [0, 1, 31, 32, 1023, 1024, 65535, 1 << 20]
    assert varints("".join(vint(v) for v in vals)) == vals, "varint round-trip failed"
    signed = [0, -1, 1, -1000, 1000, -(1 << 19)]
    assert [unzig(v) for v in varints("".join(vint(zig(v)) for v in signed))] == signed, "zigzag failed"
    q = [0, 12.5, 50, 87.3, 100]
    assert max(abs(a - b) for a, b in zip(qual_values(qual_string(q)), q)) < 0.6, "quality string failed"
    assert cluster_num("CFR-014064") == 14064 and cluster_module("CFR-014064.m2") == 2
    assert island_num("CFRI046929") == 46929 and type_num("CFRIT000123") == 123
    print("codec round-trips: varint, zigzag, quality string, id helpers - all pass")
