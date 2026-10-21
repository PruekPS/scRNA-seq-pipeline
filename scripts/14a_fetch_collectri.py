"""Step 14a (optional, run once): download the CollecTRI human TF -> target
network from OmniPath and write it as the source,target,weight CSV that
scripts/14_tf_activity_ulm.py reads.

Mirrors decoupler's get_collectri(): AP-1 and NF-kB complexes are collapsed to
the pseudo-TFs 'AP1' and 'NFKB'; weight is +1 (activation) or -1 (repression),
preferring the consensus sign across CollecTRI's sources.

Needs internet. Run it where downloads are allowed (an HPC login node is fine:
it is a single small download, not compute), then point script 14 at the CSV.
Record the download date with the CSV -- OmniPath updates over time.

Usage:
    python 14a_fetch_collectri.py references/collectri_human.csv
"""
import io
import sys

import numpy as np
import pandas as pd
import requests

URL = ("https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&organisms=9606"
       "&fields=sources,references,curation_effort&license=academic")
AP1 = {"JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2"}
NFKB = {"NFKB1", "NFKB2", "RELA", "RELB", "REL"}


def norm_source(s):
    s = str(s)
    if s.startswith("COMPLEX:"):
        members = set(s.replace("COMPLEX:", "").split("_"))
        if members & AP1:
            return "AP1"
        if members & NFKB:
            return "NFKB"
        return None
    return s


def sign(row):
    if row.get("consensus_inhibition", 0) == 1 and row.get("consensus_stimulation", 0) != 1:
        return -1
    if row.get("is_inhibition", 0) == 1 and row.get("is_stimulation", 0) != 1:
        return -1
    return 1


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    out = sys.argv[1]
    r = requests.get(URL, timeout=300)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text), sep="\t")
    df["source"] = df["source_genesymbol"].map(norm_source)
    df["target"] = df["target_genesymbol"].astype(str)
    df = df[df["source"].notna() & ~df["target"].str.startswith("COMPLEX:")]
    df["weight"] = df.apply(sign, axis=1)
    net = df.groupby(["source", "target"], as_index=False)["weight"].mean()
    net["weight"] = np.sign(net["weight"]).replace(0, 1)
    net.to_csv(out, index=False)
    print(f"CollecTRI: {net.source.nunique()} TFs, {len(net):,} edges, "
          f"{int((net.weight < 0).sum()):,} repressive -> {out}")


if __name__ == "__main__":
    main()
