"""Step 14a (optional, run once): download the CollecTRI human TF -> target
network from OmniPath and write it as the source,target,weight CSV that
scripts/14_tf_activity_ulm.py reads.

Mirrors decoupler's get_collectri(): AP-1 and NF-kB complexes are collapsed to
the pseudo-TFs 'AP1' and 'NFKB'; weight is +1 (activation) or -1 (repression),
preferring the consensus sign across CollecTRI's sources.

Needs internet. Run it where downloads are allowed (an HPC login node is fine:
it is a single small download, not compute), then point script 14 at the CSV.
Record the download date with the CSV -- OmniPath updates over time.

Mouse data (--organism mouse): CollecTRI is curated in human. Its genes are
converted to mouse with MGI's curated human-mouse homology table
(HOM_MouseHumanSequence.rpt), keeping one-to-one orthologs only. Human TF
footprints transferred this way recover mouse perturbations (Holland, Szalai &
Saez-Rodriguez, BBA Gene Regul Mech 2020, doi:10.1016/j.bbagrm.2019.194431).
Do not just change the letter case (FOXP3 -> Foxp3): on a real mouse dataset
that mapped ~6% of one-to-one genes to the wrong or a non-orthologous symbol
and added ~8,000 edges through one-to-many families (chemokines, Cdkn2a/b);
TF scores barely moved, but PROGENy-style pathway scores did (TNFa r = 0.40).
The pseudo-TFs AP1 and NFKB keep their names.

Usage:
    python 14a_fetch_collectri.py references/collectri_human.csv
    python 14a_fetch_collectri.py references/collectri_mouse.csv --organism mouse
    # offline: reuse a downloaded MGI table
    python 14a_fetch_collectri.py references/collectri_mouse.csv --organism mouse \
        --orthologs references/HOM_MouseHumanSequence.rpt
"""
import argparse
import io
import pathlib

import numpy as np
import pandas as pd
import requests

URL = ("https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&organisms=9606"
       "&fields=sources,references,curation_effort&license=academic")
MGI_URL = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
AP1 = {"JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2"}
NFKB = {"NFKB1", "NFKB2", "RELA", "RELB", "REL"}


def norm_source(uniprot_id, symbol):
    """OmniPath marks complexes in the UniProt column ('COMPLEX:P05412_P17535'); the
    gene-symbol column then holds the members joined by '_' ('JUN_JUND'). AP-1 and
    NF-kB complexes become the pseudo-TFs AP1 / NFKB (as decoupler does); other
    complexes are dropped."""
    symbol = str(symbol)
    if str(uniprot_id).startswith("COMPLEX:"):
        members = set(symbol.split("_"))
        if members & AP1:
            return "AP1"
        if members & NFKB:
            return "NFKB"
        return None
    return symbol


def sign(row):
    if row.get("consensus_inhibition", 0) == 1 and row.get("consensus_stimulation", 0) != 1:
        return -1
    if row.get("is_inhibition", 0) == 1 and row.get("is_stimulation", 0) != 1:
        return -1
    return 1


def parse_collectri(df):
    """OmniPath interactions table -> source,target,weight (complexes collapsed, signs +-1)."""
    df = df.copy()
    df["src"] = [norm_source(u, g) for u, g in zip(df["source"], df["source_genesymbol"])]
    df["tgt"] = df["target_genesymbol"].astype(str)
    df = df[df["src"].notna() & ~df["target"].astype(str).str.startswith("COMPLEX:")]
    df["weight"] = df.apply(sign, axis=1)
    net = df.groupby(["src", "tgt"], as_index=False)["weight"].mean().rename(columns={"src": "source", "tgt": "target"})
    net["weight"] = np.sign(net["weight"]).replace(0, 1)
    return net


def mgi_one_to_one(path_or_none):
    """human symbol -> mouse symbol for MGI homology classes with exactly one gene per species."""
    if path_or_none:
        hom = pd.read_csv(path_or_none, sep="\t", dtype=str)
    else:
        r = requests.get(MGI_URL, timeout=300)
        r.raise_for_status()
        hom = pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str)
    hom = hom.rename(columns={"DB Class Key": "cls", "Common Organism Name": "org", "Symbol": "sym"})
    hs = hom[hom["org"] == "human"][["cls", "sym"]]
    mm = hom[hom["org"].str.startswith("mouse")][["cls", "sym"]]
    pairs = hs.merge(mm, on="cls", suffixes=("_h", "_m"))
    one = pairs.groupby("cls").filter(lambda d: len(d) == 1)
    one = one[~one["sym_h"].duplicated(keep=False) & ~one["sym_m"].duplicated(keep=False)]
    return one.set_index("sym_h")["sym_m"]


def to_mouse(net, o2o):
    keep_src = {"AP1", "NFKB"}
    src = net["source"].where(net["source"].isin(keep_src), net["source"].map(o2o))
    tgt = net["target"].map(o2o)
    out = net.assign(source=src, target=tgt).dropna(subset=["source", "target"])
    return out.groupby(["source", "target"], as_index=False)["weight"].first()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("output", help="CSV to write (source,target,weight)")
    ap.add_argument("--organism", choices=["human", "mouse"], default="human")
    ap.add_argument("--orthologs", default=None,
                    help="Local MGI HOM_MouseHumanSequence.rpt (default: download it from MGI)")
    args = ap.parse_args()
    out = args.output
    r = requests.get(URL, timeout=300)
    r.raise_for_status()
    net = parse_collectri(pd.read_csv(io.StringIO(r.text), sep="\t"))
    if args.organism == "mouse":
        o2o = mgi_one_to_one(args.orthologs)
        n_h = len(net)
        net = to_mouse(net, o2o)
        print(f"MGI one-to-one orthologs: {len(o2o):,}; kept {len(net):,} of {n_h:,} human edges")
    pathlib.Path(out).parent.mkdir(parents=True, exist_ok=True)
    net.to_csv(out, index=False)
    print(f"CollecTRI: {net.source.nunique()} TFs, {len(net):,} edges, "
          f"{int((net.weight < 0).sum()):,} repressive -> {out}")


if __name__ == "__main__":
    main()
