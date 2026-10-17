#!/usr/bin/env python
"""Baseline quantitativo da reconstrução (issue #2 do fork).

A digitalização não é reproduzível (ruído do PulseGenerator com TRandom3(0)),
então o baseline parte de um ESD fixo (Zee.ESD.root, neste diretório) e só a
reconstrução é refeita. Subcomandos:

    baseline.py run      --esd Zee.ESD.root --outdir DIR    # reco_trf + ntuple_trf
    baseline.py extract  --aod A --ntuple N --outdir DIR    # gera eventos.csv e clusters.csv
    baseline.py compare  REF_DIR NEW_DIR                    # compara com tolerâncias
    baseline.py check    --outdir DIR                       # run + extract + compare com este diretório

Precisa do ambiente do Lorenzetti (source build/lzt_setup.sh).
"""
import argparse
import csv
import math
import os
import subprocess
import sys

import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))

# Tolerâncias iniciais (ver README.md): inteiros e booleanos exatos; floats com
# |a - b| <= ATOL + RTOL * |ref|.
RTOL = 1e-5
ATOL = 1e-6

# Variáveis por cluster lidas da NTUPLE (uma linha por cluster).
CLUSTER_VARS = [
    "cl_eta", "cl_phi", "cl_e", "cl_et", "cl_e0", "cl_e1", "cl_e2", "cl_e3",
    "cl_ehad1", "cl_ehad2", "cl_ehad3", "cl_etot", "cl_e233", "cl_e237", "cl_e277",
    "cl_emaxs1", "cl_emaxs2", "cl_e2tsts1", "cl_reta", "cl_rphi", "cl_rhad",
    "cl_rhad1", "cl_eratio", "cl_f0", "cl_f1", "cl_f2", "cl_f3", "cl_weta2",
    "truth_cl_e", "truth_cl_et", "truth_cl_eta", "truth_cl_phi",
    "el_et", "el_eta", "el_phi",
]
PID_VARS = ["el_tight", "el_medium", "el_loose", "el_vloose"]
RING_VARS = ["cl_rings_Rings", "cl_ringsL0", "truth_cl_rings_TruthRings"]
# Fora do baseline até a issue #19: em clusters centrais esses membros não são
# inicializados e saem com lixo de memória.
EXCLUDED = ["cl_secondR", "cl_lambdaCenter", "cl_secondLambda", "cl_fracMax",
            "cl_lateralMom", "cl_longitudinalMom"]

# Contagens por evento lidas do AOD. O contador "<Container>_" não pode ser lido
# sem dicionário (devolve 0 ou lixo); o número de objetos é o comprimento de uma
# folha-membro.
EVENT_COUNTS = {
    "n_cells":         "CaloCellContainer_Cells.e",
    "n_truth_cells":   "CaloCellContainer_TruthCells.e",
    "n_clusters":      "CaloClusterContainer_Clusters.e",
    "n_truth_clusters": "CaloClusterContainer_TruthClusters.e",
    "n_electrons":     "ElectronContainer_Electrons.eta",
    "n_rings":         "CaloRingsContainer_Rings.cluster_link",
    "n_seeds":         "SeedContainer_Seeds.e",
}
EVENT_SUMS = {
    "sum_cells_e":       "CaloCellContainer_Cells.e",
    "sum_truth_cells_e": "CaloCellContainer_TruthCells.e",
}


def _leaf_values(tree, name):
    leaf = tree.GetLeaf(name)
    return [leaf.GetValue(j) for j in range(leaf.GetLen())]


def extract_events(aod_path):
    tf = ROOT.TFile.Open(aod_path)
    tree = tf.Get("CollectionTree")
    rows = []
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        row = {"event": int(tree.GetLeaf("EventInfoContainer_Events.eventNumber").GetValue())}
        for col, leaf in EVENT_COUNTS.items():
            row[col] = tree.GetLeaf(leaf).GetLen()
        for col, leaf in EVENT_SUMS.items():
            row[col] = math.fsum(_leaf_values(tree, leaf))
        rows.append(row)
    tf.Close()
    # a simulação multithread grava os eventos fora de ordem
    return sorted(rows, key=lambda r: r["event"])


def extract_clusters(ntuple_path):
    tf = ROOT.TFile.Open(ntuple_path)
    tree = tf.Get("physics")
    rows = []
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        row = {"event": tree.EventNumber}
        for var in CLUSTER_VARS:
            row[var] = getattr(tree, var)
        for var in PID_VARS:
            row[var] = int(bool(getattr(tree, var)))
        for var in RING_VARS:
            rings = list(getattr(tree, var))
            row[var + "_n"] = len(rings)
            row[var + "_sum"] = math.fsum(rings)
            row[var] = " ".join(repr(r) for r in rings)
        rows.append(row)
    tf.Close()
    # ordem estável: evento, depois eta (dentro do evento a ordem segue a reco)
    rows.sort(key=lambda r: (r["event"], r["cl_eta"]))
    counter = {}
    for row in rows:
        row["cluster"] = counter.get(row["event"], 0)
        counter[row["event"]] = row["cluster"] + 1
    return rows


def write_csv(rows, path, first=("event",)):
    cols = list(first) + [c for c in rows[0] if c not in first]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in row.items()})


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _is_int(text):
    try:
        int(text)
        return True
    except ValueError:
        return False


def _close(ref, new):
    a, b = float(ref), float(new)
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    return abs(a - b) <= ATOL + RTOL * abs(a)


def compare_tables(ref_rows, new_rows, key):
    """Lista de diferenças entre duas tabelas casadas por `key`."""
    diffs = []
    index = {tuple(r[k] for k in key): r for r in new_rows}
    if len(new_rows) != len(ref_rows):
        diffs.append(f"número de linhas: ref={len(ref_rows)} novo={len(new_rows)}")
    for ref in ref_rows:
        k = tuple(ref[c] for c in key)
        new = index.get(k)
        if new is None:
            diffs.append(f"{dict(zip(key, k))}: ausente")
            continue
        for col, rv in ref.items():
            nv = new.get(col)
            if nv is None:
                diffs.append(f"{dict(zip(key, k))} {col}: coluna ausente")
            elif col in RING_VARS:  # vetor de rings (pode ser vazio)
                rl, nl = rv.split(), nv.split()
                if len(rl) != len(nl) or not all(_close(x, y) for x, y in zip(rl, nl)):
                    diffs.append(f"{dict(zip(key, k))} {col}: vetor difere")
            elif _is_int(rv) and _is_int(nv):
                if rv != nv:
                    diffs.append(f"{dict(zip(key, k))} {col}: ref={rv} novo={nv}")
            elif not _close(rv, nv):
                diffs.append(f"{dict(zip(key, k))} {col}: ref={rv} novo={nv}")
    return diffs


def cmd_run(args):
    os.makedirs(args.outdir, exist_ok=True)
    aod = os.path.join(args.outdir, "Zee.AOD.root")
    ntuple = os.path.join(args.outdir, "Zee.NTUPLE.root")
    for cmd, log in [(["reco_trf.py", "-i", os.path.abspath(args.esd), "-o", aod], "reco.log"),
                     (["ntuple_trf.py", "-i", aod, "-o", ntuple], "ntuple.log")]:
        with open(os.path.join(args.outdir, log), "w") as fh:
            if subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=args.outdir) != 0:
                sys.exit(f"falhou: {' '.join(cmd)} (ver {fh.name})")
    return aod, ntuple


def cmd_extract(args):
    os.makedirs(args.outdir, exist_ok=True)
    write_csv(extract_events(args.aod), os.path.join(args.outdir, "eventos.csv"))
    write_csv(extract_clusters(args.ntuple), os.path.join(args.outdir, "clusters.csv"),
              first=("event", "cluster"))
    print(f"tabelas escritas em {args.outdir}")


def cmd_compare(args):
    ok = True
    for name, key in [("eventos.csv", ("event",)), ("clusters.csv", ("event", "cluster"))]:
        diffs = compare_tables(read_csv(os.path.join(args.ref, name)),
                               read_csv(os.path.join(args.new, name)), key)
        print(f"{name}: {'OK' if not diffs else f'{len(diffs)} diferença(s)'}")
        for d in diffs[:20]:
            print("   ", d)
        ok = ok and not diffs
    print(f"tolerâncias: inteiros exatos, floats |a-b| <= {ATOL} + {RTOL}*|ref|")
    return 0 if ok else 1


def cmd_check(args):
    args.esd = os.path.join(HERE, "Zee.ESD.root")
    args.aod, args.ntuple = cmd_run(args)
    cmd_extract(args)
    args.ref, args.new = HERE, args.outdir
    return cmd_compare(args)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="reco_trf + ntuple_trf sobre um ESD")
    p.add_argument("--esd", default=os.path.join(HERE, "Zee.ESD.root"))
    p.add_argument("--outdir", required=True)
    p = sub.add_parser("extract", help="gera eventos.csv e clusters.csv")
    p.add_argument("--aod", required=True)
    p.add_argument("--ntuple", required=True)
    p.add_argument("--outdir", required=True)
    p = sub.add_parser("compare", help="compara dois diretórios de tabelas")
    p.add_argument("ref")
    p.add_argument("new")
    p = sub.add_parser("check", help="run + extract + compare com o baseline versionado")
    p.add_argument("--outdir", required=True)
    args = parser.parse_args()
    ROOT.gErrorIgnoreLevel = ROOT.kError
    rc = {"run": lambda a: (cmd_run(a), 0)[1], "extract": lambda a: (cmd_extract(a), 0)[1],
          "compare": cmd_compare, "check": cmd_check}[args.cmd](args)
    sys.stdout.flush()
    # sair sem passar pelo teardown do PyROOT, que às vezes dá segfault com estas libs
    os._exit(rc or 0)


if __name__ == "__main__":
    main()
