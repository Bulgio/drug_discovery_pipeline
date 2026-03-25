"""
parse_gromacs_xvg.py
--------------------
Legge i file XVG prodotti da GROMACS per un singolo ligando e genera un JSON
con tutti i parametri rilevanti per la caratterizzazione del complesso.

Struttura attesa delle cartelle:
    <base_dir>/
        <ligand_id>/
            results/
                energy_components.xvg
                rmsd_ligand.xvg
                rmsd_protein.xvg
                rmsf_residues.xvg
                gyrate.xvg
                hbonds_num.xvg
                contacts.xvg
                sasa_protein.xvg
                sasa_residue.xvg
                ss_count.xvg
                eigenvalues.xvg

Utilizzo:
    # Singolo ligando
    python parse_gromacs_xvg.py --base_dir /path/to/output --ligand_id 972703

    # Tutti i ligandi in binding_energies_ranked.txt
    python parse_gromacs_xvg.py --base_dir /path/to/output \
                                 --ranking_file binding_energies_ranked.txt \
                                 --output_json results.json
"""

import os
import re
import json
import argparse
import numpy as np
from pathlib import Path


# ─────────────────────────────────────────────
#  Parser generico XVG
# ─────────────────────────────────────────────

def parse_xvg(filepath: str) -> dict:
    """
    Legge un file XVG di GROMACS.
    Restituisce:
        {
            "legends": [str, ...],   # etichette delle colonne dati
            "data":    np.ndarray    # shape (n_frames, n_cols)
                                     # col 0 = asse X (tempo o residuo)
        }
    """
    legends = []
    data_rows = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Legenda delle serie
            m = re.search(r'@ s\d+ legend "(.+)"', line)
            if m:
                legends.append(m.group(1))
                continue
            # Salta commenti e direttive
            if line.startswith("#") or line.startswith("@"):
                continue
            # Riga dati
            try:
                values = [float(v) for v in line.split()]
                if values:
                    data_rows.append(values)
            except ValueError:
                continue

    data = np.array(data_rows) if data_rows else np.empty((0, 0))
    return {"legends": legends, "data": data}


# ─────────────────────────────────────────────
#  Utility statistica
# ─────────────────────────────────────────────

def stats(arr: np.ndarray) -> dict:
    """Ritorna mean, std, min, max arrotondati a 4 decimali."""
    if arr.size == 0:
        return {"mean": None, "std": None, "min": None, "max": None}
    return {
        "mean": round(float(np.mean(arr)), 4),
        "std":  round(float(np.std(arr)),  4),
        "min":  round(float(np.min(arr)),  4),
        "max":  round(float(np.max(arr)),  4),
    }


# ─────────────────────────────────────────────
#  Estrattori per ciascun file XVG
# ─────────────────────────────────────────────

def extract_energy(filepath: str) -> dict:
    """
    energy_components.xvg
    Colonne attese: Time | Coulomb(SR) | Total Energy | Pressure | Box-Z
    """
    xvg = parse_xvg(filepath)
    d = xvg["data"]
    legends = xvg["legends"]
    result = {}

    col_map = {legend: i + 1 for i, legend in enumerate(legends)}

    for key, label in [
        ("coulomb_SR_kJ_mol",  "Coulomb (SR)"),
        ("total_energy_kJ_mol","Total Energy"),
        ("pressure_bar",       "Pressure"),
        ("box_z_nm",           "Box-Z"),
    ]:
        if label in col_map and d.shape[1] > col_map[label]:
            result[key] = stats(d[:, col_map[label]])
        else:
            result[key] = stats(np.array([]))

    return result


def extract_rmsd(filepath: str, label: str) -> dict:
    """
    rmsd_ligand.xvg / rmsd_protein.xvg
    Colonne: Time(ns) | RMSD(nm)
    Converte nm → Å (* 10).
    """
    xvg = parse_xvg(filepath)
    d = xvg["data"]
    if d.size == 0 or d.shape[1] < 2:
        return {label: stats(np.array([]))}
    rmsd_nm = d[:, 1]
    rmsd_A  = rmsd_nm * 10.0
    return {
        label: {
            **stats(rmsd_A),
            "unit": "Angstrom",
            "simulation_time_ns": round(float(d[-1, 0]), 3) if d.size else None,
        }
    }


def extract_rmsf(filepath: str,
                  pocket_residues: list = None) -> dict:
    """
    rmsf_residues.xvg
    Colonne: Residue | RMSF(nm)
    Converte nm → Å.
    Salva anche i valori per-residuo e, se pocket_residues è fornito,
    le stats focalizzate sul pocket.
    """
    xvg = parse_xvg(filepath)
    d = xvg["data"]
    if d.size == 0 or d.shape[1] < 2:
        return {"rmsf_residues_A": stats(np.array([]))}

    rmsf_A = d[:, 1] * 10.0
    res_ids = d[:, 0].astype(int)

    # Mappa residuo → RMSF (Å)
    per_residue = {int(res_ids[i]): round(float(rmsf_A[i]), 4)
                   for i in range(len(res_ids))}

    result = {
        "rmsf_residues_A": {
            **stats(rmsf_A),
            "unit": "Angstrom",
            "n_residues": int(d.shape[0]),
            "per_residue": per_residue,
        }
    }

    if pocket_residues:
        pocket_vals = np.array([per_residue[r] for r in pocket_residues
                                 if r in per_residue])
        pocket_map  = {r: per_residue[r] for r in pocket_residues
                       if r in per_residue}
        result["rmsf_pocket_A"] = {
            **stats(pocket_vals),
            "unit": "Angstrom",
            "residues": pocket_map,
        }

    return result


def extract_gyrate(filepath: str) -> dict:
    """
    gyrate.xvg
    Colonne: Time | Rg | Rg/sX | Rg/sY | Rg/sZ
    """
    xvg = parse_xvg(filepath)
    d = xvg["data"]
    if d.size == 0 or d.shape[1] < 2:
        return {"radius_of_gyration_nm": stats(np.array([]))}
    result = {"radius_of_gyration_nm": {**stats(d[:, 1]), "unit": "nm"}}
    if d.shape[1] >= 5:
        result["rg_x_nm"] = stats(d[:, 2])
        result["rg_y_nm"] = stats(d[:, 3])
        result["rg_z_nm"] = stats(d[:, 4])
    return result


def extract_hbonds(filepath: str) -> dict:
    """
    hbonds_num.xvg
    Colonne: Time | H-bonds | Pairs within 0.35nm
    """
    xvg = parse_xvg(filepath)
    d = xvg["data"]
    if d.size == 0 or d.shape[1] < 2:
        return {"hbonds": stats(np.array([]))}
    result = {"hbonds": {**stats(d[:, 1]), "unit": "count"}}
    if d.shape[1] >= 3:
        result["hbond_pairs_035nm"] = stats(d[:, 2])
    return result


def extract_contacts(filepath: str) -> dict:
    """
    contacts.xvg  (gmx mindist)
    Colonne: Time | MinDist_Protein-LIG (nm)
    """
    xvg = parse_xvg(filepath)
    d = xvg["data"]
    if d.size == 0 or d.shape[1] < 2:
        return {"min_distance_protein_lig_nm": stats(np.array([]))}
    return {
        "min_distance_protein_lig_nm": {
            **stats(d[:, 1]),
            "unit": "nm",
        }
    }


def extract_sasa_protein(filepath: str) -> dict:
    """
    sasa_protein.xvg
    Colonne: Time | Total SASA (nm²)
    """
    xvg = parse_xvg(filepath)
    d = xvg["data"]
    if d.size == 0 or d.shape[1] < 2:
        return {"sasa_total_nm2": stats(np.array([]))}
    return {"sasa_total_nm2": {**stats(d[:, 1]), "unit": "nm^2"}}


def extract_sasa_residue(filepath: str,
                          pocket_residues: list = None) -> dict:
    """
    sasa_residue.xvg
    Colonne: Residue | Avg SASA | Std SASA (nm²)
    Salva anche i valori per-residuo e, se pocket_residues è fornito,
    le stats focalizzate sul pocket.
    """
    xvg = parse_xvg(filepath)
    d = xvg["data"]
    if d.size == 0 or d.shape[1] < 2:
        return {"sasa_per_residue_nm2": stats(np.array([]))}

    res_ids = d[:, 0].astype(int)
    avg_vals = d[:, 1]

    per_residue = {int(res_ids[i]): round(float(avg_vals[i]), 4)
                   for i in range(len(res_ids))}

    result = {
        "sasa_per_residue_nm2": {
            **stats(avg_vals),
            "unit": "nm^2",
            "n_residues": int(d.shape[0]),
            "per_residue": per_residue,
        }
    }

    if pocket_residues:
        pocket_vals = np.array([per_residue[r] for r in pocket_residues
                                 if r in per_residue])
        pocket_map  = {r: per_residue[r] for r in pocket_residues
                       if r in per_residue}
        result["sasa_pocket_nm2"] = {
            **stats(pocket_vals),
            "unit": "nm^2",
            "residues": pocket_map,
        }

    return result


def extract_ss_count(filepath: str) -> dict:
    """
    ss_count.xvg  (gmx dssp)
    Colonne: Time | Loops | Breaks | Bends | Turns | PP-II | pi | 3_10 | beta | bridges | alpha
    Calcola la media temporale di ciascun elemento di struttura secondaria.
    """
    xvg = parse_xvg(filepath)
    d = xvg["data"]
    legends = xvg["legends"]

    if d.size == 0:
        return {"secondary_structure": {}}

    # Mappa legend → colonna (col 0 = Time)
    ss_map = {
        "Loops":       "loops",
        "Breaks":      "breaks",
        "Bends":       "bends",
        "Turns":       "turns",
        "PP-Helices":  "pp2_helices",
        "pi-Helices":  "pi_helices",
        "3_10-Helices":"310_helices",
        "beta-Strands": "beta_strands",
        "beta-Bridges": "beta_bridges",
        "alpha-Helices":"alpha_helices",
    }

    # Normalizza i legend names per il match
    def normalize(s):
        s = re.sub(r'\\[a-z]\\f\{\}', '', s)
        s = re.sub(r'\\[sSnN]', '', s)
        s = s.replace('\\x', '').strip()
        return s

    result = {}
    for i, legend in enumerate(legends):
        col_idx = i + 1
        if col_idx >= d.shape[1]:
            break
        norm = normalize(legend)
        # Cerca corrispondenza parziale
        key = None
        for k, v in ss_map.items():
            if k.lower() in norm.lower() or norm.lower() in k.lower():
                key = v
                break
        if key is None:
            key = f"ss_col_{col_idx}"
        result[key] = stats(d[:, col_idx])

    return {"secondary_structure": result}


def extract_mmpbsa_csv(filepath: str) -> dict:
    """
    MMPBSA.csv  (gmx_MMPBSA)
    Struttura: blocchi separati da riga vuota, ciascuno con header e dati.
    Sezioni attese:
        - Complex Energy Terms
        - Receptor Energy Terms
        - Ligand Energy Terms
        - Delta Energy Terms   ← quella che ci interessa di più
    Colonne: Frame#, BOND, ANGLE, DIHED, VDWAALS, EEL,
             1-4 VDW, 1-4 EEL, EPB, ENPOLAR, EDISPER, GGAS, GSOLV, TOTAL

    Restituisce stats per ogni sezione, con focus su Delta (binding free energy).
    """
    import csv, io

    sections = {}
    current_section = None
    current_rows = []
    current_cols = []

    def finalize_section():
        if current_section and current_rows:
            data = np.array(current_rows, dtype=float)
            result = {}
            for i, col in enumerate(current_cols[1:], start=1):   # skip Frame#
                key = col.strip().replace(" ", "_").replace("-", "_").replace("/", "_")
                result[key] = stats(data[:, i])
            sections[current_section] = result

    with open(filepath, newline='', encoding='utf-8-sig') as f:
        content = f.read()

    # Normalizza line endings
    content = content.replace('\r\n', '\n').replace('\r', '\n')
    lines = content.split('\n')

    for line in lines:
        line = line.strip()
        if not line:
            finalize_section()
            current_section = None
            current_rows = []
            current_cols = []
            continue

        # Intestazione di sezione (non inizia con numero o "Frame")
        if line.startswith("Frame #") or line.startswith("Frame#"):
            current_cols = [c.strip() for c in line.split(',')]
            continue

        # Riga dati (inizia con un numero intero)
        if current_cols and line[0].isdigit():
            try:
                vals = [float(v) for v in line.split(',')]
                if len(vals) == len(current_cols):
                    current_rows.append(vals)
            except ValueError:
                pass
            continue

        # Titolo sezione
        if not current_cols:
            # Pulizia: rimuovi "POISSON BOLTZMANN:" e simili
            clean = line.replace("POISSON BOLTZMANN:", "").strip().rstrip(':').strip()
            if clean:
                finalize_section()
                current_section = clean.lower().replace(" ", "_")
                current_rows = []
                current_cols = []

    finalize_section()  # flush ultimo blocco

    # Rinomina le sezioni per chiarezza
    rename = {
        "complex_energy_terms": "complex",
        "receptor_energy_terms": "receptor",
        "ligand_energy_terms": "ligand",
        "delta_energy_terms": "delta",
    }
    result = {}
    for k, v in sections.items():
        result[rename.get(k, k)] = v

    # Estrai un summary flat delle delta (le più importanti)
    summary = {}
    if "delta" in result:
        d = result["delta"]
        for key in ["VDWAALS", "EEL", "GGAS", "GSOLV", "TOTAL",
                    "EPB", "ENPOLAR"]:
            k2 = key.replace("-", "_")
            if k2 in d:
                summary[f"delta_{k2}"] = d[k2]

    return {
        "sections": result,
        "summary_kJ_mol": summary,
        "method": "PB",
    }


def extract_hbmap(filepath: str) -> dict:
    """
    hbmap.xpm  (gmx hbond -hbm)
    Matrice binaria: righe = H-bond index, colonne = frame temporali.
    Carattere 'o' (o qualsiasi non-spazio dopo i colori) = H-bond presente.

    Estrae per ogni H-bond (riga):
        - frazione di tempo in cui è presente (occupancy)
    E globalmente:
        - n_hbonds_unique   : numero di H-bond unici osservati almeno una volta
        - mean_occupancy    : occupancy media su tutti gli H-bond osservati
        - max_occupancy     : H-bond più persistente
        - hbond_occupancies : lista ordinata delle occupancy (0-1) per H-bond
    """
    present_char = None
    n_cols = None
    n_rows = None
    rows = []

    with open(filepath, "r") as f:
        for line in f:
            # Dimensioni: "n_cols n_rows n_colors n_chars"
            if n_cols is None:
                m = re.search(r'"(\d+)\s+(\d+)\s+\d+\s+\d+"', line)
                if m:
                    n_cols = int(m.group(1))
                    n_rows = int(m.group(2))
                    continue

            # Colori: trova il carattere che mappa su "Present"
            if present_char is None:
                m = re.search(r'"(.)\s+c\s+#\S+\s*"\s*/\*\s*"Present"', line)
                if m:
                    present_char = m.group(1)
                    continue

            # Righe dati (stringhe tra virgolette dopo i commenti header)
            if n_cols and present_char:
                m = re.match(r'^"(.+)"[,;]?\s*$', line)
                if m:
                    row_str = m.group(1)
                    if len(row_str) == n_cols:
                        rows.append(row_str)

    if not rows or present_char is None:
        return {"hbmap": {"error": "parse failed"}}

    # Calcola occupancy per ogni H-bond (riga)
    occupancies = []
    for row in rows:
        count = sum(1 for c in row if c == present_char)
        occupancies.append(count / n_cols)

    occupancies = np.array(occupancies)
    # H-bond osservati almeno una volta
    active = occupancies[occupancies > 0]

    return {
        "hbmap": {
            "n_hbond_slots":      int(len(occupancies)),
            "n_hbonds_unique":    int(len(active)),
            "mean_occupancy":     round(float(np.mean(active)), 4) if active.size else 0.0,
            "max_occupancy":      round(float(np.max(active)), 4) if active.size else 0.0,
            "hbond_occupancies":  [round(float(o), 4) for o in sorted(occupancies[occupancies > 0], reverse=True)],
            "n_frames":           int(n_cols),
        }
    }


def extract_eigenvalues(filepath: str) -> dict:
    """
    eigenvalues.xvg  (gmx covar – PCA)
    Colonne: EigIndex | Eigenvalue (nm²)
    Estrae i primi 3 autovalori e la varianza spiegata cumulativa.
    """
    xvg = parse_xvg(filepath)
    d = xvg["data"]
    if d.size == 0 or d.shape[1] < 2:
        return {"pca_eigenvalues": {}}

    eigenvalues = d[:, 1]
    total_var   = float(np.sum(eigenvalues))
    top3        = eigenvalues[:3].tolist() if len(eigenvalues) >= 3 else eigenvalues.tolist()
    cumvar_top3 = float(np.sum(eigenvalues[:3])) / total_var * 100 if total_var > 0 else 0

    return {
        "pca_eigenvalues": {
            "top3_nm2":          [round(v, 4) for v in top3],
            "cumulative_var_top3_pct": round(cumvar_top3, 2),
            "total_variance_nm2": round(total_var, 4),
            "unit": "nm^2",
        }
    }


# ─────────────────────────────────────────────
#  Funzione principale per un ligando
# ─────────────────────────────────────────────

def process_ligand(ligand_id: str, base_ligand_dir: str,
                   binding_energy: float = None,
                   rank: int = None,
                   pocket_residues: list = None) -> dict:
    """
    Struttura attesa:
        <base_dir>/<ligand_id>/
            results/
                xvg/     <- file .xvg
                xpm/     <- file .xpm (hbmap, ecc.)
            MMPBSA.csv   <- direttamente nella cartella del ligando
            MMPBSA.dat

    pocket_residues: lista di interi con i residui del sito di binding,
                     es. [46,47,48,49,50,242,243,244,245,246]
                     Se fornita, aggiunge sezioni *_pocket focalizzate.
    """
    base     = Path(base_ligand_dir)
    xvg_dir  = base / "results" / "xvg"
    xpm_dir  = base / "results" / "xpm"
    mmdir    = base                         # MMPBSA files stanno qui

    record = {
        "ligand_id":      str(ligand_id),
        "rank":           rank,
        "binding_energy": {
            "mean_kJ_mol": binding_energy,
            "source":      "binding_energies_ranked.txt",
        },
    }

    # (filename, extractor, section, directory)
    extractors = [
        ("energy_components.xvg", extract_energy,                              "energy_components", xvg_dir),
        ("rmsd_ligand.xvg",       lambda p: extract_rmsd(p, "rmsd_ligand_A"),  "rmsd_ligand",       xvg_dir),
        ("rmsd_protein.xvg",      lambda p: extract_rmsd(p, "rmsd_protein_A"), "rmsd_protein",      xvg_dir),
        ("rmsf_residues.xvg",     lambda p: extract_rmsf(p, pocket_residues),  "rmsf",              xvg_dir),
        ("gyrate.xvg",            extract_gyrate,                              "gyrate",            xvg_dir),
        ("hbonds_num.xvg",        extract_hbonds,                              "hbonds",            xvg_dir),
        ("contacts.xvg",          extract_contacts,                            "contacts",          xvg_dir),
        ("sasa_protein.xvg",      extract_sasa_protein,                        "sasa_protein",      xvg_dir),
        ("sasa_residue.xvg",      lambda p: extract_sasa_residue(p, pocket_residues), "sasa_residue", xvg_dir),
        ("ss_count.xvg",          extract_ss_count,                            "ss_count",          xvg_dir),
        ("eigenvalues.xvg",       extract_eigenvalues,                         "pca",               xvg_dir),
        ("hbmap.xpm",             extract_hbmap,                               "hbmap",             xpm_dir),
        ("MMPBSA.csv",            extract_mmpbsa_csv,                          "mmpbsa",            mmdir),
    ]

    for filename, extractor, section, directory in extractors:
        filepath = directory / filename
        if filepath.exists():
            try:
                record[section] = extractor(str(filepath))
            except Exception as e:
                record[section] = {"error": str(e)}
        else:
            record[section] = {"error": f"{filename} not found in {directory}"}

    return record


# ─────────────────────────────────────────────
#  CLI
# ─────────────────────────────────────────────

def load_ranking(ranking_file: str) -> list:
    """Legge binding_energies_ranked.txt → [(ligand_id, energy), ...]"""
    entries = []
    with open(ranking_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                entries.append((parts[0], float(parts[1])))
    return entries


def merge_gromacs_into_selection(selection_json: str, base_dir: str,
                                  ranking_file: str = None,
                                  pocket_residues: list = None) -> list:
    """
    Legge il JSON di selezione (con id, smiles, scores, ecc.) e aggiunge
    i dati GROMACS estratti dagli XVG come chiave 'gromacs' in ogni entry.

    Se ranking_file è fornito, aggiunge anche binding_energy dal file di ranking.
    Se la cartella results di un ligando non esiste, aggiunge gromacs: null.
    """
    with open(selection_json) as f:
        selection = json.load(f)

    # Carica binding energies se disponibili
    energy_map = {}
    rank_map = {}
    if ranking_file and Path(ranking_file).exists():
        entries = load_ranking(ranking_file)
        for rank, (ligand_id, energy) in enumerate(entries, start=1):
            energy_map[ligand_id] = energy
            rank_map[ligand_id]   = rank

    base = Path(base_dir)
    total = len(selection)

    for i, entry in enumerate(selection, start=1):
        ligand_id   = str(entry["id"])
        base_ligand = base / ligand_id
        print(f"[{i:>3}/{total}] {ligand_id} ...", end=" ")

        if not (base_ligand / "results").exists():
            print("⚠ cartella results non trovata, skipping")
            entry["gromacs"] = None
            continue

        gromacs_data = process_ligand(
            ligand_id       = ligand_id,
            base_ligand_dir = str(base_ligand),
            binding_energy  = energy_map.get(ligand_id),
            rank            = rank_map.get(ligand_id),
            pocket_residues = pocket_residues,
        )
        # Rimuovi ligand_id ridondante (già in entry["id"])
        gromacs_data.pop("ligand_id", None)

        entry["gromacs"] = gromacs_data
        print("✓")

    return selection


def main():
    parser = argparse.ArgumentParser(
        description="Estrae dati XVG GROMACS e li aggiunge a un JSON di selezione."
    )
    parser.add_argument("--base_dir",        required=True,
                        help="Directory base con le cartelle <ligand_id>/results/")
    parser.add_argument("--output_json",     default="gromacs_results.json",
                        help="File JSON di output (default: gromacs_results.json)")

    # Modalità A: merge nel JSON di selezione (raccomandato)
    parser.add_argument("--selection_json",  default=None,
                        help="JSON di selezione (updated_selection.json) da arricchire")
    parser.add_argument("--ranking_file",    default=None,
                        help="binding_energies_ranked.txt per aggiungere binding energy e rank")

    parser.add_argument("--pocket_residues", default=None,
                        help="Residui del pocket, es. '46,47,48,49,50,242,243,244,245,246'")

    # Modalità B: singolo ligando standalone
    parser.add_argument("--ligand_id",       default=None,
                        help="Processa un solo ligando senza JSON di selezione")

    args = parser.parse_args()
    base_dir = Path(args.base_dir)
    pocket_residues = (
        [int(r.strip()) for r in args.pocket_residues.split(",")]
        if args.pocket_residues else None
    )

    # ── Modalità A: arricchisce il JSON di selezione ──────────────────────────
    if args.selection_json:
        result = merge_gromacs_into_selection(
            selection_json  = args.selection_json,
            base_dir        = str(base_dir),
            ranking_file    = args.ranking_file,
            pocket_residues = pocket_residues,
        )
        output_path = Path(args.output_json)
        with open(output_path, "w") as f:
            json.dump(result, f, indent=2)
        print(f"\n✓ JSON arricchito salvato in: {output_path}")
        print(f"  Entry totali: {len(result)}")
        print(f"  Con dati GROMACS: {sum(1 for e in result if e.get('gromacs'))}")

    # ── Modalità B: singolo ligando standalone ────────────────────────────────
    elif args.ligand_id:
        record = process_ligand(
            ligand_id       = args.ligand_id,
            base_ligand_dir = str(base_dir / args.ligand_id),
            binding_energy  = None,
            rank            = None,
            pocket_residues = pocket_residues,
        )
        output_path = Path(args.output_json)
        with open(output_path, "w") as f:
            json.dump([record], f, indent=2)
        print(f"\n✓ JSON salvato in: {output_path}")

    # ── Modalità C: tutti i ligandi dal ranking (legacy) ─────────────────────
    elif args.ranking_file:
        entries = load_ranking(args.ranking_file)
        all_records = []
        for rank, (ligand_id, energy) in enumerate(entries, start=1):
            print(f"[{rank:>3}/{len(entries)}] {ligand_id} ...")
            record = process_ligand(ligand_id, str(base_dir / ligand_id), energy, rank, pocket_residues)
            all_records.append(record)
        output_path = Path(args.output_json)
        with open(output_path, "w") as f:
            json.dump(all_records, f, indent=2)
        print(f"\n✓ JSON salvato in: {output_path}")

    else:
        parser.error("Specifica --selection_json, --ligand_id, oppure --ranking_file.")


if __name__ == "__main__":
    main()
