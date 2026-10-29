#!/usr/bin/env python3
"""Build the data for the chloroplast hit browser: Palomar et al. 2024 (Wierzbicki lab, PNAS,
GSE228230) tracks + our rnahub hits, Mehta 2025 Rfam families, Rfam cmscan hits, TransTermHP terminators.

Tracks (from palomar2024/cov_mq0, made by map_run.sh), replicates RPM-normalised and averaged:
  IPOD, RpoB ptChIP, IPOD/ChIP 4do mutants  coverage ratio sample/input
  ptATAC                                    Tn5 insertion (read 5' end) ratio nucleoid/naked
  SOLINA                                    coverage ratio membrane/stroma
Writes palomar2024/browser_data.json, palomar2024/hits_palomar.tsv (per-hit log2 ratios and
genome-wide percentile of each value among 50-bp bins) and palomar2024/browser.html (from browser_template.html).

Run from the project dir:  python3 palomar2024/build_browser.py
"""
import bisect
import collections
import csv
import glob
import json
import os
import re
import sys

import numpy as np
from Bio import SeqIO

sys.path.insert(0, 'mehta2025')
from mehta_pipeline import helices_genome, regions  # noqa: E402

D = 'palomar2024'
L = 154478
BIN = 50
PSEUDO = 0.5  # RPM pseudocount per nt before ratios

# track id: (label, group, numerator condition, denominator condition, signal kind)
TRACKS = [
    ('ipod', 'IPOD, 14 d Col-0', 'Protein occupancy', 'IPOD 14do Col0', 'input 14do Col0', 'cov'),
    ('chip_day', 'RpoB ptChIP, day', 'PEP binding', 'ptChIP Day', 'input Day', 'cov'),
    ('chip_night', 'RpoB ptChIP, night', 'PEP binding', 'ptChIP Night', 'input Night', 'cov'),
    ('chip_4do_col0', 'RpoB ptChIP, 4 d Col-0', 'PEP binding', 'ChIP 4do Col0', 'input 4do Col0', 'cov'),
    ('chip_4do_sig2', 'RpoB ptChIP, 4 d sig2', 'PEP binding', 'ChIP 4do sig2-2', 'input 4do sig2-2', 'cov'),
    ('chip_4do_chlm', 'RpoB ptChIP, 4 d chlm', 'PEP binding', 'ChIP 4do chlm', 'input 4do chlm', 'cov'),
    ('chip_4do_chl27', 'RpoB ptChIP, 4 d chl27', 'PEP binding', 'ChIP 4do chl27', 'input 4do chl27', 'cov'),
    ('atac', 'ptATAC, nucleoid/naked', 'DNA accessibility', 'ATAC Nucleoid', 'ATAC Naked', '5p'),
    ('ipod_4do_col0', 'IPOD, 4 d Col-0', 'Protein occupancy', 'IPOD 4do Col0', 'IPOD input 4do Col0', 'cov'),
    ('ipod_4do_ptac3', 'IPOD, 4 d ptac3', 'Protein occupancy', 'IPOD 4do pTAC3', 'IPOD input 4do pTAC3', 'cov'),
    ('solina_4do_col0', 'SOLINA, 4 d Col-0', 'Membrane association', 'SOLINA 4do Col0 Membrane', 'SOLINA 4do Col0 Stroma', 'cov'),
    ('solina_day', 'SOLINA, 14 d day', 'Membrane association', 'SOLINA 14do Day Membrane', 'SOLINA 14do Day Stroma', 'cov'),
    ('solina_night', 'SOLINA, 14 d night', 'Membrane association', 'SOLINA 14do Night Membrane', 'SOLINA 14do Night Stroma', 'cov'),
    ('solina_sig2', 'SOLINA, 4 d sig2', 'Membrane association', 'SOLINA 4do sig2-2 Membrane', 'SOLINA 4do sig2-2 Stroma', 'cov'),
    ('solina_sig6', 'SOLINA, 4 d sig6', 'Membrane association', 'SOLINA 4do sig6-1 Membrane', 'SOLINA 4do sig6-1 Stroma', 'cov'),
    ('solina_chlm', 'SOLINA, 4 d chlm', 'Membrane association', 'SOLINA 4do chlm Membrane', 'SOLINA 4do chlm Stroma', 'cov'),
    ('solina_chl27', 'SOLINA, 4 d chl27', 'Membrane association', 'SOLINA 4do chl27 Membrane', 'SOLINA 4do chl27 Stroma', 'cov'),
]


def condition(title):
    """GEO sample title -> condition key used in TRACKS (replicate dropped, spellings unified)."""
    t = re.sub(r',?\s*Rep\d+$', '', title).replace('sig22-', 'sig2-2')
    f = [x.strip() for x in t.split(',')]
    if f[0] == 'Chloroplast' and f[1] == 'ATAC':
        return f'ATAC {f[2]}'
    if f[0] == 'Chloroplast' and f[1] == 'IPOD':
        return 'IPOD 14do Col0'
    if f[0] == 'Chloroplast' and f[1] == 'input':
        return 'input 14do Col0'
    if f[0] == 'RpoB':
        return f'{f[1]} {f[2]}'
    if f[0] == 'ChIP':
        return f'{f[3]} {f[1]} {f[2]}' if f[3] == 'input' else f'ChIP {f[1]} {f[2]}'
    if f[0] == 'IPOD':
        return f'IPOD input {f[1]} {f[2]}' if f[3] == 'input' else f'IPOD {f[1]} {f[2]}'
    if f[0] == 'SOLINA':
        return ' '.join(f)
    return t  # Chloroplast, Membrane, EGS/nox: no stroma partner, unused


def load_conditions():
    runs = list(csv.DictReader(open(f'{D}/ena_runs.tsv'), delimiter='\t'))
    cond = collections.defaultdict(lambda: {'cov': [], '5p': []})
    reps = collections.Counter()
    for r in runs:
        run, c = r['run_accession'], condition(r['sample_title'])
        cov, p5 = f'{D}/cov_mq0/{run}.cov.txt', f'{D}/cov_mq0/{run}.5p.txt'
        if not (os.path.exists(cov) and os.path.exists(p5)):
            continue
        y = np.loadtxt(p5)  # written last by map_run.sh, so a full one means the run is finished
        if len(y) != L or y.sum() == 0:
            continue
        x = np.loadtxt(cov)
        n = y.sum()  # mapped reads (each read has one 5' end)
        cond[c]['cov'].append(x / n * 1e6)
        cond[c]['5p'].append(y / n * 1e6)
        reps[c] += 1
    return {c: {k: np.mean(v, axis=0) for k, v in d.items() if v} for c, d in cond.items()}, reps


def bins(a):
    return a[: L // BIN * BIN].reshape(-1, BIN).sum(axis=1)


def genes():
    rec = SeqIO.read('NC_000932.1_Athaliana_chloroplast.gb', 'genbank')
    out = []
    for f in rec.features:
        if f.type not in ('CDS', 'tRNA', 'rRNA'):
            continue
        q = f.qualifiers
        n = q.get('gene', [f.type])[0]
        geneid = next((x.split(':')[1] for x in q.get('db_xref', []) if x.startswith('GeneID:')), '')
        for p in f.location.parts:
            out.append({'name': n, 'start': int(p.start) + 1, 'end': int(p.end), 'strand': '+' if p.strand == 1 else '-', 'type': f.type,
                        'product': q.get('product', [''])[0], 'locus': q.get('locus_tag', [''])[0], 'geneid': geneid,
                        'protein': q.get('protein_id', [''])[0]})
    seen, uniq = set(), []
    for g in out:
        k = (g['name'], g['start'], g['end'])
        if k not in seen:
            seen.add(k)
            uniq.append(g)
    return uniq


def hits_and_annotations():
    regs = regions()
    summ = {k: {r['region']: r for r in csv.DictReader(open(f'rnahub_{k}_summary.tsv'), delimiter='\t')} for k in ('genus', 'full')}
    region_rows, helices = [], []
    for name, (s, e) in sorted(regs.items(), key=lambda x: x[1]):
        g, f = summ['genus'].get(name, {}), summ['full'].get(name, {})
        region_rows.append({'name': name, 'start': s, 'end': e,
                            'genus_cov': g.get('n_covarying_pairs', ''), 'genus_canon': g.get('n_cov_pairs_canonical_ge50', ''),
                            'full_cov': f.get('n_covarying_pairs', ''), 'full_canon': f.get('n_cov_pairs_canonical_ge50', '')})
        for run, jobs in (('genus', 'rnahub_jobs_genus'), ('full', 'rnahub_jobs_full')):
            for lo, hi, nbp, ncov in helices_genome(jobs, name, s) or []:
                helices.append({'region': name, 'run': run, 'start': lo, 'end': hi, 'nbp': nbp, 'ncov': ncov})
    mehta = []
    for r in csv.DictReader(open('mehta2025/mehta_vs_rnahub_loci.tsv'), delimiter='\t'):
        if r['start']:
            mehta.append({'acc': r['acc'], 'id': r['id'], 'type': r['type'].strip(';'), 'start': int(r['start']), 'end': int(r['end']),
                          'strand': r['strand'], 'evalue': r['evalue'], 'genus_helix': r['genus_helix_on_locus'], 'full_helix': r['full_helix_on_locus'],
                          'region': r['our_region'] if r['our_region'] != '-' else ''})
    fams = {r['acc']: r for r in csv.DictReader(open('mehta2025/mehta_families.tsv'), delimiter='\t')}
    mehta_nohit = [{'acc': a, 'id': r['id']} for a, r in fams.items() if a not in {m['acc'] for m in mehta}]
    rfam = []
    for l in open('rfam_cmscan.tblout'):
        if l.startswith('#'):
            continue
        f = l.split()
        if f[19] == '=':
            continue
        s, e = int(f[9]), int(f[10])
        rfam.append({'id': f[1], 'acc': f[2], 'start': min(s, e), 'end': max(s, e), 'strand': f[11], 'evalue': f[17]})
    term = [dict(r, start=int(r['start']), end=int(r['end']), confidence=int(r['confidence']))
            for r in csv.DictReader(open('terminators/transterm_terminators.tsv'), delimiter='\t')]
    return region_rows, helices, mehta, mehta_nohit, rfam, term


def main():
    cond, reps = load_conditions()
    print('conditions with data:', dict(reps))
    tracks, binned = [], {}
    for tid, label, group, num, den, kind in TRACKS:
        if num not in cond or den not in cond:
            print('skip track', tid, '(missing', num if num not in cond else den, ')')
            continue
        a, b = cond[num][kind], cond[den][kind]
        bn = np.log2((bins(a) + PSEUDO * BIN) / (bins(b) + PSEUDO * BIN))
        binned[tid] = (a, b, np.sort(bn))
        tracks.append({'id': tid, 'label': label, 'group': group, 'kind': kind, 'reps': [reps[num], reps[den]],
                       'values': [round(float(v), 3) for v in bn]})

    def locus_values(s, e):
        out = {}
        for tid, (a, b, srt) in binned.items():
            n = e - s + 1
            v = float(np.log2((a[s - 1:e].sum() + PSEUDO * n) / (b[s - 1:e].sum() + PSEUDO * n)))
            out[tid] = [round(v, 3), round(100 * bisect.bisect_left(srt, v) / len(srt), 1)]
        return out

    region_rows, helices, mehta, mehta_nohit, rfam, term = hits_and_annotations()
    hits = []
    hits.append({'kind': 'ours', 'name': 'ndhD–psaC hairpin', 'region': 'intergenic_117168-117317_ndhD_psaC', 'start': 117195, 'end': 117257, 'strand': '-',
                 'note': 'rnahub: 13-bp stem, 4 cov bp (genus) / 2 (full); RNAfold -16.9 kcal/mol, z=-5.1 vs 1000 dinucleotide shuffles; no TransTermHP terminator'})
    for h in helices:
        if h['run'] == 'genus' and h['ncov'] >= 2 and h['nbp'] >= 4:
            hits.append({'kind': 'helix', 'name': h['region'], 'region': h['region'], 'start': h['start'], 'end': h['end'], 'strand': '',
                         'note': f"genus helix {h['nbp']} bp, {h['ncov']} covarying"})
    for m in mehta:
        hits.append({'kind': 'mehta', 'name': f"{m['id']} ({m['acc']})", 'acc': m['acc'], 'region': m['region'], 'start': m['start'], 'end': m['end'],
                     'strand': m['strand'], 'note': f"{m['type']}; E={m['evalue']}; our helix on locus: genus {m['genus_helix']}, full {m['full_helix']}"})
    for t in term:
        if t['confidence'] >= 90:
            hits.append({'kind': 'term', 'name': f"terminator after {t['upstream_gene']}", 'start': t['start'], 'end': t['end'],
                         'strand': t['strand'], 'note': f"TransTermHP conf {t['confidence']}, hairpin {t['hairpin_score']}, tail {t['tail_score']}"})
    for h in hits:
        h['values'] = locus_values(h['start'], h['end'])
    # whole-genome nhmmer conservation (genome_cons/merge.py) if present, else the non-coding rnahub alignments
    cf = 'genome_cons/conservation_genome.npz' if os.path.exists('genome_cons/conservation_genome.npz') else f'{D}/conservation_genus.npz'
    c = np.load(cf)['identity'] if os.path.exists(cf) else np.full(L, np.nan)
    cb = c[: L // 10 * 10].reshape(-1, 10)
    cons = [None if np.isnan(r).all() else round(float(np.nanmean(r)), 2) for r in cb]  # 10-nt bins, identity across genera
    data = {'length': L, 'bin': BIN, 'tracks': tracks, 'conservation': cons, 'genes': genes(), 'regions': region_rows, 'helices': helices,
            'mehta': mehta, 'mehta_nohit': mehta_nohit, 'rfam': rfam, 'terminators': term, 'hits': hits,
            # R-scape CaCoFold R2R diagrams published next to the page as r2r/{run}/{region}.svg
            'r2r': {run: sorted(os.path.basename(f)[:-4] for f in glob.glob(f'{D}/r2r/{run}/*.svg')) for run in ('genus', 'full')}}
    json.dump(data, open(f'{D}/browser_data.json', 'w'), separators=(',', ':'))
    with open(f'{D}/hits_palomar.tsv', 'w') as o:
        tids = [t['id'] for t in tracks]
        o.write('kind\tname\tstart\tend\tstrand\tnote\t' + '\t'.join(f'{t}_log2\t{t}_pct' for t in tids) + '\n')
        for h in hits:
            o.write(f"{h['kind']}\t{h['name']}\t{h['start']}\t{h['end']}\t{h['strand']}\t{h['note']}\t"
                    + '\t'.join(f"{h['values'][t][0]}\t{h['values'][t][1]}" for t in tids) + '\n')
    # stranded RNA-seq (rnaseq/map_run.sh; reads map to the transcript strand): coverage per 10-nt bin,
    # scaled to depth per 1e8 total coverage so samples are comparable
    data['rnaseq'] = []
    for f in sorted(glob.glob('rnaseq/cov/*.tsv')):
        a = np.loadtxt(f)
        k = 1e8 / (a[:, 1].sum() + a[:, 2].sum())
        binm = lambda v: [round(float(x), 1) for x in (v[: L // 10 * 10] * k).reshape(-1, 10).mean(1)]
        title = os.path.basename(f).split('.')[0].replace('_nonpolyA', '')
        data['rnaseq'].append({'label': f'RNA-seq, {title}', 'plus': binm(a[:, 1]), 'minus': binm(a[:, 2])})
    if os.path.exists('terminome/ends.json'):  # Terminome-seq RNA 5'/3' ends and TSS (terminome/build_ends.py)
        data['ends'] = json.load(open('terminome/ends.json'))
    if os.path.exists(f'{D}/promoters.json'):  # promoters_table.py: PEP promoters with sigma annotation
        data['promoters'] = json.load(open(f'{D}/promoters.json'))
    if os.path.exists(f'{D}/helix_expression.json'):  # helix_expression.py: RNA-seq on each helix, for the click card
        data['helix_expr'] = json.load(open(f'{D}/helix_expression.json'))
    if os.path.exists(f'{D}/sigma_pred.json'):  # sigma_predict.py
        data['extra'] = [json.load(open(f'{D}/sigma_pred.json'))]
    if os.path.exists(f'{D}/rscape_codon.json'):  # genome_cons/rscape_codon.py: codon-like vs RNA-like covariation
        data.setdefault('extra', []).append(json.load(open(f'{D}/rscape_codon.json')))
    for k in ('prodigal', 'orfs', 'blastx'):  # genefind/gene_tracks.py: gene predictions vs the GenBank CDS
        if os.path.exists(f'{D}/genefind_{k}.json'):
            data.setdefault('extra', []).append(json.load(open(f'{D}/genefind_{k}.json')))
    page = open(f'{D}/browser_template.html').read().replace('/*DATA*/null', json.dumps(data, separators=(',', ':')))
    open(f'{D}/browser.html', 'w').write(page)
    print(len(tracks), 'tracks,', len(hits), 'hits ->', f'{D}/browser_data.json', f'{D}/hits_palomar.tsv', f'{D}/browser.html')


if __name__ == '__main__':
    main()
