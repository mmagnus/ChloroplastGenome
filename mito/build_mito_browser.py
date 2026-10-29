#!/usr/bin/env python3
"""Genome browser page for the A. thaliana mitochondrion (NC_037304.1, Col-0, 367,808 bp), built from the chloroplast
viewer template (palomar2024/browser_template.html; genome-specific text and rows come from DATA.page / DATA.extra).

Rows: GenBank genes (CDS, tRNA, rRNA); plastid-derived insertions (MTPT, blastn of NC_000932.1 vs the mitochondrion,
mtpt_blast.tsv) as the top map row; Rfam cmscan (rfam/run_rfam.sh); tRNAscan-SE -O (trna/trnascan_O.out) marked
against the GenBank tRNAs; cosRNAs (GenBank ncRNA: clustered organellar short RNAs, PPR-protein footprints);
C-to-U editing sites (GenBank misc_feature); repeats (GenBank repeat_region); stranded RNA-seq (map_rnaseq.sh);
Terminome-seq ends and TSS (map_terminome.sh, then terminome/build_ends.py mito/ends mito/ends.json 367808).
The Palomar et al. 2024 libraries are from purified chloroplasts; their mitochondrial reads are contamination, so they
are not shown. Writes mito/browser.html and mito/browser_data.json. Run from the project dir:
    python3 mito/build_mito_browser.py
"""
import glob
import json
import os
import re

import numpy as np
from Bio import SeqIO

D = 'mito'
L = 367808
GB = SeqIO.read(f'{D}/NC_037304.1.gb', 'genbank')


def genes():
    out, seen = [], set()
    for f in GB.features:
        if f.type not in ('CDS', 'tRNA', 'rRNA'):
            continue
        q = f.qualifiers
        n = q.get('gene', [f.type])[0]
        geneid = next((x.split(':')[1] for x in q.get('db_xref', []) if x.startswith('GeneID:')), '')
        for p in f.location.parts:  # trans-spliced genes (nad1, nad2, nad5) have exons far apart
            k = (n, int(p.start), int(p.end))
            if k in seen:
                continue
            seen.add(k)
            out.append({'name': n, 'start': int(p.start) + 1, 'end': int(p.end), 'strand': '+' if p.strand == 1 else '-', 'type': f.type,
                        'product': q.get('product', [''])[0], 'locus': q.get('locus_tag', [''])[0], 'geneid': geneid,
                        'protein': q.get('protein_id', [''])[0]})
    return out


def rfam():
    out = []
    for line in open(f'{D}/rfam/mito_cmscan.tblout'):
        if line.startswith('#'):
            continue
        f = line.split()
        if f[19] == '=':  # lost the clan competition
            continue
        s, e = int(f[9]), int(f[10])
        out.append({'id': f[1], 'acc': f[2], 'start': min(s, e), 'end': max(s, e), 'strand': f[11], 'evalue': f[17]})
    return out


def trnascan(gl):
    trna = [g for g in gl if g['type'] == 'tRNA']
    out = []
    for line in open(f'{D}/trna/trnascan_O.out').readlines()[3:]:
        w = line.split()
        a, b, sc = int(w[2]), int(w[3]), float(w[8])
        s, e, st = min(a, b), max(a, b), '+' if a < b else '-'
        ov = [t for t in trna if t['start'] <= e and t['end'] >= s]
        same = [t for t in ov if t['strand'] == st]
        status = f"GenBank {same[0]['name']}" if same else (f"GenBank {ov[0]['name']} on the OTHER strand" if ov else 'not in GenBank')
        col = '--c6' if same else '--c2'
        out.append({'start': s, 'end': e, 'strand': st, 'name': f'{w[4]}-{w[5]}', 'color': col, 'w': 1 if sc >= 40 else .4, 'match': bool(same),
                    'info': f'<b>tRNAscan-SE: {w[4]}-{w[5]}</b> ({st}), {sc} bits<br>{status}' + (f'<br>{" ".join(w[9:])}' if len(w) > 9 else '')})
    return out


def gb_feats(ftype, color, name_fn, info_fn):
    out = []
    for f in GB.features:
        if f.type != ftype:
            continue
        for p in f.location.parts:
            out.append({'start': int(p.start) + 1, 'end': int(p.end), 'strand': '+' if p.strand == 1 else '-' if p.strand == -1 else '',
                        'name': name_fn(f), 'color': color, 'info': info_fn(f, p)})
    return out


def mtpt():
    """Plastid-derived segments: blastn hits of the chloroplast genome on the mitochondrion (>= 100 nt), merged."""
    iv = []
    for line in open(f'{D}/mtpt_blast.tsv'):
        qs, qe, ss, se, pid, ln, ev = line.split()
        if int(ln) >= 100:
            iv.append((min(int(ss), int(se)), max(int(ss), int(se)), float(pid), int(qs), int(qe)))
    iv.sort()
    merged = []
    for s, e, pid, qs, qe in iv:
        if merged and s <= merged[-1][1]:
            m = merged[-1]
            m[1] = max(m[1], e)
            m[2] = max(m[2], pid)
        else:
            merged.append([s, e, pid, qs, qe])
    return merged


def main():
    gl = genes()
    page_rfam = rfam()
    mt = mtpt()
    extra = [
        {'key': 'trnascan', 'label': 'tRNAscan-SE (-O)', 'color': '--c6', 'items': trnascan(gl)},
        {'key': 'cosrna', 'label': 'cosRNA (PPR footprints)', 'color': '--c5',
         'items': gb_feats('ncRNA', '--c5', lambda f: f.qualifiers.get('note', [''])[0].replace('cosRNA ', ''),
                           lambda f, p: f"<b>{f.qualifiers.get('note', [''])[0]}</b> ({int(p.end) - int(p.start)} nt)<br>clustered organellar short RNA: a small RNA protected by an RNA-binding (PPR-like) protein")},
        {'key': 'editing', 'label': 'C→U editing sites', 'color': '--c4',
         'items': gb_feats('misc_feature', '--c4', lambda f: '',
                           lambda f, p: f"<b>RNA editing</b> {int(p.start) + 1:,}<br>{f.qualifiers.get('note', [''])[0]}")},
        {'key': 'repeats', 'label': 'Repeats', 'color': '--gene',
         'items': gb_feats('repeat_region', '--gene', lambda f: f.qualifiers.get('note', [''])[0],
                           lambda f, p: f"<b>Repeat {f.qualifiers.get('note', [''])[0]}</b><br>{int(p.start) + 1:,}–{int(p.end):,} ({int(p.end) - int(p.start)} nt)")},
        {'key': 'mtpt', 'label': 'Plastid-derived (MTPT)', 'color': '--c3',
         'items': [{'start': s, 'end': e, 'strand': '', 'name': f'cp {qs:,}–{qe:,}', 'color': '--c3',
                    'info': f'<b>Plastid-derived insertion</b><br>{s:,}–{e:,} ({e - s + 1} nt), up to {pid:.0f}% identity to chloroplast {qs:,}–{qe:,}<br>reads here may come from the chloroplast'}
                   for s, e, pid, qs, qe in mt]},
    ]
    data = {'length': L, 'bin': 50, 'tracks': [], 'conservation': None, 'genes': gl, 'rfam': page_rfam, 'extra': extra, 'r2r': {}}
    data['rnaseq'] = []
    for f in sorted(glob.glob(f'{D}/rnaseq/*.tsv')):
        a = np.loadtxt(f)
        k = 1e8 / (a[:, 1].sum() + a[:, 2].sum())
        binm = lambda v: [round(float(x), 1) for x in (np.r_[v, np.zeros((-len(v)) % 10)] * k).reshape(-1, 10).mean(1)]
        title = os.path.basename(f).split('.')[0].replace('_nonpolyA', '')
        data['rnaseq'].append({'label': f'RNA-seq, {title}', 'plus': binm(a[:, 1]), 'minus': binm(a[:, 2])})
    if os.path.exists(f'{D}/ends.json'):
        data['ends'] = json.load(open(f'{D}/ends.json'))
    tr = extra[0]['items']
    n_disc = sum(not t['match'] for t in tr if t['w'] == 1)
    ntrna = sum(g['type'] == 'tRNA' for g in gl)
    sub = ('Arabidopsis thaliana mitochondrion (NC_037304.1, Col-0, 367,808 bp). Stranded RNA-seq (Di et al. 2014, SRP028304) and '
           'Terminome-seq RNA ends (Castandet et al. 2019, PRJNA533962) were remapped to the chloroplast and mitochondrion together, '
           'so reads from plastid-derived insertions (MTPT) go to their best match. Rows also show GenBank cosRNAs (PPR footprints), '
           'C→U editing sites and repeats, Rfam cmscan and tRNAscan-SE. '
           'The Palomar et al. 2024 libraries come from purified chloroplasts, so mitochondrial reads in them are contamination and are not shown. '
           '<a href="https://claude.ai/artifact/UeZy1saThZ1DuhPQN97wGy" target="_blank" rel="noopener">Chloroplast page</a>.')
    data['page'] = {
        'h1': 'ChloroplastGenome · mitochondrion', 'title': 'Arabidopsis Mitochondrion Browser', 'sub': sub,
        'url': 'https://claude.ai/artifact/VY8BpJ8nRKAYdxtMmUwHHQ', 'mapLabel': 'Plastid-derived (MTPT)', 'mapNoLabels': True, 'tick': 50000, 'jump': None,
        'map': [[s, e, f'MTPT {s:,}–{e:,}'] for s, e, *_ in mt] or [[1, L, 'mitochondrion']],
        'kv': [[len({g['name'] for g in gl if g['type'] == 'CDS'}), 'protein-coding genes'], [ntrna, 'GenBank tRNAs'],
               [n_disc, 'confident tRNAscan-SE calls not matching a GenBank tRNA'], [len(page_rfam), 'Rfam cmscan hits'],
               [len(extra[1]['items']), 'cosRNAs (PPR footprints)'], [f"{sum(e - s + 1 for s, e, *_ in mt) / 1000:.1f} kb", 'plastid-derived DNA']],
        'hint': 'Drag to pan, scroll or pinch to zoom, click a feature for links. Drag a track name (⋮⋮) to reorder. RNA-seq: + strand up (green), − strand down (purple), log scale. '
                'tRNAscan-SE calls are green when they match a GenBank tRNA and orange when they do not; faded = under 40 bits.',
        'foot': 'Reads mapped with bowtie2 --local to chloroplast + mitochondrion, multi-mapping reads kept (no MAPQ filter). '
                'Terminome 3′-end libraries read antisense, 5′-end libraries sense (as for the chloroplast).',
        'noUniprot': False,
    }
    page = open('palomar2024/browser_template.html').read().replace('/*DATA*/null', json.dumps(data, separators=(',', ':')))
    page = re.sub(r'<title>.*?</title>', '<title>Arabidopsis Mitochondrion Browser</title>', page, count=1)
    open(f'{D}/browser.html', 'w').write(page)
    json.dump(data, open(f'{D}/browser_data.json', 'w'), separators=(',', ':'))
    print(f'{len(gl)} gene parts, {len(page_rfam)} Rfam, {len(tr)} tRNAscan ({n_disc} not matching GenBank), {len(extra[1]["items"])} cosRNA, '
          f'{len(extra[2]["items"])} editing, {len(mt)} MTPT, {len(data["rnaseq"])} RNA-seq -> {D}/browser.html')


if __name__ == '__main__':
    main()
