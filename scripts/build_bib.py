"""Assemble paper/lrd.bib from ADS exports: every citation key mapped to its
validated ADS bibcode (scripts/bibcodes.json, written by ads_validate.py, plus
the corrections and additions below).

    export $(cat .env | tr -d '\r')
    python scripts/build_bib.py

Only keys that are actually cited in paper/lrd_engine.tex are written.
"""
import json
import os
import re

import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEX = os.path.join(ROOT, 'paper', 'lrd_engine.tex')
BIB = os.path.join(ROOT, 'paper', 'lrd.bib')
TOK = os.environ['ADS_TOKEN']
H = {'Authorization': 'Bearer ' + TOK, 'Content-Type': 'application/json'}

# Keys -> bibcodes found by ads_validate.py
keys = json.load(open(os.path.join(ROOT, 'scripts', 'bibcodes.json')))
# Corrections to the automatic match, and references ads_validate.py could not
# resolve (checked by hand on ADS).
keys.update({
    'Sedov:1946a': '1946JApMM..10..241S',     # Prikl. Mat. Mekh. 10, 241; ADS title "Propagation of strong shock waves"
    'Cioffi:1988a': '1988ApJ...334..252C',    # radiative SNR: t_sf ∝ E^{3/14} n^{-4/7} zeta_m^{-5/14}
})


# Fields missing from the ADS export, added only when absent
PATCH = {
    'Merritt:2013a': {'publisher': 'Princeton University Press'},
    'Netzer:2006a': {'publisher': 'Springer'},
}


def cited_keys(texpath):
    src = open(texpath, encoding='utf-8').read()
    src = re.sub(r'(?<!\\)%.*', '', src)
    out, seen = [], set()
    for m in re.finditer(r'\\[Cc]ite[a-zA-Z]*\s*(?:\[[^\]]*\]\s*)*\{([^}]+)\}', src):
        for k in m.group(1).split(','):
            k = k.strip()
            if k and k not in seen:
                seen.add(k)
                out.append(k)
    return out


cited = cited_keys(TEX)
missing_keys = [k for k in cited if k not in keys]
if missing_keys:
    raise SystemExit('cited keys with no bibcode: %s' % missing_keys)
unused = [k for k in keys if k not in cited]
if unused:
    print('not cited (dropped from the bib):', unused)
order = [(k, keys[k]) for k in cited]

bibcodes = [bc for _, bc in order]
r = requests.post('https://api.adsabs.harvard.edu/v1/export/bibtex', headers=H,
                  data=json.dumps({'bibcode': bibcodes, 'keyformat': '%R', 'maxauthor': 0, 'journalformat': 1}))
r.raise_for_status()
export = r.json()['export']
entries = {}
for m in re.finditer(r'@(\w+)\{([^,]+),(.*?)\n\}\n', export, flags=re.S):
    entries[m.group(2).strip()] = (m.group(1), m.group(3))
missing = [b for b in bibcodes if b not in entries]
print('missing from export:', missing)

out = ['%% lrd.bib -- references for "A Metal-Poor Tidal Disruption Engine: The Warm',
       '%% Molecular Equilibrium and Its Application to Little Red Dots".',
       '%% Every entry below was exported from NASA/ADS by scripts/build_bib.py and keyed',
       '%% to the bibcodes listed in scripts/bibcodes.json. Do not edit by hand.', '']
for key, bc in order:
    if bc not in entries:
        continue
    kind, body = entries[bc]
    if 'keywords = {' in body:
        body = body.replace('\n     keywords = {' + body.split('keywords = {')[1].split('},')[0] + '},', '')
    body = re.sub(r'\n\s*adsnote = \{[^}]*\},?', '', body)
    # fields the ADS export leaves empty (bibtex warns on a book without a publisher)
    for field, value in PATCH.get(key, {}).items():
        if not re.search(r'\n\s*%s\s*=' % field, body):
            body = body.rstrip(',') + ',\n    %s = {%s}' % (field, value)
    out.append('@%s{%s,%s\n}\n' % (kind, key, body.rstrip(',')))
open(BIB, 'w', encoding='utf-8').write('\n'.join(out))
print('wrote', len(order) - len(missing), 'entries to', os.path.relpath(BIB, ROOT))
