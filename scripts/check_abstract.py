"""Check the abstract of paper/lrd_engine.tex against the arXiv 1920-character
limit and the no-custom-macro rule (.cursor/rules/abstract-*.mdc)."""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEX = os.path.join(ROOT, 'paper', 'lrd_engine.tex')

t = open(TEX, encoding='utf-8').read()
# every \newcommand in the preamble is forbidden in the abstract, plus the cite commands
MACROS = set(re.findall(r'\\newcommand\{\\([A-Za-z]+)\}', t)) | {'citep', 'citet', 'citetalias', 'citealt'}

a = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', t, re.S).group(1)
s = ' '.join(a.split())
found = sorted({m for m in re.findall(r'\\([A-Za-z]+)', s) if m in MACROS})
leads = s.lstrip().startswith(('\\cite', 'Guillochon'))
print('abstract characters (collapsed, incl. TeX):', len(s))
print('forbidden macros:', found or 'none')
print('leads with citation:', leads)
sys.exit(0 if len(s) <= 1920 and not found and not leads else 1)
