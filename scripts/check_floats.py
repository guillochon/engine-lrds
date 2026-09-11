"""Report, for each figure and table, the page of its first mention in the body and the
page on which its caption appears (from the compiled PDF), and flag any float that
lands more than one page after its first reference or shares a page with another figure.

    python scripts/check_floats.py [paper/lrd_engine.pdf]
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PDF = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'paper', 'lrd_engine.pdf')


def pages(pdf):
    txt = subprocess.run(['pdftotext', pdf, '-'], capture_output=True, text=True,
                         encoding='utf-8', errors='replace').stdout
    return txt.split('\f')


def main():
    pg = pages(PDF)
    ok = True
    fig_pages = {}
    # tables deliberately placed in the appendix, after the reference list, are exempt
    ref_list_p = next((i + 1 for i, t in enumerate(pg) if 'REFERENCES' in t), len(pg))
    for kind in ['Figure', 'Table']:
        n = 1
        while True:
            cap = re.compile(r'^\s*%s %d\. ' % (kind, n), re.M)
            # a body mention: not at the start of a line (that is the caption) and not "Figure 21"
            ref = re.compile(r'(?<!\n)%s %d(?![0-9])' % (kind, n))
            cap_p = next((i + 1 for i, t in enumerate(pg) if cap.search(t)), None)
            if cap_p is None:
                break

            # a mention that is not the caption itself and not another paper's figure ("their Figure 3", "Figure 3 of")
            def mentions(t):
                return [m for m in ref.finditer(t)
                        if not t[max(0, m.start() - 6):m.start()].endswith('their ')
                        and not re.match(r'\s+of\s', t[m.end():m.end() + 6])]
            ref_p = next((i + 1 for i, t in enumerate(pg) if mentions(t)), None)
            flag = ''
            if kind == 'Table' and cap_p > ref_list_p:
                flag = '  (appendix)'
            elif ref_p is not None and cap_p > ref_p + 1:
                flag = '  <-- %d pages after first reference' % (cap_p - ref_p)
                ok = False
            elif ref_p is not None and cap_p < ref_p:
                flag = '  <-- appears before first reference'
                ok = False
            print('%s %d: first reference p.%s, caption p.%d%s' % (kind, n, ref_p, cap_p, flag))
            if kind == 'Figure':
                fig_pages.setdefault(cap_p, []).append(n)
            n += 1
    for p, fs in fig_pages.items():
        if len(fs) > 1:
            print('two figures on page %d: %s' % (p, fs))
            ok = False
    print('OK' if ok else 'PLACEMENT PROBLEMS')


if __name__ == '__main__':
    main()
