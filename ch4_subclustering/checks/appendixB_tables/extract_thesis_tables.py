#!/usr/bin/env python3
"""Extract Appendix B tables (B2-B95) from the thesis docx into one CSV per table.

In the thesis docx each table's caption is the paragraph that FOLLOWS the table.
Usage: extract_thesis_tables.py <thesis.docx> <out_dir>
Writes <out_dir>/B<nn>.csv plus <out_dir>/index.csv (table, caption, n_rows).
"""
import csv, os, re, sys
import docx
from docx.table import Table
from docx.text.paragraph import Paragraph

docx_path, out_dir = sys.argv[1], sys.argv[2]
os.makedirs(out_dir, exist_ok=True)
d = docx.Document(docx_path)

items = []
for el in d.element.body.iterchildren():
    tag = el.tag.split('}')[1]
    if tag == 'p':
        t = Paragraph(el, d).text.strip()
        if t:
            items.append(('p', t))
    elif tag == 'tbl':
        items.append(('t', Table(el, d)))

index = []
for i, (ty, x) in enumerate(items):
    if ty != 't' or i + 1 >= len(items) or items[i + 1][0] != 'p':
        continue
    cap = items[i + 1][1]
    m = re.match(r'Table\s*B(\d+)\.', cap)
    if not m or int(m.group(1)) < 2:
        continue
    nn = int(m.group(1))
    rows = [[c.text.strip() for c in r.cells] for r in x.rows]
    with open(os.path.join(out_dir, f'B{nn:02d}.csv'), 'w', newline='') as fh:
        csv.writer(fh).writerows(rows)
    index.append((f'B{nn:02d}', cap, len(rows) - 1))

with open(os.path.join(out_dir, 'index.csv'), 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['table', 'caption', 'n_data_rows'])
    w.writerows(index)
print(f'extracted {len(index)} tables to {out_dir}')
