import os, re
from pypdf import PdfReader, PdfWriter, PageObject, Transformation
from pypdf.generic import RectangleObject

input_dir  = "."
out_kaon   = "KaonExclComp.pdf"
out_proton = "ProtonExclComp.pdf"
peak_path = "PeakAreas_Kaon_Proton.pdf"
pattern    = re.compile(r"nSigmaTPC_e_(KaonExclComp|ProtonExclComp)-(\d+\.\d+)\.pdf$")

def collect_files():
    kaon, proton = [], []
    for fn in sorted(os.listdir(input_dir)):
        m = pattern.match(fn)
        if not m:
            continue
        comp, sigma = m.group(1), float(m.group(2))
        path = os.path.join(input_dir, fn)
        (kaon if comp == "KaonExclComp" else proton).append((sigma, path))
    kaon.sort(key=lambda x: x[0])
    proton.sort(key=lambda x: x[0])
    return kaon, proton

def smallest_defined_box(page):
    boxes = []
    for name in ("artbox", "trimbox", "bleedbox", "cropbox", "mediabox"):
        box = getattr(page, name, None)
        if box is None:
            continue
        l, b, r, t = float(box.left), float(box.bottom), float(box.right), float(box.top)
        boxes.append((l, b, r, t, (r-l)*(t-b)))
    if not boxes:
        mb = page.mediabox
        return float(mb.left), float(mb.bottom), float(mb.right), float(mb.top)
    l, b, r, t, _ = min(boxes, key=lambda x: x[4])
    return l, b, r, t

def render_page_to_size(src, W, H, margin_frac=0.02):
    l, b, r, t = smallest_defined_box(src)
    w, h = r - l, t - b
    cx_src, cy_src = (l + r) / 2.0, (b + t) / 2.0

    effW = W * (1 - 2*margin_frac)
    effH = H * (1 - 2*margin_frac)
    s = min(effW / w, effH / h)  

    new = PageObject.create_blank_page(width=W, height=H)
    tf = (Transformation()
          .translate(-cx_src, -cy_src)
          .scale(s)
          .translate(W/2.0, H/2.0))
    new.merge_transformed_page(src, tf)
    return new

def merge(file_list, out_pdf, cols=2, rows=3, trim=(0,0,0,0),
          safe=0.997, prepend_pages=None, title_margin=0.02):
    writer = PdfWriter()

    # --- Eingabeseiten sammeln ---
    pages = []
    for _, path in file_list:
        pages.extend(PdfReader(path).pages)

    # --- prepend_pages in Liste verwandeln (None / single / list/tuple) ---
    if prepend_pages is None:
        prepend_pages = []
    elif not isinstance(prepend_pages, (list, tuple)):
        prepend_pages = [prepend_pages]
    else:
        prepend_pages = list(prepend_pages)

    # --- Seitenlayout bestimmen (W,H) ---
    if pages:
        l0, b0, r0, t0 = smallest_defined_box(pages[0])
        w0, h0 = r0 - l0, t0 - b0
        tl, tr, tb, tt = trim
        l0 += tl*w0; r0 -= tr*w0; b0 += tb*h0; t0 -= tt*h0
        cw0, ch0 = r0 - l0, t0 - b0
        cell_w, cell_h = ch0, cw0        # wegen 90°-Rotation
        W = cols * cell_w
        H = rows * cell_h
    elif prepend_pages:
        # Fallback: orientiere dich an der ersten Titelseite
        l1, b1, r1, t1 = smallest_defined_box(prepend_pages[0])
        W, H = (r1 - l1), (t1 - b1)
    else:
        # letzter Fallback: A4
        W, H = 595, 842

    # --- Titelseiten voranstellen ---
    for p in prepend_pages:
        writer.add_page(render_page_to_size(p, W, H, margin_frac=title_margin))

    # --- Wenn keine Mosaik-Seiten: direkt schreiben ---
    if not pages:
        with open(out_pdf, "wb") as f:
            writer.write(f)
        print(f"Erstellt: {out_pdf}")
        return

    # --- Mosaik-Seiten bauen (wie gehabt) ---
    tl, tr, tb, tt = trim
    cell_w = W / cols
    cell_h = H / rows

    for i in range(0, len(pages), cols*rows):
        new = PageObject.create_blank_page(width=W, height=H)
        chunk = pages[i:i + cols*rows]

        for idx, src in enumerate(chunk):
            l, b, r, t = smallest_defined_box(src)
            w, h = r - l, t - b
            l += tl*w; r -= tr*w; b += tb*h; t -= tt*h
            w, h = r - l, t - b

            src.cropbox = RectangleObject([l, b, r, t])

            cx_src, cy_src = (l + r)/2.0, (b + t)/2.0
            rw, rh = h, w  # wegen Rotation
            s = min(cell_w / rw, cell_h / rh) * safe

            col, row = idx % cols, idx // cols
            cx = (col + 0.5) * cell_w
            cy = H - (row + 0.5) * cell_h

            tf = (Transformation()
                  .translate(-cx_src, -cy_src)
                  .rotate(270)
                  .scale(s)
                  .translate(cx, cy))
            new.merge_transformed_page(src, tf)

        writer.add_page(new)

    with open(out_pdf, "wb") as f:
        writer.write(f)
    print(f"Erstellt: {out_pdf}")

peak_kaon_pages   = []
peak_proton_pages = []

if os.path.exists(peak_path):
    try:
        pr = PdfReader(peak_path)
        pages = pr.pages
        if len(pages) > 0: peak_kaon_pages.append(pages[0])   # Seite 1
        if len(pages) > 1: peak_proton_pages.append(pages[1]) # Seite 2
        if len(pages) > 2: peak_kaon_pages.append(pages[2])   # Seite 3
        if len(pages) > 3: peak_proton_pages.append(pages[3]) # Seite 4
    except Exception as e:
        print("Warnung: Peak PDF konnte nicht gelesen werden:", e)

kaon_files, proton_files = collect_files()
merge(kaon_files,   out_kaon,   trim=(0,0,0,0), safe=0.997, prepend_pages=peak_kaon_pages,   title_margin=0.02)
merge(proton_files, out_proton, trim=(0,0,0,0), safe=0.997, prepend_pages=peak_proton_pages, title_margin=0.02)
