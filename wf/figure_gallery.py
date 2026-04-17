from __future__ import annotations

import html
import re

from pathlib import Path
from typing import Sequence


_PAGINATED_IMAGE_RE = re.compile(r"^(?P<stem>.+)_(?P<page>\d{3})\.png$")


def _title_from_stem(stem: str) -> str:
    return stem.replace("_", " ").strip().title()


def write_html_gallery(
    html_path: Path,
    title: str,
    image_paths: Sequence[Path],
) -> None:
    if len(image_paths) == 0:
        return

    html_path.parent.mkdir(parents=True, exist_ok=True)
    relative_paths = [path.relative_to(html_path.parent) for path in image_paths]

    blocks = []
    for idx, rel_path in enumerate(relative_paths, start=1):
        caption = f"Page {idx}"
        rel_str = html.escape(rel_path.as_posix(), quote=True)
        blocks.append(
            "<section class=\"page\">"
            f"<h2>{caption}</h2>"
            f"<img loading=\"lazy\" src=\"{rel_str}\" alt=\"{caption}\" />"
            "</section>"
        )

    body = "\n".join(blocks)
    page_title = html.escape(title)
    html_doc = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{page_title}</title>
  <style>
    body {{
      margin: 24px auto;
      max-width: 1600px;
      padding: 0 16px 32px;
      background: #f4f4f4;
      color: #1f1f1f;
      font-family: Arial, sans-serif;
    }}
    h1 {{
      margin: 0 0 8px;
    }}
    p {{
      margin: 0 0 20px;
      color: #404040;
    }}
    .page {{
      margin: 0 0 24px;
      background: #fff;
      border: 1px solid #ddd;
      border-radius: 8px;
      padding: 12px;
    }}
    .page h2 {{
      margin: 0 0 10px;
      font-size: 18px;
    }}
    img {{
      width: 100%;
      height: auto;
      display: block;
    }}
  </style>
</head>
<body>
  <h1>{page_title}</h1>
  <p>{len(image_paths)} pages</p>
  {body}
</body>
</html>
"""
    html_path.write_text(html_doc, encoding="utf-8")


def build_figure_site(results_dir: Path, figures_dir: Path) -> None:
    png_files = sorted(figures_dir.glob("*.png"))
    pdf_files = sorted(figures_dir.glob("*.pdf"))

    paginated_groups: dict[str, list[tuple[int, Path]]] = {}
    single_images: list[Path] = []

    for path in png_files:
        match = _PAGINATED_IMAGE_RE.match(path.name)
        if match is None:
            single_images.append(path)
            continue

        stem = match.group("stem")
        page = int(match.group("page"))
        paginated_groups.setdefault(stem, []).append((page, path))

    gallery_entries: list[tuple[str, Path, int]] = []
    for stem, pages in sorted(paginated_groups.items()):
        image_paths = [path for _, path in sorted(pages)]
        gallery_path = results_dir / f"{stem}.html"
        write_html_gallery(gallery_path, _title_from_stem(stem), image_paths)
        gallery_entries.append((stem, gallery_path, len(image_paths)))

    gallery_items = []
    for stem, gallery_path, page_count in gallery_entries:
        href = html.escape(gallery_path.relative_to(results_dir).as_posix(), quote=True)
        gallery_items.append(
            "<li>"
            f"<a href=\"{href}\">{html.escape(_title_from_stem(stem))}</a>"
            f" ({page_count} pages)"
            "</li>"
        )

    single_items = []
    for image_path in single_images:
        rel_path = image_path.relative_to(results_dir).as_posix()
        title = _title_from_stem(image_path.stem)
        rel_href = html.escape(rel_path, quote=True)
        alt = html.escape(title, quote=True)
        title_html = html.escape(title)
        single_items.append(
            "<section class=\"single-figure\">"
            f"<h2>{title_html}</h2>"
            f"<a href=\"{rel_href}\"><img loading=\"lazy\" src=\"{rel_href}\" alt=\"{alt}\" /></a>"
            "</section>"
        )

    legacy_pdf_items = []
    for pdf_path in pdf_files:
        rel_path = html.escape(pdf_path.relative_to(results_dir).as_posix(), quote=True)
        legacy_pdf_items.append(
            f"<li><a href=\"{rel_path}\">{html.escape(pdf_path.name)}</a></li>"
        )

    index_html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Figure Gallery</title>
  <style>
    body {{
      margin: 24px auto;
      max-width: 1600px;
      padding: 0 16px 32px;
      background: #f7f7f7;
      color: #1f1f1f;
      font-family: Arial, sans-serif;
    }}
    h1, h2 {{
      margin-bottom: 12px;
    }}
    .section {{
      margin-bottom: 32px;
      background: #fff;
      border: 1px solid #ddd;
      border-radius: 8px;
      padding: 16px;
    }}
    ul {{
      margin: 0;
      padding-left: 20px;
    }}
    .single-figure {{
      margin-bottom: 24px;
    }}
    .single-figure img {{
      width: 100%;
      height: auto;
      display: block;
      border: 1px solid #ddd;
      border-radius: 8px;
      background: #fff;
    }}
  </style>
</head>
<body>
  <h1>Figure Gallery</h1>
  <div class="section">
    <h2>Multi-Page Galleries</h2>
    <ul>
      {"".join(gallery_items) if gallery_items else "<li>No paginated figure galleries generated.</li>"}
    </ul>
  </div>
  <div class="section">
    <h2>Single Images</h2>
    {"".join(single_items) if single_items else "<p>No single-image figures generated.</p>"}
  </div>
  <div class="section">
    <h2>Legacy PDFs</h2>
    <ul>
      {"".join(legacy_pdf_items) if legacy_pdf_items else "<li>No PDF figures generated.</li>"}
    </ul>
  </div>
</body>
</html>
"""
    (results_dir / "figures.html").write_text(index_html, encoding="utf-8")
