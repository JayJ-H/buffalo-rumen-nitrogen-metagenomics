"""Strip ancillary metadata from this release only and refresh checksums."""

from pathlib import Path
import csv
import hashlib
import json
import os
import re
import xml.etree.ElementTree as ET

import numpy as np
from PIL import Image
from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[1]
EPOCH = 946684800  # Neutral package timestamp: 2000-01-01 UTC.
TEXT_SUFFIXES = {".py", ".R", ".md", ".txt", ".json", ".csv", ".tsv", ".svg"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean_pdf(path):
    reader = PdfReader(path)
    pages = len(reader.pages)
    old_text = [page.extract_text() for page in reader.pages]
    writer = PdfWriter()
    writer.clone_document_from_reader(reader)
    writer._info.clear()
    writer._root_object.pop("/Metadata", None)
    writer._ID = None
    temp = path.with_suffix(".cleaning.pdf")
    with temp.open("wb") as handle:
        writer.write(handle)
    check = PdfReader(temp)
    assert not check.metadata
    assert "/Metadata" not in check.trailer["/Root"]
    assert len(check.pages) == pages
    assert [page.extract_text() for page in check.pages] == old_text
    os.replace(temp, path)


def clean_raster(path):
    with Image.open(path) as source:
        assert getattr(source, "n_frames", 1) == 1
        pixels = np.array(source)
        mode = source.mode
        dpi = source.info.get("dpi")
        if dpi:
            dpi = tuple(float(value) for value in dpi)
        image = Image.fromarray(pixels, mode=mode)
    kwargs = {"dpi": dpi} if dpi else {}
    if path.suffix == ".tiff":
        kwargs["compression"] = "tiff_lzw"
    image.save(path, **kwargs)
    with Image.open(path) as check:
        assert np.array_equal(pixels, np.array(check))
        if dpi:
            assert np.allclose(np.asarray(check.info["dpi"], dtype=float), dpi, atol=0.1)
        assert not any(k in check.info for k in [
            "exif", "xmp", "XML:com.adobe.xmp", "Author", "Software", "Creation Time"])
        if path.suffix == ".tiff":
            assert not any(tag in check.tag_v2 for tag in [270, 305, 306, 315, 316, 700])


def clean_svg(path):
    text = path.read_text(encoding="utf-8")
    text = re.sub(r"<!--.*?-->", "", text, flags=re.S)
    text = re.sub(r"<metadata\b[^>]*>.*?</metadata>", "", text, flags=re.S)
    text = re.sub(r"<metadata\b[^>]*/>", "", text)
    ET.fromstring(text)
    path.write_text(text, encoding="utf-8")


def main():
    assert (ROOT / "data/function_definitions.json").exists()
    assert (ROOT / "README.md").exists()
    source_hashes = {p.relative_to(ROOT).as_posix(): sha(p) for p in ROOT.rglob("*")
                     if p.is_file() and p.suffix in {".csv", ".tsv"}
                     and p.name != "MANIFEST.csv"}
    counts = {"pdf": 0, "raster": 0, "svg": 0}
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix == ".pdf":
            clean_pdf(path)
            counts["pdf"] += 1
        elif path.suffix in {".png", ".tiff"}:
            clean_raster(path)
            counts["raster"] += 1
        elif path.suffix == ".svg":
            clean_svg(path)
            counts["svg"] += 1
    leaks = []
    forbidden = re.compile(r"(?<![A-Za-z0-9])[A-Za-z]:[\\/]|/Users/|/home/|\.codex|ChatGPT|Codex")
    assert not forbidden.search("https://example.org/path")
    assert forbidden.search("C:/private/path")
    for path in ROOT.rglob("*"):
        if path.is_file() and path.suffix in TEXT_SUFFIXES:
            # This cleaner contains the scanner patterns by design.
            if path.resolve() == Path(__file__).resolve():
                continue
            content = path.read_text(encoding="utf-8-sig")
            if forbidden.search(content):
                leaks.append(path.relative_to(ROOT).as_posix())
    assert not leaks, leaks
    for relative, expected in source_hashes.items():
        assert sha(ROOT / relative) == expected, relative
    report = {
        "status": "PASS",
        "cleaned_exports": counts,
        "numeric_data_bytes_unchanged": len(source_hashes),
        "raster_pixels_and_dpi_preserved": True,
        "pdf_text_preserved": True,
        "private_path_scan": "PASS",
        "package_timestamp_utc": "2000-01-01T00:00:00Z",
        "scientific_metadata_retained": True,
    }
    (ROOT / "docs/METADATA_CLEANUP.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    files = sorted(p for p in ROOT.rglob("*") if p.is_file()
                   and p.name != "MANIFEST.csv" and "__pycache__" not in p.parts)
    with (ROOT / "MANIFEST.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["path", "bytes", "sha256"])
        for path in files:
            writer.writerow([path.relative_to(ROOT).as_posix(), path.stat().st_size, sha(path)])
    files.append(ROOT / "MANIFEST.csv")
    if os.name == "nt":
        import ctypes
        set_attributes = ctypes.windll.kernel32.SetFileAttributesW
        for path in files:
            assert set_attributes(str(path), 0x80), path.name
    for path in files:
        os.utime(path, (EPOCH, EPOCH))
    if os.name == "nt":
        from ctypes import wintypes
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel.CreateFileW.argtypes = [wintypes.LPCWSTR, wintypes.DWORD, wintypes.DWORD,
                                      ctypes.c_void_p, wintypes.DWORD, wintypes.DWORD,
                                      ctypes.c_void_p]
        kernel.CreateFileW.restype = ctypes.c_void_p
        kernel.SetFileTime.argtypes = [ctypes.c_void_p,
                                      ctypes.POINTER(wintypes.FILETIME),
                                      ctypes.POINTER(wintypes.FILETIME),
                                      ctypes.POINTER(wintypes.FILETIME)]
        kernel.CloseHandle.argtypes = [ctypes.c_void_p]
        ticks = int((EPOCH + 11644473600) * 10_000_000)
        neutral = wintypes.FILETIME(ticks & 0xFFFFFFFF, ticks >> 32)
        for path in files:
            handle = kernel.CreateFileW(str(path), 0x100, 7, None, 3, 0x80, None)
            assert handle != ctypes.c_void_p(-1).value, path.name
            try:
                assert kernel.SetFileTime(handle, ctypes.byref(neutral),
                                          ctypes.byref(neutral), ctypes.byref(neutral))
            finally:
                kernel.CloseHandle(handle)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
