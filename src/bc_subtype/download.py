"""Download a pinned public DataHub snapshot with fail-closed integrity checks."""

import hashlib
import json
from importlib.resources import files
from pathlib import Path
from urllib.request import Request, urlopen


def source_manifest():
    return json.loads(files("bc_subtype").joinpath("resources/sources.json").read_text())


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_raw(raw_dir):
    manifest = source_manifest()
    for name, expected in manifest["files"].items():
        path = Path(raw_dir) / name
        if not path.is_file():
            raise FileNotFoundError(f"Missing {path}; run bc-subtype download first")
        if path.stat().st_size != expected["bytes"] or sha256(path) != expected["sha256"]:
            raise ValueError(f"Integrity check failed: {path}; move it aside and download again")
    return manifest


def download(raw_dir):
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    manifest = source_manifest()
    base = (
        f"https://media.githubusercontent.com/media/{manifest['repository']}/"
        f"{manifest['revision']}/public/{manifest['study_id']}"
    )
    for name, expected in manifest["files"].items():
        target = raw_dir / name
        if target.exists():
            if target.stat().st_size == expected["bytes"] and sha256(target) == expected["sha256"]:
                print(f"Verified cached {name}", flush=True)
                continue
            raise ValueError(f"Existing file failed integrity check: {target}; move it aside")
        temporary = target.with_suffix(target.suffix + ".part")
        request = Request(
            f"{base}/{name}?download=true", headers={"User-Agent": "bc-subtype-classifier/0.1"}
        )
        print(f"Downloading {name} ({expected['bytes']:,} bytes)", flush=True)
        with urlopen(request, timeout=60) as response, temporary.open("wb") as stream:
            while block := response.read(1024 * 1024):
                stream.write(block)
        if temporary.stat().st_size != expected["bytes"] or sha256(temporary) != expected["sha256"]:
            raise ValueError(f"Downloaded file failed integrity check: {temporary}")
        temporary.replace(target)
    verify_raw(raw_dir)
    return manifest
