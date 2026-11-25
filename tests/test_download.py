import hashlib
import io

import pytest

from bc_subtype import download as module


def test_download_verifies_bytes_and_rejects_corrupt_cache(tmp_path, monkeypatch):
    contents = b"real-data-not-an-lfs-pointer\n"
    manifest = {
        "repository": "example/data",
        "revision": "pinned",
        "study_id": "study",
        "files": {
            "file.txt": {
                "bytes": len(contents),
                "sha256": hashlib.sha256(contents).hexdigest(),
            }
        },
    }
    monkeypatch.setattr(module, "source_manifest", lambda: manifest)
    monkeypatch.setattr(module, "urlopen", lambda *a, **k: io.BytesIO(contents))
    module.download(tmp_path)
    assert (tmp_path / "file.txt").read_bytes() == contents
    (tmp_path / "file.txt").write_bytes(b"corruption")
    with pytest.raises(ValueError, match="integrity check"):
        module.download(tmp_path)


def test_corrupt_transfer_never_becomes_final_source_file(tmp_path, monkeypatch):
    manifest = {
        "repository": "example/data",
        "revision": "pinned",
        "study_id": "study",
        "files": {"file.txt": {"bytes": 100, "sha256": "0" * 64}},
    }
    monkeypatch.setattr(module, "source_manifest", lambda: manifest)
    monkeypatch.setattr(module, "urlopen", lambda *a, **k: io.BytesIO(b"bad"))
    with pytest.raises(ValueError, match="integrity check"):
        module.download(tmp_path)
    assert not (tmp_path / "file.txt").exists()
