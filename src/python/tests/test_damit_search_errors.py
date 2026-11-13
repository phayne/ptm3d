#!/usr/bin/env python3
"""A DAMIT search that cannot reach the site reports the failure, not "no model".

The search walks DAMIT's paginated browse list. When a page answered with a
server error (HTTP 500 or 504), the page was treated as empty, and a search
for an asteroid that DAMIT does hold returned an empty list. These checks
serve the pages from memory. No network.

Run:  python3 src/python/tests/test_damit_search_errors.py
"""
from __future__ import annotations

import io
import sys
import urllib.error
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))

PAGE_1 = ('<a href="/projects/damit/asteroids/view/101">(1) Ceres</a>'
          '<a href="/projects/damit/asteroids/view/102">(2) Pallas</a>')


def _client(status):
    """A client whose browse pages answer with `status` (200 serves PAGE_1)."""
    from damit.damit_client import DAMITClient

    class _Client(DAMITClient):
        def _fetch(self, url):
            if status == 200:
                return PAGE_1
            err = urllib.error.HTTPError(url, status, "server says no", {}, io.BytesIO())
            raise RuntimeError(f"Failed to fetch {url}: {err}") from err

    return _Client()


def test_a_server_error_is_reported():
    for status in (500, 504):
        try:
            result = _client(status).search_asteroid("433")
        except RuntimeError as exc:
            assert "DAMIT did not answer" in str(exc) and str(status) in str(exc), exc
        else:
            raise AssertionError(f"HTTP {status} returned {result!r} instead of an error")
        try:
            _client(status).get_model_by_number(433)
        except RuntimeError:
            pass
        else:
            raise AssertionError(f"get_model_by_number hid HTTP {status}")
    print("    [search] HTTP 500 and 504 raise RuntimeError naming the status")


def test_a_missing_page_is_still_the_end_of_the_list():
    assert _client(404)._browse_page(9999) == []
    assert _client(200).search_asteroid("2") == [
        {"damit_id": 102, "number": 2, "name": "Pallas"}]
    print("    [search] HTTP 404 is an empty page; a page that answers is parsed as before")


CHECKS = [v for k, v in sorted(globals().items()) if k.startswith("test_")]


def main():
    failed = 0
    for check in CHECKS:
        try:
            check()
        except Exception as exc:                             # noqa: BLE001
            failed += 1
            print(f"FAIL {check.__name__}: {exc!r}")
        else:
            print(f"PASS {check.__name__}")
    print(f"\n{len(CHECKS) - failed}/{len(CHECKS)} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
