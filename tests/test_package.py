from __future__ import annotations

import importlib.metadata

import my_utils as m


def test_version():
    assert importlib.metadata.version("my_utils") == m.__version__
