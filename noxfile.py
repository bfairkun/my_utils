#!/usr/bin/env -S uv run --script
# ^ Makes this file directly executable: ./noxfile.py
# uv run --script reads the inline dependencies below and runs with them

# /// script
# dependencies = ["nox>=2025.2.9"]
# ///
# ^ Inline script metadata: uv uses this to know nox is required to run this file

"""Nox runner."""

# Nox is a task runner (like `make` but Python).
# Each @nox.session below is a task you can run with: nox -s <session_name>
# Running just `nox` runs all sessions that don't have default=False.

from __future__ import annotations  # Allows modern type hint syntax on Python 3.10

import shutil
from pathlib import Path

import nox

DIR = Path(__file__).parent.resolve()  # Absolute path to the repo root
PROJECT = nox.project.load_toml()  # Reads pyproject.toml (used to get dependency lists)

nox.needs_version = ">=2025.2.9"  # Minimum nox version required
nox.options.default_venv_backend = (
    "uv|virtualenv"  # Use uv for venvs if available, else virtualenv
)


# ── Session: lint ─────────────────────────────────────────────────────────────
# Runs pre-commit hooks via a tool called `prek` (a pre-commit wrapper).
# Usage: nox -s lint
# Note: We don't use this in CI — pre-commit runs locally via git hooks instead.
@nox.session
def lint(session: nox.Session) -> None:
    """
    Run the linter.
    """
    session.install("prek")  # Install prek into a temporary isolated venv
    session.run(
        "prek",
        "run",
        "--all-files",
        "--show-diff-on-failure",
        *session.posargs,
        # --all-files: check every file, not just staged ones
        # --show-diff-on-failure: print the diff when a hook fails
        # *session.posargs: pass any extra args you provide on the command line
    )


# ── Session: pylint ───────────────────────────────────────────────────────────
# Runs pylint (a thorough but slow Python linter).
# Usage: nox -s pylint
# Note: Removed from CI because it requires building pysam from source.
@nox.session
def pylint(session: nox.Session) -> None:
    """
    Run Pylint.
    """
    # This needs to be installed into the package environment, and is slower
    # than a pre-commit check
    session.install(
        "-e.", "pylint>=3.2"
    )  # Install this package + pylint into a temp venv
    session.run(
        "pylint", "my_utils", *session.posargs
    )  # Run pylint on the my_utils package


# ── Session: tests ────────────────────────────────────────────────────────────
# Runs pytest. This is what CI runs (but CI calls `uv run pytest` directly).
# Usage: nox -s tests
@nox.session
def tests(session: nox.Session) -> None:
    """
    Run the unit and regular tests.
    """
    test_deps = nox.project.dependency_groups(
        PROJECT, "test"
    )  # Reads [dependency-groups] test from pyproject.toml
    session.install("-e.", *test_deps)  # Install this package + pytest + pytest-cov
    session.run("pytest", *session.posargs)  # Run pytest (pass extra args if provided)


# ── Session: docs ─────────────────────────────────────────────────────────────
# Builds or serves the documentation website (if you set one up with mkdocs).
# default=False means this does NOT run with plain `nox` — must explicitly call: nox -s docs
# reuse_venv=True means it reuses the venv between runs (faster for iterating on docs)
@nox.session(reuse_venv=True, default=False)
def docs(session: nox.Session) -> None:
    """
    Make or serve the docs. Pass --non-interactive to avoid serving.
    """

    doc_deps = nox.project.dependency_groups(
        PROJECT, "docs"
    )  # Reads [dependency-groups] docs from pyproject.toml
    session.install("-e.", *doc_deps)  # Install this package + mkdocs + plugins

    if session.interactive:
        session.run(
            "mkdocs", "serve", "--clean", *session.posargs
        )  # Live-reload server at localhost:8000
    else:
        session.run(
            "mkdocs", "build", "--clean", *session.posargs
        )  # Build static HTML to site/ directory


# ── Session: build ────────────────────────────────────────────────────────────
# Builds the distributable package files (wheel + sdist).
# default=False means it does NOT run with plain `nox`.
# Usage: nox -s build
@nox.session(default=False)
def build(session: nox.Session) -> None:
    """
    Build an SDist and wheel.
    """

    build_path = DIR.joinpath("build")
    if build_path.exists():
        shutil.rmtree(build_path)  # Clean old build artifacts first

    session.install(
        "build"
    )  # Install the `build` package (Python's standard build frontend)
    session.run("python", "-m", "build")  # Creates dist/*.whl and dist/*.tar.gz


if __name__ == "__main__":
    nox.main()  # Allows running this file directly: ./noxfile.py
