#!/usr/bin/env python3
"""Assemble one version of the Alamo GitHub Pages site.

The documentation generator owns docs/build/html. This script deliberately
only knows about its published outputs, so it works with both the current
Sphinx build and the input-builder documentation on solidsgroup/alamo's
development branch.
"""

import argparse
import html
import re
import shutil
from pathlib import Path


SLUG_RE = re.compile(r"^(master|development|pr-[0-9]+)$")


def copy_tree(source, destination, *, required=False):
    source = Path(source)
    destination = Path(destination)
    if not source.is_dir():
        if required:
            raise FileNotFoundError(f"Required Pages input does not exist: {source}")
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source, destination, dirs_exist_ok=True)
    return True


def write_report_index(reports_root):
    reports_root = Path(reports_root)
    entries = []
    for directory in sorted(path for path in reports_root.iterdir()
                            if path.is_dir()):
        entry = directory / "report.html"
        if not entry.is_file():
            entry = directory / "report" / "report.html"
        if entry.is_file():
            entries.append((directory.name, entry.relative_to(reports_root)))
    if not entries:
        return
    links = "\n".join(
        f'<li><a href="{html.escape(path.as_posix())}">'
        f"{html.escape(name)}</a></li>"
        for name, path in entries
    )
    (reports_root / "index.html").write_text(
        "<!doctype html>\n"
        '<meta charset="utf-8">\n'
        "<title>Alamo CI reports</title>\n"
        "<h1>Alamo CI reports</h1>\n"
        f"<ul>{links}</ul>\n",
        encoding="utf-8",
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--slug", required=True,
                        help="master, development, or pr-N")
    parser.add_argument("--output", default="githubpages")
    parser.add_argument("--docs", default="docs/build/html")
    parser.add_argument("--coverage", default="cov")
    parser.add_argument("--reports", default=None,
                        help="Directory containing named report directories")
    parser.add_argument("--input-builders",
                        default="docs/source/_static/input-builders")
    parser.add_argument("--input-schemas",
                        default="docs/source/_static/input-schemas")
    args = parser.parse_args()

    if not SLUG_RE.fullmatch(args.slug):
        parser.error("--slug must be master, development, or pr-N")

    output = Path(args.output)
    if output.exists():
        shutil.rmtree(output)

    version_root = output / "docs" / args.slug
    copy_tree(args.docs, version_root, required=True)
    copy_tree(args.coverage, version_root / "cov")

    if args.reports:
        reports_root = version_root / "reports"
        if copy_tree(args.reports, reports_root):
            write_report_index(reports_root)

    # These directories are introduced by the input-builder documentation
    # work on the upstream development branch. They are optional here and
    # become part of the bundle automatically after that work is merged.
    input_builders = Path(args.input_builders)
    input_schemas = Path(args.input_schemas)
    if not input_builders.is_dir():
        input_builders = Path(args.docs) / "_static" / "input-builders"
    if not input_schemas.is_dir():
        input_schemas = Path(args.docs) / "_static" / "input-schemas"
    copy_tree(input_builders, output / "inputs" / args.slug)
    copy_tree(input_schemas, output / "input-schemas" / args.slug)

    # Development is the continuously updated default site. Master and PR
    # prefix uploads deliberately do not replace the root redirect.
    if args.slug == "development":
        output.mkdir(parents=True, exist_ok=True)
        destination = f"docs/{args.slug}/"
        (output / "index.html").write_text(
            "<!doctype html>\n"
            '<meta charset="utf-8">\n'
            f'<meta http-equiv="refresh" content="0;url={html.escape(destination)}">\n'
            "<title>Alamo documentation</title>\n"
            f'<a href="{html.escape(destination)}">Open Alamo documentation</a>\n',
            encoding="utf-8",
        )

    print(f"Packaged Pages content for {args.slug} in {output}")


if __name__ == "__main__":
    main()
