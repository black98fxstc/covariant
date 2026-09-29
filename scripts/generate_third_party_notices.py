#!/usr/bin/env python3
"""Collect vcpkg package copyright notices for Leonard release artifacts."""

import argparse
from pathlib import Path


def find_copyright_file(package_dir: Path):
    for candidate in package_dir.iterdir():
        if candidate.name.lower() == "copyright" and candidate.is_file():
            return candidate
    return None


def collect_notices(vcpkg_root: Path):
    triplet_roots = []
    if (vcpkg_root / "share").is_dir():
        triplet_roots.append(vcpkg_root)
    else:
        triplet_roots.extend(
            child for child in vcpkg_root.iterdir()
            if child.is_dir() and (child / "share").is_dir()
        )

    notices = []
    missing = []
    for triplet_root in sorted(triplet_roots, key=lambda path: path.name.lower()):
        share_root = triplet_root / "share"
        for package_dir in sorted(share_root.iterdir(), key=lambda path: path.name.lower()):
            if not package_dir.is_dir():
                continue
            copyright_file = find_copyright_file(package_dir)
            label = f"{triplet_root.name}/{package_dir.name}"
            if copyright_file is None:
                missing.append(label)
                continue
            notices.append((label, copyright_file))

    return notices, missing


def write_notices(output_path: Path, notices, missing):
    lines = [
        "Third-party software notices for Leonard",
        "",
        "This file was generated from the vcpkg packages used to build this release.",
        "",
    ]
    for label, copyright_file in notices:
        lines.extend(
            [
                "=" * 72,
                label,
                "=" * 72,
                copyright_file.read_text(encoding="utf-8", errors="replace").rstrip(),
                "",
            ]
        )

    if missing:
        lines.extend(
            [
                "=" * 72,
                "Packages without a vcpkg copyright file",
                "=" * 72,
                "These entries are commonly vcpkg aliases or build-tool packages; review them before release:",
                *missing,
                "",
            ]
        )

    output_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcpkg-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if not args.vcpkg_root.is_dir():
        parser.error(f"vcpkg installation not found: {args.vcpkg_root}")

    notices, missing = collect_notices(args.vcpkg_root)
    if not notices:
        parser.error(f"no vcpkg copyright files found under {args.vcpkg_root}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_notices(args.output, notices, missing)
    print(f"Collected {len(notices)} package notices into {args.output}")
    if missing:
        print(f"Review required for {len(missing)} packages without copyright metadata:")
        for package in missing:
            print(f"  {package}")


if __name__ == "__main__":
    main()
