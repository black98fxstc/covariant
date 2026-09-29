#!/usr/bin/env python3
"""Validate Leonard's release version and derive platform metadata."""

import argparse
import json
import os
import re
from pathlib import Path


VERSION_PATTERN = re.compile(
    r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)"
    r"(?:-(alpha|beta|rc)\.([1-9]\d*))?$"
)
MAC_SUFFIX = {"alpha": "a", "beta": "b", "rc": "fc"}


def parse_version(value):
    match = VERSION_PATTERN.fullmatch(value)
    if not match:
        raise ValueError(
            f"invalid version {value!r}; expected MAJOR.MINOR.PATCH"
            " or MAJOR.MINOR.PATCH-(alpha|beta|rc).NUMBER"
        )

    major, minor, patch, prerelease, prerelease_number = match.groups()
    if len(major) > 4 or len(minor) > 2 or len(patch) > 2:
        raise ValueError(f"version {value!r} exceeds macOS bundle component limits")
    if prerelease_number and int(prerelease_number) > 255:
        raise ValueError(f"version {value!r} exceeds the macOS prerelease limit of 255")

    short_version = f"{major}.{minor}.{patch}"
    bundle_version = short_version
    if prerelease:
        bundle_version += f"{MAC_SUFFIX[prerelease]}{prerelease_number}"
    return {
        "version": value,
        "short_version": short_version,
        "mac_bundle_version": bundle_version,
    }


def read_version(path):
    with path.open(encoding="utf-8") as source:
        manifest = json.load(source)
    value = manifest.get("version-string")
    if not isinstance(value, str):
        raise ValueError(f"{path} must contain a string version-string")
    return parse_version(value)


def validate_tag(tag, version):
    if not tag.startswith("v"):
        raise ValueError(f"invalid release tag {tag!r}; expected v{version}")
    tag_version = parse_version(tag[1:])["version"]
    if tag_version != version:
        raise ValueError(
            f"release tag {tag!r} does not match vcpkg.json ({version}); "
            "update version-string before tagging"
        )


def main():
    repository = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--version-file", type=Path, default=repository / "vcpkg.json")
    parser.add_argument("--tag", help="release tag to validate (for example, v0.1.0-alpha.3)")
    parser.add_argument(
        "--print",
        dest="field",
        choices=("version", "short-version", "mac-bundle-version"),
        default="version",
    )
    parser.add_argument(
        "--github-output",
        action="store_true",
        help="write all version fields to the path in GITHUB_OUTPUT",
    )
    args = parser.parse_args()

    try:
        metadata = read_version(args.version_file)
        tag = args.tag
        if tag is None and os.environ.get("GITHUB_REF_TYPE") == "tag":
            tag = os.environ.get("GITHUB_REF_NAME")
            if not tag:
                parser.error("GITHUB_REF_NAME is required for a tag build")
        if tag:
            validate_tag(tag, metadata["version"])
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.error(str(error))

    if args.github_output:
        output_path = os.environ.get("GITHUB_OUTPUT")
        if not output_path:
            parser.error("GITHUB_OUTPUT is required with --github-output")
        with open(output_path, "a", encoding="utf-8") as output:
            for key, value in metadata.items():
                output.write(f"{key}={value}\n")
    else:
        print(metadata[args.field.replace("-", "_")])


if __name__ == "__main__":
    main()
