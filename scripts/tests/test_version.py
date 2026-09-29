import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[2]
SCRIPT = REPOSITORY / "scripts" / "version.py"
SPEC = importlib.util.spec_from_file_location("leonard_version", SCRIPT)
VERSION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERSION)


class VersionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metadata = VERSION.read_version(REPOSITORY / "vcpkg.json")

    def test_prerelease_metadata(self):
        self.assertEqual(
            VERSION.parse_version("0.1.0-alpha.2"),
            {
                "version": "0.1.0-alpha.2",
                "short_version": "0.1.0",
                "mac_bundle_version": "0.1.0a2",
            },
        )

    def test_invalid_versions_are_rejected(self):
        for value in ("v0.1.0", "0.1", "01.1.0", "0.1.0-preview.1", "0.1.0-alpha"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                VERSION.parse_version(value)

    def test_mismatched_tag_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "does not match vcpkg.json"):
            VERSION.validate_tag("v0.1.0-alpha.1", "0.1.0-alpha.2")

    def test_tag_build_writes_validated_github_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "output"
            env = os.environ | {
                "GITHUB_REF_TYPE": "tag",
                "GITHUB_REF_NAME": f"v{self.metadata['version']}",
                "GITHUB_OUTPUT": str(output),
            }
            subprocess.run(
                [sys.executable, str(SCRIPT), "--github-output"],
                cwd=REPOSITORY,
                env=env,
                check=True,
            )
            self.assertEqual(
                output.read_text(encoding="utf-8").splitlines(),
                [
                    f"version={self.metadata['version']}",
                    f"short_version={self.metadata['short_version']}",
                    f"mac_bundle_version={self.metadata['mac_bundle_version']}",
                ],
            )

    def test_tag_build_rejects_mismatch(self):
        mismatched_version = (
            "0.0.1" if self.metadata["version"] == "0.0.0" else "0.0.0"
        )
        env = os.environ | {
            "GITHUB_REF_TYPE": "tag",
            "GITHUB_REF_NAME": f"v{mismatched_version}",
        }
        result = subprocess.run(
            [sys.executable, str(SCRIPT)],
            cwd=REPOSITORY,
            env=env,
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("does not match vcpkg.json", result.stderr)


if __name__ == "__main__":
    unittest.main()
