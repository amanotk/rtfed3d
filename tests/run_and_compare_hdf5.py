#!/usr/bin/env python3

import argparse
import pathlib
import subprocess
import sys
import tempfile

import h5py
import numpy as np


def compare_values(path, expected, actual, rtol, atol):
    if isinstance(expected, np.ndarray) or isinstance(actual, np.ndarray):
        expected_arr = np.asarray(expected)
        actual_arr = np.asarray(actual)
        if expected_arr.shape != actual_arr.shape:
            raise AssertionError(
                f"Shape mismatch at {path}: expected {expected_arr.shape}, got {actual_arr.shape}"
            )
        if not np.allclose(
            expected_arr, actual_arr, rtol=rtol, atol=atol, equal_nan=True
        ):
            diff = np.max(np.abs(expected_arr - actual_arr))
            raise AssertionError(f"Value mismatch at {path}: max abs diff {diff}")
        return

    if isinstance(expected, (bytes, str)) or isinstance(actual, (bytes, str)):
        if expected != actual:
            raise AssertionError(
                f"Value mismatch at {path}: expected {expected!r}, got {actual!r}"
            )
        return

    if np.issubdtype(type(expected), np.integer) and np.issubdtype(
        type(actual), np.integer
    ):
        if expected != actual:
            raise AssertionError(
                f"Value mismatch at {path}: expected {expected}, got {actual}"
            )
        return

    if not np.isclose(expected, actual, rtol=rtol, atol=atol, equal_nan=True):
        raise AssertionError(
            f"Value mismatch at {path}: expected {expected}, got {actual}"
        )


def compare_attrs(expected_obj, actual_obj, path, rtol, atol):
    expected_keys = sorted(expected_obj.attrs.keys())
    actual_keys = sorted(actual_obj.attrs.keys())
    if expected_keys != actual_keys:
        raise AssertionError(
            f"Attribute keys mismatch at {path}: expected {expected_keys}, got {actual_keys}"
        )
    for key in expected_keys:
        compare_values(
            f"{path}@{key}", expected_obj.attrs[key], actual_obj.attrs[key], rtol, atol
        )


def compare_group(expected_group, actual_group, path, rtol, atol):
    compare_attrs(expected_group, actual_group, path, rtol, atol)

    expected_keys = sorted(expected_group.keys())
    actual_keys = sorted(actual_group.keys())
    if expected_keys != actual_keys:
        raise AssertionError(
            f"Object keys mismatch at {path}: expected {expected_keys}, got {actual_keys}"
        )

    for key in expected_keys:
        next_path = f"{path}/{key}" if path != "/" else f"/{key}"
        expected_obj = expected_group[key]
        actual_obj = actual_group[key]

        if isinstance(expected_obj, h5py.Dataset):
            if not isinstance(actual_obj, h5py.Dataset):
                raise AssertionError(f"Object type mismatch at {next_path}")
            compare_attrs(expected_obj, actual_obj, next_path, rtol, atol)
            compare_values(next_path, expected_obj[...], actual_obj[...], rtol, atol)
        else:
            if not isinstance(actual_obj, h5py.Group):
                raise AssertionError(f"Object type mismatch at {next_path}")
            compare_group(expected_obj, actual_obj, next_path, rtol, atol)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output-name", required=True)
    parser.add_argument("--rtol", type=float, default=1.0e-12)
    parser.add_argument("--atol", type=float, default=1.0e-12)
    args = parser.parse_args()

    executable = pathlib.Path(args.executable).resolve()
    config = pathlib.Path(args.config).resolve()
    reference = pathlib.Path(args.reference).resolve()

    with tempfile.TemporaryDirectory(prefix="rtfed3d-regression-") as tmpdir:
        tmpdir_path = pathlib.Path(tmpdir)
        subprocess.run(
            [str(executable), "-c", str(config)],
            cwd=tmpdir,
            check=True,
        )

        actual = tmpdir_path / args.output_name
        if not actual.exists():
            raise FileNotFoundError(f"Expected output file was not created: {actual}")

        with (
            h5py.File(reference, "r") as expected_h5,
            h5py.File(actual, "r") as actual_h5,
        ):
            compare_group(expected_h5, actual_h5, "/", args.rtol, args.atol)

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
