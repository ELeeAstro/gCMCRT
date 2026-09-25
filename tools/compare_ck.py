#!/usr/bin/env python3
"""Compare CPU and GPU correlated-k output for a 1D atmosphere."""

import argparse
import array
import math
import sys
from pathlib import Path

from compare_rayleigh import profile_layer_count, wavelength_grid


def read_gordinates(path):
    """Read the g ordinates and weights written by optools."""
    try:
        with path.open("r", encoding="utf-8") as stream:
            lines = [
                line.strip()
                for line in stream
                if line.strip() and not line.lstrip().startswith(("#", "!"))
            ]
        n_g = int(lines[0].split()[0])
        if len(lines) != n_g + 1:
            raise ValueError(
                "contains {} g-ordinate rows; expected {}".format(len(lines) - 1, n_g)
            )
        values = []
        for index in range(n_g):
            fields = lines[index + 1].split()
            values.append((float(fields[0]), float(fields[1])))
        return n_g, values
    except (OSError, ValueError, IndexError) as exc:
        raise ValueError("cannot read {}: {}".format(path, exc))


def read_ck_table(path, n_wavelength, n_layer, n_g):
    """Read native-endian float32 CK direct-access records."""
    expected_values = n_wavelength * n_layer * n_g
    expected_bytes = expected_values * 4
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise ValueError("cannot read {}: {}".format(path, exc))

    if len(raw) != expected_bytes:
        raise ValueError(
            "{} contains {} bytes; expected {} "
            "({} wavelengths x {} layers x {} g ordinates x 4 bytes)".format(
                path,
                len(raw),
                expected_bytes,
                n_wavelength,
                n_layer,
                n_g,
            )
        )

    values = array.array("f")
    values.frombytes(raw)
    if values.itemsize != 4:
        raise ValueError("this Python build does not use four-byte native floats")
    return values, raw


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Compare CPU and GPU CK.cmcrt files as native-endian float32 "
            "direct-access records."
        )
    )
    parser.add_argument("cpu", type=Path, help="CK.cmcrt produced by goptools_cpu")
    parser.add_argument("gpu", type=Path, help="CK.cmcrt produced by goptools_gpu")
    parser.add_argument("--cpu-gord", required=True, type=Path)
    parser.add_argument("--gpu-gord", required=True, type=Path)
    parser.add_argument("--profile", required=True, type=Path)
    parser.add_argument("--wavelengths", required=True, type=Path)
    parser.add_argument("--rtol", type=float, default=5.0e-6)
    parser.add_argument("--atol", type=float, default=0.0)
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.rtol < 0.0 or args.atol < 0.0:
        print("ERROR: tolerances must be non-negative", file=sys.stderr)
        return 2

    try:
        n_layer = profile_layer_count(args.profile)
        n_wavelength, wavelengths = wavelength_grid(args.wavelengths)
        cpu_n_g, cpu_gord = read_gordinates(args.cpu_gord)
        gpu_n_g, gpu_gord = read_gordinates(args.gpu_gord)
        if cpu_n_g != gpu_n_g:
            raise ValueError(
                "g-ordinate counts differ: CPU {} versus GPU {}".format(
                    cpu_n_g, gpu_n_g
                )
            )
        cpu, cpu_raw = read_ck_table(
            args.cpu, n_wavelength, n_layer, cpu_n_g
        )
        gpu, gpu_raw = read_ck_table(
            args.gpu, n_wavelength, n_layer, cpu_n_g
        )
    except ValueError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2

    gord_failures = []
    for index, (cpu_pair, gpu_pair) in enumerate(zip(cpu_gord, gpu_gord), start=1):
        for label, cpu_value, gpu_value in (
            ("ordinate", cpu_pair[0], gpu_pair[0]),
            ("weight", cpu_pair[1], gpu_pair[1]),
        ):
            if not math.isfinite(cpu_value) or not math.isfinite(gpu_value):
                gord_failures.append((index, label, cpu_value, gpu_value))
            elif not math.isclose(
                cpu_value, gpu_value, rel_tol=args.rtol, abs_tol=args.atol
            ):
                gord_failures.append((index, label, cpu_value, gpu_value))

    cpu_nonfinite = sum(not math.isfinite(value) for value in cpu)
    gpu_nonfinite = sum(not math.isfinite(value) for value in gpu)
    if cpu_nonfinite or gpu_nonfinite:
        print(
            "FAIL: found {} CPU and {} GPU non-finite opacity values".format(
                cpu_nonfinite, gpu_nonfinite
            )
        )
        return 1

    maximum_absolute = -1.0
    maximum_relative = -1.0
    maximum_index = 0
    failures = 0
    for index, (cpu_value, gpu_value) in enumerate(zip(cpu, gpu)):
        absolute = abs(cpu_value - gpu_value)
        scale = max(abs(cpu_value), abs(gpu_value))
        relative = absolute / scale if scale > 0.0 else 0.0
        if absolute > maximum_absolute:
            maximum_absolute = absolute
            maximum_index = index
        maximum_relative = max(maximum_relative, relative)
        if absolute > args.atol + args.rtol * abs(cpu_value):
            failures += 1

    record_size = n_layer * cpu_n_g
    wavelength_index, within_record = divmod(maximum_index, record_size)
    layer_index, g_index = divmod(within_record, cpu_n_g)

    print(
        "Compared {} values ({} wavelengths x {} layers x {} g ordinates)".format(
            len(cpu), n_wavelength, n_layer, cpu_n_g
        )
    )
    print("Maximum absolute error: {:.9e}".format(maximum_absolute))
    print("Maximum relative error: {:.9e}".format(maximum_relative))
    print(
        "Largest absolute error at wavelength {} ({:.9g} um), layer {}, g {}".format(
            wavelength_index + 1,
            wavelengths[wavelength_index],
            layer_index + 1,
            g_index + 1,
        )
    )

    if gord_failures:
        index, label, cpu_value, gpu_value = gord_failures[0]
        print(
            "FAIL: {} g-ordinate values differ; first at g {} {}: CPU {} GPU {}".format(
                len(gord_failures), index, label, cpu_value, gpu_value
            )
        )
        return 1
    print("gord.cmcrt values match within tolerance")

    if failures:
        print(
            "FAIL: {} opacity values exceed rtol={} and atol={}".format(
                failures, args.rtol, args.atol
            )
        )
        return 1

    if cpu_raw == gpu_raw:
        print("PASS: CK.cmcrt files are bitwise identical")
    else:
        print("PASS: all opacity values are within tolerance")
    return 0


if __name__ == "__main__":
    sys.exit(main())
