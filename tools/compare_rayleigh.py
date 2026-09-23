#!/usr/bin/env python3
"""Compare CPU and GPU optools output tables for a 1D atmosphere."""

import argparse
import array
import math
import sys
from pathlib import Path


def profile_layer_count(path):
    """Read nlay using the same three-record header layout as read_prf."""
    try:
        with path.open("r", encoding="utf-8") as stream:
            next(stream)
            next(stream)
            line = next(stream)
        return int(line.split()[0])
    except (OSError, StopIteration, ValueError, IndexError) as exc:
        raise ValueError("cannot read the layer count from {}: {}".format(path, exc))


def wavelength_grid(path):
    """Read the wavelength count and wavelength values from wavelengths.wl."""
    try:
        with path.open("r", encoding="utf-8") as stream:
            lines = [
                line.strip()
                for line in stream
                if line.strip() and not line.lstrip().startswith(("#", "!"))
            ]

        nwl = int(lines[0].split()[0])
        values = [float(lines[index + 1].split()[1]) for index in range(nwl)]
        return nwl, values
    except (OSError, ValueError, IndexError) as exc:
        raise ValueError("cannot read the wavelength grid from {}: {}".format(path, exc))


def read_table(path, nwl, nlay):
    """Read native-endian float32 direct-access records."""
    expected_values = nwl * nlay
    expected_bytes = expected_values * 4

    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise ValueError("cannot read {}: {}".format(path, exc))

    if len(raw) != expected_bytes:
        raise ValueError(
            "{} contains {} bytes; expected {} ({} wavelengths x {} layers x 4 bytes)".format(
                path, len(raw), expected_bytes, nwl, nlay
            )
        )

    values = array.array("f")
    values.frombytes(raw)
    if values.itemsize != 4:
        raise ValueError("this Python build does not use four-byte native floats")

    return values, raw


def parse_arguments(quantity="Rayleigh"):
    parser = argparse.ArgumentParser(
        description=(
            "Compare CPU and GPU {}.cmcrt files. Files are interpreted as "
            "native-endian float32 direct-access records.".format(quantity)
        )
    )
    parser.add_argument(
        "cpu", type=Path, help="{}.cmcrt produced by goptools_cpu".format(quantity)
    )
    parser.add_argument(
        "gpu", type=Path, help="{}.cmcrt produced by goptools_gpu".format(quantity)
    )
    parser.add_argument("--profile", required=True, type=Path, help="1D atmospheric .prf file")
    parser.add_argument(
        "--wavelengths", required=True, type=Path, help="wavelengths.wl used for both runs"
    )
    parser.add_argument("--rtol", type=float, default=5.0e-6, help="relative tolerance")
    parser.add_argument("--atol", type=float, default=0.0, help="absolute tolerance")
    return parser.parse_args()


def location(index, nlay, wavelengths):
    wavelength_index, layer_index = divmod(index, nlay)
    return wavelength_index, layer_index, wavelengths[wavelength_index]


def main(quantity="Rayleigh"):
    args = parse_arguments(quantity)

    if args.rtol < 0.0 or args.atol < 0.0:
        print("ERROR: tolerances must be non-negative", file=sys.stderr)
        return 2

    try:
        nlay = profile_layer_count(args.profile)
        nwl, wavelengths = wavelength_grid(args.wavelengths)
        cpu, cpu_raw = read_table(args.cpu, nwl, nlay)
        gpu, gpu_raw = read_table(args.gpu, nwl, nlay)
    except ValueError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2

    cpu_nonfinite = sum(not math.isfinite(value) for value in cpu)
    gpu_nonfinite = sum(not math.isfinite(value) for value in gpu)
    if cpu_nonfinite or gpu_nonfinite:
        print(
            "FAIL: found {} CPU and {} GPU non-finite values".format(
                cpu_nonfinite, gpu_nonfinite
            )
        )
        return 1

    maximum_absolute = -1.0
    maximum_relative = -1.0
    maximum_absolute_index = 0
    maximum_relative_index = 0
    failure_count = 0
    first_failure = None

    for index, (cpu_value, gpu_value) in enumerate(zip(cpu, gpu)):
        absolute_error = abs(float(cpu_value) - float(gpu_value))
        magnitude = max(abs(float(cpu_value)), abs(float(gpu_value)))
        relative_error = absolute_error / max(magnitude, sys.float_info.min)

        if absolute_error > maximum_absolute:
            maximum_absolute = absolute_error
            maximum_absolute_index = index
        if relative_error > maximum_relative:
            maximum_relative = relative_error
            maximum_relative_index = index

        if absolute_error > args.atol + args.rtol * magnitude:
            failure_count += 1
            if first_failure is None:
                first_failure = (index, cpu_value, gpu_value, absolute_error)

    bitwise_equal = sum(
        cpu_raw[offset : offset + 4] == gpu_raw[offset : offset + 4]
        for offset in range(0, len(cpu_raw), 4)
    )

    abs_wl, abs_layer, abs_wavelength = location(
        maximum_absolute_index, nlay, wavelengths
    )
    rel_wl, rel_layer, rel_wavelength = location(
        maximum_relative_index, nlay, wavelengths
    )

    print("{} CPU/GPU comparison".format(quantity))
    print("  shape:                 {} wavelengths x {} layers".format(nwl, nlay))
    print("  tolerances:            rtol={:.3e}, atol={:.3e}".format(args.rtol, args.atol))
    print(
        "  bitwise-equal values:  {} / {} ({:.2f}%)".format(
            bitwise_equal, len(cpu), 100.0 * bitwise_equal / len(cpu)
        )
    )
    print("  maximum absolute error: {:.8e}".format(maximum_absolute))
    print(
        "    at wavelength {} ({:.8g} um), layer {}".format(
            abs_wl + 1, abs_wavelength, abs_layer + 1
        )
    )
    print("  maximum relative error: {:.8e}".format(maximum_relative))
    print(
        "    at wavelength {} ({:.8g} um), layer {}".format(
            rel_wl + 1, rel_wavelength, rel_layer + 1
        )
    )
    print("  values outside tolerance: {}".format(failure_count))

    if first_failure is not None:
        index, cpu_value, gpu_value, absolute_error = first_failure
        wl_index, layer_index, wavelength = location(index, nlay, wavelengths)
        print(
            "FAIL: first mismatch at wavelength {} ({:.8g} um), layer {}: "
            "CPU={:.8e}, GPU={:.8e}, abs_error={:.8e}".format(
                wl_index + 1,
                wavelength,
                layer_index + 1,
                cpu_value,
                gpu_value,
                absolute_error,
            )
        )
        return 1

    print("PASS: CPU and GPU {} outputs agree within tolerance".format(quantity))
    return 0


if __name__ == "__main__":
    sys.exit(main())
