"""Check that a cross-compiled wheel is for the platform it claims.

A cross-build that silently produced a host binary passes every other test in
the pipeline: the wheel exists, its name carries the right tag, and on an
x86-64 runner it would even import and run. Nothing would catch it until a user
on the target platform tried to install it. So the machine type is read out of
the compiled extension itself, from the object file's own header, and compared
against what was asked for.

For Linux the glibc versions the extension needs are checked too. The wheel's
manylinux tag is a promise about which systems it will install on, and pip
believes it; the promise is only kept if the symbols the binary references stay
below the version named. `maturin --zig` sets that floor deliberately, which is
the whole reason for building this way, and this is what says it worked.

    python3 ci/check_wheel.py <wheel> <x86_64|aarch64>
"""

from __future__ import annotations

import glob
import re
import struct
import sys
import zipfile

# ELF `e_machine` and PE `Machine`, for the architectures this project builds
# for. Two formats because a Windows extension is a PE .pyd and a Linux one an
# ELF .so, and the point of the check is not to trust the file name.
ELF_MACHINES = {0x3E: "x86_64", 0xB7: "aarch64"}
PE_MACHINES = {0x8664: "x86_64", 0xAA64: "aarch64"}

# The highest glibc a wheel may require, by the tag it carries. `maturin --zig`
# targets 2.17; the manylinux container route reached 2.28. Both are recorded
# so that a change of build strategy shows up here as a failure rather than as
# a wheel that quietly stops installing on older systems.
MANYLINUX_CEILINGS = {
    "manylinux_2_17": (2, 17),
    "manylinux2014": (2, 17),
    "manylinux_2_28": (2, 28),
    "manylinux_2_34": (2, 34),
}


def extension_of(wheel_path):
    """The compiled extension inside the wheel, as raw bytes."""

    with zipfile.ZipFile(wheel_path) as archive:
        names = [n for n in archive.namelist() if n.endswith((".so", ".pyd", ".dylib"))]

        if len(names) != 1:
            raise SystemExit(
                f"expected exactly one compiled extension in {wheel_path}, found {names}"
            )

        return names[0], archive.read(names[0])


def machine_of(blob):
    """The architecture the object file was built for, read from its header."""

    if blob[:4] == b"\x7fELF":
        # e_machine is a half-word at offset 18, in the endianness byte 5 names.
        endian = "<" if blob[5] == 1 else ">"
        (value,) = struct.unpack_from(endian + "H", blob, 18)
        return ELF_MACHINES.get(value, f"unknown ELF machine 0x{value:x}")

    if blob[:2] == b"MZ":
        # The PE header's offset lives at 0x3C, and Machine follows its signature.
        (pe_offset,) = struct.unpack_from("<I", blob, 0x3C)
        if blob[pe_offset:pe_offset + 4] != b"PE\0\0":
            raise SystemExit("MZ header without a PE signature")
        (value,) = struct.unpack_from("<H", blob, pe_offset + 4)
        return PE_MACHINES.get(value, f"unknown PE machine 0x{value:x}")

    raise SystemExit("the extension is neither ELF nor PE")


def glibc_versions(blob):
    """Every `GLIBC_x.y` the extension names.

    Scanned out of the raw bytes rather than parsed out of `.gnu.version_r`:
    the version strings are stored literally, so a scan finds all of them, and
    a full ELF section parser would be a lot of code to reach the same list.
    """

    found = set()
    for match in re.finditer(rb"GLIBC_(\d+)\.(\d+)", blob):
        found.add((int(match.group(1)), int(match.group(2))))

    return found


def main(argv):

    if len(argv) != 3:
        raise SystemExit(__doc__)

    # The caller passes a glob, since the version is in the file name.
    matches = sorted(glob.glob(argv[1]))
    if len(matches) != 1:
        raise SystemExit(f"expected one wheel matching {argv[1]!r}, found {matches}")

    wheel_path, expected = matches[0], argv[2]

    name, blob = extension_of(wheel_path)
    machine = machine_of(blob)

    print(f"{wheel_path}")
    print(f"  extension: {name}")
    print(f"  machine:   {machine}")

    if machine != expected:
        raise SystemExit(f"FAIL: built for {machine}, expected {expected}")

    versions = glibc_versions(blob)
    if versions:
        highest = max(versions)
        print(f"  glibc:     up to {highest[0]}.{highest[1]}")

        ceiling = next(
            (v for tag, v in MANYLINUX_CEILINGS.items() if tag in wheel_path.replace("-", "_")),
            None,
        )
        if ceiling is None:
            raise SystemExit(f"FAIL: no manylinux tag recognised in {wheel_path}")
        if highest > ceiling:
            raise SystemExit(
                f"FAIL: needs glibc {highest[0]}.{highest[1]} but the tag promises "
                f"no more than {ceiling[0]}.{ceiling[1]}"
            )
        print(f"  tag allows up to {ceiling[0]}.{ceiling[1]}: ok")

    print("  ok")


if __name__ == "__main__":
    main(sys.argv)
