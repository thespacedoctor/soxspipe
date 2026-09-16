"""Print a data-only digest of every reduced product, so two runs can be compared.

FITS headers carry timestamps and paths that differ on every run, so only the data
is hashed: image arrays, and every column of every table, rendered at full precision.
"""

import hashlib
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits


def digest_array(values):
    """Return a digest of an array's exact bytes, in its own memory order."""
    array = np.ascontiguousarray(values)
    return hashlib.sha256(array.tobytes()).hexdigest()[:16]


def digest_file(path):
    """Return one digest line per HDU that carries data."""
    lines = []
    with fits.open(path) as hdus:
        for index, hdu in enumerate(hdus):
            if hdu.data is None:
                continue
            if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                for name in hdu.data.names:
                    lines.append(f"{path.name} hdu{index} {name} {digest_array(hdu.data[name])}")
            else:
                lines.append(f"{path.name} hdu{index} DATA {digest_array(hdu.data)}")
    return lines


def main(workspace):
    """Print every product digest, sorted, so the output is comparable line by line."""
    products = sorted(Path(workspace).rglob("*.fits"))
    lines = []
    for product in products:
        try:
            lines.extend(digest_file(product))
        except Exception as error:  # A PROBE MUST NEVER FAIL THE RUN
            lines.append(f"{product.name} UNREADABLE {error}")
    for line in sorted(lines):
        print(f"DIGEST {line}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
