"""Shared helpers for the G4GeoLokiTests test scripts (test-only code).

Conventions used by the tests:
  * All lengths are in mm (Geant4 / dgcode internal units).
  * Pixel ids are zero based, as in the Geant4 geometry and the analysis
    programs (the ICD/IDF ids are pixelId+1).
  * tube numbering (GeoBCSBanks):
       tubeId = layer * (2*nPacks) + row,  row = 2*pack + (in-pack tube >= 4)
    (row/pack order reversed for the upside-down banks 1,2,5,6).
"""

import contextlib
import os
import sys
import tempfile

NBANKS = 9
NSTRAWS = 7

def fmt(v, prec=4):
    """Format float with fixed precision, avoiding '-0.0000'."""
    s = f'{v:.{prec}f}'
    if s.startswith('-') and float(s) == 0.0:
        s = s[1:]
    return s

def fmtvec(v, prec=4):
    return '(' + ', '.join(fmt(e, prec) for e in v) + ')'

@contextlib.contextmanager
def captured_fd_output():
    """Capture everything written to the C-level stdout/stderr file
    descriptors (e.g. G4cout from C++ code) into a list of lines, which is
    available as the context value after the block."""
    sys.stdout.flush()
    sys.stderr.flush()
    lines = []
    with tempfile.TemporaryFile(mode='w+b') as tmp:
        saved_out, saved_err = os.dup(1), os.dup(2)
        os.dup2(tmp.fileno(), 1)
        os.dup2(tmp.fileno(), 2)
        try:
            yield lines
        finally:
            try:
                import G4Utils
                G4Utils.flush()
            except Exception:
                pass
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved_out, 1)
            os.dup2(saved_err, 2)
            os.close(saved_out)
            os.close(saved_err)
            tmp.seek(0)
            lines.extend(tmp.read().decode('utf-8', errors='replace').splitlines())

def tube_id_new(nPacks, layer, row):
    return layer * 2 * nPacks + row

def sample_tubes_new(bankId, nPacks):
    """Tube ids covering all 4 layers (front 0,1 and back 2,3)
    and the first, second, middle and last rows of the bank."""
    rows = 2 * nPacks
    tubes = []
    for layer in range(4):
        for row in sorted({0, 1, rows // 2, rows - 1}):
            tubes.append(tube_id_new(nPacks, layer, row))
    return tubes

def sample_pixels(aimHelper, nPixelsPerStraw, banks=range(NBANKS)):
    """Well-chosen sample of (bankId, tubeId, strawId, inStrawPixel, pixelId).

    For each bank: all 7 straws of the first sampled tube, and for the other
    sampled tubes one straw each (cycling through the straws), with pixel
    positions alternating between the first, last and a middle pixel of the
    straw. The very first and very last pixels of each bank are always included.
    """
    import G4GeoLokiTests.TestUtils as TU
    N = nPixelsPerStraw
    res = []
    for b in banks:
        nPacks = TU.getNumberOfPacksByBankId(b)
        offset = aimHelper.getBankPixelOffset(b)
        tubes = sample_tubes_new(b, nPacks)
        inStrawChoices = [0, N - 1, N // 2 - 1, N // 2]
        k = 0
        for it, t in enumerate(tubes):
            straws = range(NSTRAWS) if it == 0 else [(t + it) % NSTRAWS]
            for s in straws:
                j = inStrawChoices[k % len(inStrawChoices)]
                k += 1
                res.append((b, t, s, j, offset + (t * NSTRAWS + s) * N + j))
        # first and last pixel of the bank:
        nTubes = TU.getNumberOfTubes(b)
        first = (b, 0, 0, 0, offset)
        last = (b, nTubes - 1, NSTRAWS - 1, N - 1, offset + nTubes * NSTRAWS * N - 1)
        for e in (first, last):
            if e not in res:
                res.append(e)
    return res

def straw_axis_dir(aimHelper, pixelId, N):
    """Unit vector along the straw in the direction of increasing pixel index
    (derived from AimHelper centres of neighbouring pixels), and the pixel pitch."""
    j = pixelId % N
    if N < 2:
        raise ValueError('need at least 2 pixels per straw')
    p0, p1 = (pixelId, pixelId + 1) if j < N - 1 else (pixelId - 1, pixelId)
    c0 = aimHelper.getPixelCentreCoordinates(p0)
    c1 = aimHelper.getPixelCentreCoordinates(p1)
    d = [c1[i] - c0[i] for i in range(3)]
    norm = sum(e * e for e in d) ** 0.5
    return [e / norm for e in d], norm

def perpendicular_dirs(u):
    """Two unit vectors perpendicular to unit vector u (and each other)."""
    # pick the world axis least aligned with u
    ax = min(range(3), key=lambda i: abs(u[i]))
    e = [0.0, 0.0, 0.0]
    e[ax] = 1.0
    # v = e - (e.u)u
    dot = sum(e[i] * u[i] for i in range(3))
    v = [e[i] - dot * u[i] for i in range(3)]
    nv = sum(c * c for c in v) ** 0.5
    v = [c / nv for c in v]
    w = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
    return v, w

def add(p, d, s):
    return [p[i] + s * d[i] for i in range(3)]

class Checker:
    """Collects check results; prints failures immediately and a summary at the end."""
    def __init__(self):
        self.nchecks = 0
        self.nfail = 0
    def check(self, ok, msg):
        self.nchecks += 1
        if not ok:
            self.nfail += 1
            print('FAILURE:', msg)
        return ok
    def finish(self, title=''):
        print(f'{title}checks performed: {self.nchecks}, failures: {self.nfail}')
        if self.nfail:
            print('TEST FAILED')
            sys.exit(1)
