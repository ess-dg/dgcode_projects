"""Helpers for the LarmorTests test scripts (test-only code): samples of tubes and pixels of the rear bank in the
current and the old tube numbering (see Larmor2022Bank)."""

NSTRAWS = 7

def sample_tubes(nPacks, oldNumbering):
    """Current numbering (layer by layer): tube ids covering all 4 layers (front 0,1 and back 2,3) and the first,
    second, middle and last rows of the bank. Old numbering (pack by pack): all 8 tubes of the first and last
    packs."""
    if oldNumbering:
        return list(range(8)) + [8 * (nPacks - 1) + i for i in range(8)]
    rows = 2 * nPacks
    return [layer * rows + row for layer in range(4) for row in sorted({0, 1, rows // 2, rows - 1})]

def sample_pixels(nPacks, nPixelsPerStraw, oldNumbering):
    """Sample of (tubeId, strawId, inStrawPixel, pixelId): all 7 straws of the first sampled tube, and for the other
    sampled tubes one straw each (cycling through the straws), with pixel positions alternating between the first,
    last and a middle pixel of the straw; and the very first and very last pixels of the bank. (The sample of
    G4GeoLokiTests/testhelpers.sample_pixels for the rear bank.)"""
    N = nPixelsPerStraw
    res = []
    inStrawChoices = [0, N - 1, N // 2 - 1, N // 2]
    k = 0
    for it, t in enumerate(sample_tubes(nPacks, oldNumbering)):
        straws = range(NSTRAWS) if it == 0 else [(t + it) % NSTRAWS]
        for s in straws:
            j = inStrawChoices[k % len(inStrawChoices)]
            k += 1
            res.append((t, s, j, (t * NSTRAWS + s) * N + j))
    nTubes = 8 * nPacks
    for e in ((0, 0, 0, 0), (nTubes - 1, NSTRAWS - 1, N - 1, nTubes * NSTRAWS * N - 1)):
        if e not in res:
            res.append(e)
    return res
