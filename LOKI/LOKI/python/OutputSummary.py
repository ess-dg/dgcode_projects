"""Summary of the output files of a simulation + analysis in the current directory, for the tests (their logs).

  python3 -m LOKI.OutputSummary [--nonzero-only]

For each SimpleHists file (*.shist): the counters and, per histogram, its integral and mean (4 significant digits);
for each detection file (detectionEvents*.mcpl): the number of events and a checksum of their pixel ids; for each
mask file (*maskFile.xml): the number of masked ids and a checksum of them. --nonzero-only prints only which
histograms are filled and the number of detection events > 0, for runs with physics (whose numbers depend on the
platform), while runs with geantinos or without interactions (-lPL_Empty) are reproducible.
"""

import glob
import hashlib
import re
import sys


def _significant(value, digits=4):
    return float(f'{value:.{digits}g}') + 0.0


def shist_summary(path, nonzero_only):
    import SimpleHists as sh
    collection = sh.HistCollection(path)
    lines = [f'{path}:']
    for key in sorted(collection.keys):
        hist = collection.hist(key)
        if type(hist).__name__ == 'HistCounts':
            if nonzero_only:
                lines.append(f'  {key}: counters {", ".join(c.label for c in hist.counters)}')
            else:
                lines.append(f'  {key}: ' + ', '.join(f'{c.label} {_significant(c.value, 6)}' for c in hist.counters))
            continue
        integral = hist.integral
        if nonzero_only:
            lines.append(f'  {key}: {"filled" if integral > 0 else "empty"}')
        elif integral > 0:
            if hasattr(hist, 'getMeanY'):
                means = [_significant(hist.getMeanX()), _significant(hist.getMeanY())]
            else:
                means = [_significant(hist.getMean())]
            lines.append(f'  {key}: integral {_significant(integral, 6)}, mean {", ".join(map(str, means))}')
        else:
            lines.append(f'  {key}: empty')
    return lines


def mcpl_summary(path, nonzero_only):
    import mcpl
    file = mcpl.MCPLFile(path)
    if nonzero_only:
        return [f'{path}: detection events > 0: {file.nparticles > 0}']
    ids = []
    for block in file.particle_blocks:
        values = block.userflags if file.opt_userflags else block.ekin
        ids.extend(int(round(v)) for v in values)
    digest = hashlib.sha256(','.join(map(str, ids)).encode()).hexdigest()[:16]
    return [f'{path}: {file.nparticles} detection events, sha256 of the pixel ids {digest}']


def mask_summary(path):
    text = open(path).read()
    ids = [int(v) for v in re.search(r'<detids>(.*)</detids>', text, re.S).group(1).split(',') if v.strip()]
    digest = hashlib.sha256(','.join(map(str, ids)).encode()).hexdigest()[:16]
    return [f'{path}: {len(ids)} masked ids, sha256 {digest}']


def main(argv):
    nonzero_only = '--nonzero-only' in argv
    lines = []
    for path in sorted(glob.glob('*.shist')):
        lines += shist_summary(path, nonzero_only)
    for path in sorted(glob.glob('detectionEvents*.mcpl')):
        lines += mcpl_summary(path, nonzero_only)
    for path in sorted(glob.glob('*maskFile.xml')):
        lines += mask_summary(path)
    if not lines:
        print('no output files found')
        return 1
    print('\n'.join(lines))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
