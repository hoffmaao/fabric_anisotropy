"""Shared QC predicates for the CSARP_fabric products fabric_task writes.

Kept in one place so the figures and opr_fabric/test/test_fabric_task.m
cannot drift apart on what counts as a trustworthy interval.
"""
import numpy as np


def _as_2d(a):
    a = np.asarray(a, float)
    if a.ndim == 0:
        return a.reshape(1, 1)
    if a.ndim == 1:
        return a.reshape(-1, 1)   # single-block frame
    return a


def interpolated_intervals(d):
    """(nint x nblk) bool of intervals a fabricated node contaminates.

    dlam_interpolated flags the NODE at an interval's bottom edge whose dtau
    was interpolated across a masked gap (a waveform-combine seam, or a run
    of incoherent bins) instead of measured there. Layer stripping solves
    top-down from cumulative dtau, so interval k's dlam is dtau(k) -
    dtau(k-1): a fabricated node k corrupts interval k AND interval k+1.
    This mirrors seam_gap in opr_fabric/test/test_fabric_task.m.

    Returns None for products written before the flag existed, which the
    callers treat as "nothing to drop".
    """
    itp = d.get('dlam_interpolated')
    if itp is None or np.size(itp) == 0:
        return None
    itp = _as_2d(itp)
    bad = itp == 1
    bad[1:, :] |= itp[:-1, :] == 1
    return bad
