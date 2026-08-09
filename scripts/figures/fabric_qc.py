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


def ref_degenerate_intervals(d):
    """(nint x nblk) bool of intervals degenerate with the reference offset.

    dlam_ref_degenerate flags the interval whose dlam shares its
    information with the reference-offset nuisance of
    ptt.invertHorizontalFabricJoint (the shallowest interval when
    obs.zref is in play): quotable fabric starts below it. All false on
    the layer-stripping chain, so the predicate is safe on both chains.

    Returns None for products written before the flag existed, which the
    callers treat as "nothing to drop".
    """
    ref = d.get('dlam_ref_degenerate')
    if ref is None or np.size(ref) == 0:
        return None
    return _as_2d(ref) == 1


def dropped_intervals(d):
    """Union of every interval-drop predicate above; the figure scripts
    consume this one, so a flag added here reaches all of them at once.

    Returns None when the product predates all of the flags.
    """
    masks = [m for m in (interpolated_intervals(d),
                         ref_degenerate_intervals(d)) if m is not None]
    if not masks:
        return None
    out = masks[0]
    for m in masks[1:]:
        out = out | m
    return out
