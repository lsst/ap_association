# This file is part of ap_association.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Helpers that apply shutter-corrected, per-source mid-exposure times
(`lsst.ip.isr.shutterTiming`) to DiaSource and DiaForcedSource catalogs.

A timing failure never fails a quantum: every entry point here catches
unexpected exceptions, logs a warning and falls back to the header midpoint
(``visitInfo.date``).
"""

__all__ = ["ShutterTimingCounts", "applyShutterTiming", "computeShutterTimingSafely", "logShutterTiming",
           "recordShutterTimingMetadata"]

import dataclasses
import math

import numpy as np

from lsst.ip.isr.shutterTiming import ShutterTimingStatus, computeShutterTiming


@dataclasses.dataclass(frozen=True)
class ShutterTimingCounts:
    """How many sources got which time.  The three counts partition the
    sources.
    """

    nCorrected: int = 0
    """Sources given a corrected time with per-source status OK."""
    nDegraded: int = 0
    """Sources given a corrected time with per-source status DEGRADED."""
    nFallback: int = 0
    """Sources that kept the header midpoint (UNAVAILABLE, non-finite time,
    or a timing failure).
    """


@dataclasses.dataclass(frozen=True)
class _FailedTiming:
    """Stand-in for a `~lsst.ip.isr.shutterTiming.ShutterTiming` whose
    computation raised; only used for metadata and logging.
    """

    message: str
    status: ShutterTimingStatus = ShutterTimingStatus.UNAVAILABLE
    flags: int = 0

    def summary(self):
        return {"status": self.status.name, "flags": 0, "message": self.message,
                "centerMinusHeaderMid": math.nan}


def computeShutterTimingSafely(metadata, detector, config, log):
    """Compute the shutter timing of one detector, never raising.

    Parameters
    ----------
    metadata : `lsst.daf.base.PropertyList`
        Exposure metadata carrying the shutter cards.
    detector : `lsst.afw.cameraGeom.Detector`
        The exposure's detector.
    config : `lsst.ip.isr.shutterTiming.ShutterTimingConfig`
        Configuration of the computation.
    log : `logging.Logger`
        Logger for the warning on failure.

    Returns
    -------
    timing : `lsst.ip.isr.shutterTiming.ShutterTiming` or `None`
        The timing, or `None` if the computation raised (logged as a
        warning).
    error : `str`
        "" on success, else a description of the exception.
    """
    try:
        return computeShutterTiming(metadata, detector, config), ""
    except Exception as e:
        error = f"{type(e).__name__}: {e}"
        log.warning("Shutter timing computation failed (%s); keeping the header midpoint.",
                    error, exc_info=e)
        return None, error


def applyShutterTiming(timing, x, y, headerMid, log):
    """Per-source mid-exposure times, falling back to the header midpoint.

    Parameters
    ----------
    timing : `lsst.ip.isr.shutterTiming.ShutterTiming` or `None`
        The detector's timing; `None` gives the header midpoint everywhere.
    x, y : array-like
        Pixel positions of the sources.
    headerMid : `float`
        Header midpoint (``visitInfo.date``, MJD TAI), used where the
        per-source status is UNAVAILABLE or the corrected time is not finite.
    log : `logging.Logger`
        Logger for the warning if the evaluation raises.

    Returns
    -------
    times : `numpy.ndarray`
        Per-source mid-exposure times, MJD TAI.  DEGRADED sources get the
        corrected time.
    counts : `ShutterTimingCounts`
        How many sources got a corrected (OK or DEGRADED) or a fallback time.
    """
    n = len(x)
    times = np.full(n, headerMid, dtype=np.float64)
    if timing is None or n == 0:
        return times, ShutterTimingCounts(nFallback=n)
    try:
        if timing.status == ShutterTimingStatus.UNAVAILABLE:
            return times, ShutterTimingCounts(nFallback=n)
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        corrected = np.broadcast_to(np.asarray(timing.tMidMjdTai(x, y), dtype=np.float64), (n,))
        status = np.broadcast_to(np.asarray(timing.sourceStatus(x, y)), (n,))
        use = (status != ShutterTimingStatus.UNAVAILABLE) & np.isfinite(corrected)
        degraded = use & (status == ShutterTimingStatus.DEGRADED)
    except Exception as e:
        log.warning("Shutter timing evaluation failed (%s: %s); keeping the header midpoint.",
                    type(e).__name__, e, exc_info=e)
        return times, ShutterTimingCounts(nFallback=n)
    times[use] = corrected[use]
    nUse = int(np.count_nonzero(use))
    nDegraded = int(np.count_nonzero(degraded))
    return times, ShutterTimingCounts(nCorrected=nUse - nDegraded, nDegraded=nDegraded,
                                      nFallback=n - nUse)


def _summary(timing, error):
    """Detector-level scalars for metadata and logs, never raising."""
    if timing is None:
        timing = _FailedTiming(error or "no shutter timing")
    try:
        summary = dict(timing.summary())
    except Exception as e:
        summary = {"message": f"summary() failed: {type(e).__name__}: {e}"}
    try:
        status = ShutterTimingStatus(timing.status).name
    except Exception:
        status = str(summary.get("status", "UNAVAILABLE"))
    try:
        flags = int(timing.flags)
    except Exception:
        flags = int(summary.get("flags", 0))
    message = summary.get("message", getattr(timing, "message", ""))
    try:
        offset = float(summary.get("centerMinusHeaderMid", math.nan))
    except (TypeError, ValueError):
        offset = math.nan
    return status, flags, str(message), offset


def recordShutterTimingMetadata(metadata, timing, counts=None, error=""):
    """Write the shutter-timing summary into task metadata.

    Parameters
    ----------
    metadata : `lsst.pipe.base.TaskMetadata`
        The task's metadata (``self.metadata``).
    timing : `lsst.ip.isr.shutterTiming.ShutterTiming` or `None`
        The detector's timing; `None` if its computation failed.
    counts : `ShutterTimingCounts`, optional
        Per-source counts; if `None`, the ``nShutter*`` keys are not written.
    error : `str`, optional
        Description of the failure when ``timing`` is `None`.

    Notes
    -----
    Keys: ``shutterTimingStatus`` (name), ``shutterTimingFlags`` (int),
    ``shutterTimingMessage``, ``shutterTimingCenterMinusHeaderMid`` (s), and
    ``nShutterCorrected``, ``nShutterDegraded``, ``nShutterFallback``.
    """
    status, flags, message, offset = _summary(timing, error)
    metadata["shutterTimingStatus"] = status
    metadata["shutterTimingFlags"] = flags
    metadata["shutterTimingMessage"] = message
    metadata["shutterTimingCenterMinusHeaderMid"] = offset
    if counts is not None:
        metadata["nShutterCorrected"] = counts.nCorrected
        metadata["nShutterDegraded"] = counts.nDegraded
        metadata["nShutterFallback"] = counts.nFallback


def logShutterTiming(log, timing, counts=None, error="", what="sources"):
    """Log one INFO line with the shutter-timing summary, plus a WARNING if
    the detector-level status is UNAVAILABLE.

    Parameters
    ----------
    log : `logging.Logger`
        The task's logger.
    timing : `lsst.ip.isr.shutterTiming.ShutterTiming` or `None`
        The detector's timing; `None` if its computation failed.
    counts : `ShutterTimingCounts`, optional
        Per-source counts to include.
    error : `str`, optional
        Description of the failure when ``timing`` is `None`.
    what : `str`, optional
        Name of the sources counted, for the log line.
    """
    status, flags, message, offset = _summary(timing, error)
    countText = ""
    if counts is not None:
        countText = (f"; {what}: {counts.nCorrected} corrected, {counts.nDegraded} degraded, "
                     f"{counts.nFallback} header midpoint")
    log.info("Shutter timing: status %s, flags %#x, centre - header midpoint %.4f s%s.",
             status, flags, offset, countText)
    if status == ShutterTimingStatus.UNAVAILABLE.name:
        log.warning("Shutter timing UNAVAILABLE (%s): using the header midpoint for all %s.",
                    message, what)
