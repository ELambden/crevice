"""Approximate confidence intervals for means of serially correlated frame series.

Consecutive MD frames are correlated, so treating them as independent
samples gives confidence intervals that are too narrow. This module:

1. estimates the *statistical inefficiency* ``g`` of each series (the factor
   by which correlation inflates the variance of the mean; ``N/g`` is the
   effective number of independent frames), and
2. averages the frames into contiguous batches at least a few correlation
   times long and bootstraps those batch means with a studentized
   (percentile-t) bootstrap. Whole batches are resampled jointly across all
   observables, which preserves correlation between neighbouring positions
   and allows a simultaneous band over all of them.

A studentized bootstrap resamples whole contiguous batch means, preserving
spatial covariance and accounting for uncertainty in the estimated standard
error. It does not quantify force-field, sampling-state or geometry errors.

Assumptions
    The trajectory is stationary and representative, and the chosen batch
    length captures the serial dependence. Neither is tested here; the
    result records ``stationarity_validated=False``.
"""
from __future__ import annotations
import math


def statistical_inefficiency(values):
    """Statistical inefficiency ``g`` of a one-dimensional series.

    ``g = 1 + 2 * sum(rho_k)`` with the normalised autocorrelation ``rho_k``
    (computed by FFT). The sum uses Geyer's initial positive sequence: pairs
    ``rho_k + rho_(k+1)`` for odd ``k`` are added until the first non-positive
    pair, up to lag ``N/2``. The result is clipped to ``[1, N]``.

    Parameters
    ----------
    values : sequence of float
        Finite series of at least two values.

    Returns
    -------
    float
        ``g`` (1 means uncorrelated). A constant series returns 1.0.

    Raises
    ------
    ValueError
        If the input is not a finite one-dimensional series of length >= 2.

    Examples
    --------
    >>> from crevice.uncertainty import statistical_inefficiency
    >>> statistical_inefficiency([1.0, -1.0] * 50)      # alternating: no positive pairs
    1.0
    >>> statistical_inefficiency([3.0] * 10)
    1.0
    """
    import numpy as np
    x=np.asarray(values,dtype=float)
    if x.ndim!=1 or len(x)<2 or not np.isfinite(x).all():
        raise ValueError('Autocorrelation requires a finite one-dimensional series')
    x=x-x.mean(); n=len(x)
    if np.dot(x,x)<=np.finfo(float).eps*n:
        return 1.0
    size=1 << (2*n-1).bit_length()
    f=np.fft.rfft(x,n=size)
    ac=np.fft.irfft(f*f.conj(),n=size)[:n]
    rho=ac/ac[0]
    total=0.0
    for k in range(1,n//2,2):
        pair=float(rho[k]+rho[k+1])
        if pair<=0:break
        total+=pair
    return min(float(n),max(1.0,1+2*total))


def block_mean_confidence(matrix, *, confidence=.95, block_length=None,
                          replicates=2000, seed=20260911, min_blocks=8):
    """Pointwise percentile-t and simultaneous max-t bands for the mean.

    Only columns observed in every input frame receive intervals. Full rows are
    averaged into contiguous, balanced batches and resampled jointly. Missing
    rows are never compressed into an artificial time series. Auto target batch
    length is ``max(ceil(N**(1/3)), ceil(5*max(g)))``; at least eight batches
    are needed.

    Algorithm:

    1. For each fully observed column ``j`` compute ``g_j``
       (:func:`statistical_inefficiency`).
    2. Split the ``N`` frames into ``floor(N / L)`` contiguous batches of nearly
       equal length and average each.
    3. Estimate each column's residual batch inefficiency; columns with fewer
       than ``min_blocks`` effective batches get no interval.
    4. Resample batches with replacement ``replicates`` times and form
       studentized statistics ``(boot_mean - mean) / boot_se``.
    5. Pointwise interval: percentile-t quantiles, widened to at least the
       Student-t critical value. Simultaneous band: the ``confidence`` quantile
       of the maximum absolute studentized statistic over columns, widened to
       at least a Bonferroni-t critical value. Standard errors are inflated by
       ``sqrt(residual batch inefficiency)``.

    Parameters
    ----------
    matrix : array_like, shape (n_frames, n_observables)
        Values per frame (rows in time order). ``NaN`` marks missing
        observations; infinities are not allowed.
    confidence : float, default 0.95
        Two-sided confidence level, strictly between 0 and 1.
    block_length : int, optional
        Batch length ``L`` in frames (1 to ``n_frames``). Chosen automatically
        when omitted.
    replicates : int, default 2000
        Bootstrap replicates (integer, at least 200).
    seed : int, default 20260911
        Seed of the NumPy random generator; results are reproducible.
    min_blocks : int, default 8
        Minimum number of batches (and of effective batches per column).

    Returns
    -------
    dict
        ``status`` is ``"estimated"`` or a reason for no interval
        (``"insufficient_independent_blocks"`` when fewer than 20 frames or
        fewer than ``min_blocks`` batches, ``"no_complete_coverage"``,
        ``"insufficient_effective_batches"``,
        ``"insufficient_batch_variation"``). Lists ``pointwise_lower``,
        ``pointwise_upper``, ``simultaneous_lower`` and ``simultaneous_upper``
        have one entry per column (``None`` where no interval was estimated).
        Also reported: the batch length and its source, per-column
        ``statistical_inefficiency`` and ``effective_frames``, batch counts,
        the critical value, and the difference between second-half and
        first-half means (a crude drift indicator).

    Raises
    ------
    ValueError
        For an invalid matrix, confidence, replicates, ``min_blocks`` or
        ``block_length``.
    """
    import numpy as np
    x=np.asarray(matrix,dtype=float)
    if x.ndim!=2 or x.shape[0]<1 or x.shape[1]<1 or np.isinf(x).any():
        raise ValueError('matrix must be a nonempty frames-by-observables array without infinities')
    if not math.isfinite(confidence) or not 0<confidence<1:
        raise ValueError('confidence must lie strictly between 0 and 1')
    if isinstance(replicates,bool) or int(replicates)!=replicates or replicates<200:
        raise ValueError('replicates must be an integer >=200')
    if isinstance(min_blocks,bool) or int(min_blocks)!=min_blocks or min_blocks<2:
        raise ValueError('min_blocks must be an integer >=2')
    n,m=x.shape
    if block_length is not None and (isinstance(block_length,bool) or int(block_length)!=block_length or not 1<=block_length<=n):
        raise ValueError('block_length must be an integer between 1 and the frame count')
    eligible=np.isfinite(x).all(axis=0)
    g=np.full(m,np.nan)
    if n>=2:
        for j in np.flatnonzero(eligible):g[j]=statistical_inefficiency(x[:,j])
    length=int(block_length) if block_length is not None else max(int(math.ceil(n**(1/3))),int(math.ceil(5*np.nanmax(g))) if np.isfinite(g).any() else 1)
    result={'method':'joint_studentized_batch_bootstrap','confidence':confidence,
            'replicates':int(replicates),'seed':seed,'block_length_frames':length,
            'block_length_source':'user' if block_length is not None else 'max(ceil(N**(1/3)),ceil(5*max(statistical_inefficiency)))',
            'approximate_blocks':n//length,'minimum_blocks':min_blocks,
            'frame_count':n,'complete_coverage_columns':int(eligible.sum()),
            'statistical_inefficiency':[float(v) if np.isfinite(v) else None for v in g],
            'effective_frames':[float(n/v) if np.isfinite(v) else None for v in g],
            'pointwise_lower':[None]*m,'pointwise_upper':[None]*m,
            'simultaneous_lower':[None]*m,'simultaneous_upper':[None]*m,
            'assumptions':'stationary representative trajectory; chosen block length captures serial dependence',
            'systematic_errors_included':False,'stationarity_validated':False,
            'coverage_rule':'intervals require every original frame at a coordinate'}
    if n<20 or n//length<min_blocks:
        result['status']='insufficient_independent_blocks';return result
    if not eligible.any():
        result['status']='no_complete_coverage';return result
    y=x[:,eligible];mean=y.mean(axis=0)
    batches=np.array_split(y,n//length)
    means=np.asarray([batch.mean(axis=0) for batch in batches]);k=len(means)
    center=means.mean(axis=0);se=means.std(axis=0,ddof=1)/np.sqrt(k)
    # Guard against residual positive dependence between adjacent batch means.
    batch_g=np.asarray([statistical_inefficiency(means[:,j]) for j in range(means.shape[1])])
    effective_batches=k/batch_g
    admissible=effective_batches>=min_blocks
    result['residual_batch_inefficiency']=batch_g.tolist()
    result['effective_batch_count']=effective_batches.tolist()
    if not admissible.any():
        result['status']='insufficient_effective_batches';return result
    variable=se>np.finfo(float).eps*np.maximum(1,np.abs(mean))
    rng=np.random.default_rng(seed);student=np.zeros((int(replicates),y.shape[1]))
    degenerate=np.zeros(y.shape[1],dtype=int)
    for first in range(0,int(replicates),16):
        count=min(16,int(replicates)-first)
        sampled=means[rng.integers(0,k,size=(count,k))]
        boot_mean=sampled.mean(axis=1);boot_se=sampled.std(axis=1,ddof=1)/np.sqrt(k)
        valid=boot_se>np.finfo(float).eps*np.maximum(1,np.abs(mean))
        degenerate+=np.sum(~valid & variable,axis=0)
        student[first:first+count]=np.divide(boot_mean-center,boot_se,out=np.zeros_like(boot_mean),where=valid & variable)
    alpha=(1-confidence)/2
    if np.any(degenerate>alpha*replicates):
        result['status']='insufficient_batch_variation';return result
    lo,hi=np.quantile(student,[alpha,1-alpha],axis=0)
    from scipy.stats import t
    floor=t.ppf(1-alpha,np.maximum(effective_batches-1,1))
    lo=np.minimum(lo,-floor);hi=np.maximum(hi,floor)
    se=se*np.sqrt(batch_g)
    lower,upper=mean-hi*se,mean-lo*se
    critical=float(np.quantile(np.max(np.abs(student[:,variable]),axis=1),confidence)) if variable.any() else 0.
    for target,values in [('pointwise_lower',lower),('pointwise_upper',upper),
                          ('simultaneous_lower',mean-critical*se),('simultaneous_upper',mean+critical*se)]:
        for j,value,valid in zip(np.flatnonzero(eligible),values,admissible):
            if valid:result[target][int(j)]=float(value)
    family_size=max(1,int(np.sum(variable & admissible)))
    critical=max(critical,float(t.ppf(1-(1-confidence)/(2*family_size),np.maximum(effective_batches[admissible]-1,1)).max()))
    # Reapply the critical-value floor to the simultaneous band too.
    for key,values in [('simultaneous_lower',mean-critical*se),('simultaneous_upper',mean+critical*se)]:
        for j,value,valid in zip(np.flatnonzero(eligible),values,admissible):
            if valid:result[key][int(j)]=float(value)
    result['batch_lengths_frames']=[len(batch) for batch in batches]
    result['batch_count']=k
    result['standard_error_method']='batch-mean standard error, inflated for residual positive batch autocorrelation; Student-t pointwise floor and conservative Bonferroni-t simultaneous floor'
    result.update(status='estimated',simultaneous_family_size=int(admissible.sum()),
                  estimated_column_indices=np.flatnonzero(eligible)[admissible].tolist(),
                  simultaneous_critical_value=critical,
                  half_mean_difference=(y[n//2:].mean(axis=0)-y[:n//2].mean(axis=0)).tolist())
    return result


# ---------------------------------------------------------------------------
# Effective sample size, blocking and moving-block bootstrap of one series.
# ---------------------------------------------------------------------------

#: Methods behind :func:`sampling_uncertainty`, recorded with every result.
SAMPLING_METHODS = {
    "effective_sample_size": (
        "N / tau_int with the integrated autocorrelation time tau_int(M) = 1 + 2 * sum_{t=1..M} rho(t) "
        "(FFT autocorrelation) and Sokal's automatic window: the smallest M with M >= 5 * tau_int(M) "
        "(A. D. Sokal, 'Monte Carlo methods in statistical mechanics: foundations and new algorithms', "
        "Cargese lecture notes 1996, publ. 1997; same estimator as emcee.autocorr.integrated_time)"),
    "blocking": (
        "Flyvbjerg-Petersen blocking (H. Flyvbjerg and H. G. Petersen, J. Chem. Phys. 91, 461 (1989)): "
        "repeated pairwise averaging; the level is chosen automatically with the lag-one test of "
        "M. Jonsson, Phys. Rev. E 98, 043304 (2018), as in the author's reference implementation "
        "(github.com/computative/block, python/tictoc.py: first level k with M_k < chi2_0.99(k+1)); "
        "an odd trailing frame is dropped at each level instead of truncating to a power of two"),
    "block_bootstrap": (
        "moving-block bootstrap of the mean (H. R. Kuensch, Ann. Statist. 17, 1217 (1989)): ceil(N/L) "
        "overlapping blocks of L frames drawn with replacement; percentile interval and bootstrap SD; "
        "L = ceil(5 * tau_int) unless given, and at least 8 blocks (N/L >= 8) are required"),
}


def _autocorrelation(x):
    """Normalised autocorrelation of a centred series by zero-padded FFT."""
    import numpy as np
    n = len(x)
    size = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(x, n=size)
    ac = np.fft.irfft(f * f.conj(), n=size)[:n]
    return ac / ac[0]


def _finite_series(values):
    import numpy as np
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError("A finite one-dimensional series of at least two values is required")
    return x


def integrated_autocorrelation_time(values, *, window_constant=5.0):
    """Integrated autocorrelation time with Sokal's automatic window.

    ``tau_int(M) = 1 + 2 * sum_{t=1}^{M} rho(t)`` (so ``tau_int = 1`` for
    uncorrelated frames and the effective sample size is ``N / tau_int``); the
    window ``M`` is the smallest lag with ``M >= window_constant * tau_int(M)``
    (Sokal 1997; the estimator of ``emcee.autocorr.integrated_time``).

    Parameters
    ----------
    values : sequence of float
        Finite series of at least two values, in time order.
    window_constant : float, default 5.0
        Sokal's ``c``.

    Returns
    -------
    dict
        ``tau_int_frames`` (clipped to ``[1, N]``), ``window_frames`` ``M``,
        ``window_converged`` (``False`` when no lag satisfied the rule, so the
        last lag was used and ``tau_int`` is a lower bound) and
        ``effective_sample_size`` ``N / tau_int``. A constant series gives
        ``tau_int = 1``.

    Raises
    ------
    ValueError
        For a non-finite or too short series.

    Examples
    --------
    >>> from crevice.uncertainty import integrated_autocorrelation_time
    >>> integrated_autocorrelation_time([1.0, -1.0] * 50)["tau_int_frames"]
    1.0
    """
    import numpy as np
    x = _finite_series(values)
    n = len(x)
    x = x - x.mean()
    if np.dot(x, x) <= np.finfo(float).eps * n:
        return {"tau_int_frames": 1.0, "window_frames": 0, "window_converged": True, "effective_sample_size": float(n)}
    tau = 2.0 * np.cumsum(_autocorrelation(x)) - 1.0
    lags = np.arange(n)
    passed = lags >= window_constant * tau
    window = int(np.argmax(passed)) if passed.any() else n - 1
    value = float(min(n, max(1.0, tau[window])))
    return {"tau_int_frames": value, "window_frames": window, "window_converged": bool(passed.any()),
            "effective_sample_size": n / value}


def blocking_standard_error(values):
    """Flyvbjerg–Petersen blocking estimate of the standard error of the mean.

    The series is repeatedly averaged in adjacent pairs (block size ``2**k``
    frames at level ``k``). The level is chosen automatically with the lag-one
    autocovariance test of Jonsson (2018) as in his reference implementation:
    ``M_k = sum_{j>=k} n_j (gamma_j / s_j)**2`` and the first level with
    ``M_k < chi2.ppf(0.99, k + 1)``. The standard error is ``sqrt(s_k / n_k)``
    with the population variance ``s_k`` of the ``n_k`` block means.

    Parameters
    ----------
    values : sequence of float
        Finite series of at least two values, in time order.

    Returns
    -------
    dict
        ``standard_error``, ``block_size_frames`` (``2**k``), ``level`` ``k``,
        ``levels`` and ``converged`` (``False`` when the test never accepted a
        level before the last one: more data are needed and the standard error
        is likely an underestimate).

    Raises
    ------
    ValueError
        For a non-finite or too short series.

    Examples
    --------
    >>> from crevice.uncertainty import blocking_standard_error
    >>> r = blocking_standard_error([0.0, 1.0] * 32)   # alternating: pairs average to a constant
    >>> r["block_size_frames"], r["standard_error"]
    (2, 0.0)
    """
    import numpy as np
    from scipy.stats import chi2
    x = _finite_series(values)
    mu = x.mean()
    variances, gammas, sizes = [], [], []
    while len(x) >= 2:
        n = len(x)
        gammas.append(float(np.sum((x[:-1] - mu) * (x[1:] - mu)) / n))
        variances.append(float(x.var()))
        sizes.append(n)
        if n % 2:
            x = x[:-1]
        x = 0.5 * (x[0::2] + x[1::2])
    s, g, n_k = np.asarray(variances), np.asarray(gammas), np.asarray(sizes, dtype=float)
    levels = len(s)
    if s[0] <= np.finfo(float).eps * max(1.0, abs(mu)):
        return {"standard_error": 0.0, "block_size_frames": 1, "level": 0, "levels": levels, "converged": True}
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(s > 0, n_k * (g / s) ** 2, 0.0)
    statistic = np.cumsum(terms[::-1])[::-1]
    level = levels - 1
    for k in range(levels):
        if statistic[k] < chi2.ppf(0.99, k + 1):
            level = k
            break
    return {"standard_error": float(np.sqrt(s[level] / n_k[level])), "block_size_frames": int(2 ** level),
            "level": int(level), "levels": int(levels), "converged": bool(level < levels - 1)}


def moving_block_bootstrap(values, *, block_length, confidence=0.95, replicates=2000, seed=20261020):
    """Moving-block bootstrap of the mean of a serially correlated series.

    ``ceil(N / L)`` of the ``N - L + 1`` overlapping blocks of ``L`` frames are
    drawn with replacement; each replicate mean is the mean of the drawn block
    means (Künsch 1989).

    Parameters
    ----------
    values : sequence of float
        Finite series, in time order.
    block_length : int
        Block length ``L`` in frames, ``1 <= L <= N``.
    confidence : float, default 0.95
        Level of the two-sided percentile interval.
    replicates : int, default 2000
        Bootstrap replicates.
    seed : int, default 20261020
        Seed of the NumPy generator (reproducible).

    Returns
    -------
    dict
        ``standard_error`` (SD of the replicate means), ``lower``, ``upper``,
        ``block_length_frames``, ``block_count`` and ``replicates``.

    Examples
    --------
    >>> from crevice.uncertainty import moving_block_bootstrap
    >>> r = moving_block_bootstrap([1.0, 2.0, 3.0, 4.0] * 10, block_length=4)
    >>> r["block_count"], round(r["standard_error"], 3)
    (10, 0.0)
    """
    import numpy as np
    x = _finite_series(values)
    n = len(x)
    length = int(block_length)
    if not 1 <= length <= n:
        raise ValueError("block_length must lie between 1 and the series length")
    if not 0 < confidence < 1:
        raise ValueError("confidence must lie strictly between 0 and 1")
    cumulative = np.concatenate(([0.0], np.cumsum(x)))
    block_means = (cumulative[length:] - cumulative[:-length]) / length
    count = int(math.ceil(n / length))
    rng = np.random.default_rng(seed)
    means = np.empty(int(replicates))
    for first in range(0, int(replicates), 512):
        stop = min(int(replicates), first + 512)
        means[first:stop] = block_means[rng.integers(0, len(block_means), size=(stop - first, count))].mean(axis=1)
    alpha = (1 - confidence) / 2
    lower, upper = np.quantile(means, [alpha, 1 - alpha])
    return {"standard_error": float(means.std(ddof=1)), "lower": float(lower), "upper": float(upper),
            "block_length_frames": length, "block_count": count, "replicates": int(replicates)}


#: Fields of :func:`sampling_uncertainty` repeated as columns by the summary tables.
SAMPLING_FIELDS = (
    "status", "effective_sample_size", "integrated_autocorrelation_time_frames", "autocorrelation_window_frames",
    "blocking_standard_error", "blocking_block_size_frames", "block_bootstrap_status",
    "block_bootstrap_lower", "block_bootstrap_upper", "block_bootstrap_standard_error",
    "block_bootstrap_block_length_frames",
)


def sampling_uncertainty(values, *, confidence=0.95, replicates=2000, seed=20261020, block_length=None,
                         regular=True, min_frames=20, min_blocks=8):
    """Effective sample size, blocking standard error and block-bootstrap interval of a frame series.

    Additive companion of the batch-bootstrap ``confidence_interval``
    reported by :func:`crevice.cavity_trajectory.describe_series` (which is
    unchanged): the methods are listed in :data:`SAMPLING_METHODS`.

    Parameters
    ----------
    values : sequence of float
        One value per frame in time order; NaN marks a missing frame.
    confidence : float or None, default 0.95
        Level of the bootstrap interval; ``None`` or 0 skips the bootstrap
        (effective sample size and blocking are still reported).
    replicates : int, default 2000
        Bootstrap replicates.
    seed : int, default 20261020
        Bootstrap seed.
    block_length : int, optional
        Bootstrap block length in frames; default ``ceil(5 * tau_int)``.
    regular : bool, default True
        Whether frames are evenly spaced in time; ``False`` withholds every estimate.
    min_frames : int, default 20
        Fewer frames withhold every estimate.
    min_blocks : int, default 8
        The bootstrap needs ``N / L >= min_blocks``.

    Returns
    -------
    dict
        ``status`` (``"estimated"``, ``"unavailable_missing_frames"``,
        ``"unavailable_irregular_timestamps"``, ``"unavailable_too_few_frames"``,
        ``"unavailable_observed_constant"``), ``frame_count``,
        ``effective_sample_size``, ``integrated_autocorrelation_time_frames``,
        ``autocorrelation_window_frames``, ``autocorrelation_window_converged``,
        ``blocking_standard_error``, ``blocking_block_size_frames``,
        ``blocking_converged``, ``block_bootstrap_status`` (``"estimated"``,
        ``"insufficient_blocks"``, ``"not_requested"`` or the overall
        status), ``block_bootstrap_lower``/``block_bootstrap_upper`` (bounds
        on the mean), ``block_bootstrap_standard_error``,
        ``block_bootstrap_block_length_frames``,
        ``block_bootstrap_block_length_source`` (``"ceil(5*tau_int)"`` or
        ``"user"``), ``block_bootstrap_block_count``, ``confidence``,
        ``replicates``, ``seed`` and ``methods``. Estimates absent for a status
        are ``None``.

    Notes
    -----
    All three estimates assume a stationary series and are biased low when
    the series is short compared with its correlation time (see
    the AR(1) table in the uncertainty method note); a
    non-converged window or blocking level flags that case. They describe
    sampling uncertainty of the mean only, not force-field or model error.

    Examples
    --------
    >>> import numpy as np
    >>> from crevice.uncertainty import sampling_uncertainty
    >>> r = sampling_uncertainty(np.random.default_rng(0).normal(size=400))
    >>> r["status"], r["block_bootstrap_status"], round(r["integrated_autocorrelation_time_frames"]) == 1
    ('estimated', 'estimated', True)
    """
    import numpy as np
    x = np.asarray(values, dtype=float)
    n = int(len(x))
    result = {"status": "estimated", "frame_count": n, "effective_sample_size": None,
              "integrated_autocorrelation_time_frames": None, "autocorrelation_window_frames": None,
              "autocorrelation_window_converged": None, "blocking_standard_error": None,
              "blocking_block_size_frames": None, "blocking_converged": None,
              "block_bootstrap_status": None, "block_bootstrap_lower": None, "block_bootstrap_upper": None,
              "block_bootstrap_standard_error": None, "block_bootstrap_block_length_frames": None,
              "block_bootstrap_block_length_source": None, "block_bootstrap_block_count": None,
              "confidence": confidence or None, "replicates": int(replicates), "seed": seed,
              "methods": SAMPLING_METHODS,
              "assumptions": "stationary representative trajectory; frames evenly spaced in time"}
    if x.ndim != 1 or n == 0 or np.isinf(x).any():
        raise ValueError("values must be a one-dimensional series without infinities")
    if not np.isfinite(x).all():
        status = "unavailable_missing_frames"
    elif not regular:
        status = "unavailable_irregular_timestamps"
    elif n < min_frames:
        status = "unavailable_too_few_frames"
    elif np.ptp(x) == 0:
        status = "unavailable_observed_constant"
    else:
        status = None
    if status is not None:
        result.update(status=status, block_bootstrap_status=status)
        return result
    tau = integrated_autocorrelation_time(x)
    blocking = blocking_standard_error(x)
    result.update(effective_sample_size=tau["effective_sample_size"],
                  integrated_autocorrelation_time_frames=tau["tau_int_frames"],
                  autocorrelation_window_frames=tau["window_frames"],
                  autocorrelation_window_converged=tau["window_converged"],
                  blocking_standard_error=blocking["standard_error"],
                  blocking_block_size_frames=blocking["block_size_frames"],
                  blocking_converged=blocking["converged"])
    if not confidence:
        result["block_bootstrap_status"] = "not_requested"
        return result
    if block_length is not None:
        length, source = int(block_length), "user"
    else:
        length, source = max(1, int(math.ceil(5 * tau["tau_int_frames"]))), "ceil(5*tau_int)"
    result.update(block_bootstrap_block_length_frames=length, block_bootstrap_block_length_source=source)
    if length > n or n / length < min_blocks:
        result["block_bootstrap_status"] = "insufficient_blocks"
        return result
    boot = moving_block_bootstrap(x, block_length=length, confidence=confidence, replicates=replicates, seed=seed)
    result.update(block_bootstrap_status="estimated", block_bootstrap_lower=boot["lower"],
                  block_bootstrap_upper=boot["upper"], block_bootstrap_standard_error=boot["standard_error"],
                  block_bootstrap_block_count=boot["block_count"])
    return result


def trajectory_frame_statistics(analysis, *, confidence=0.95, block_length=None, replicates=2000, seed=20260913):
    """Per-frame summary statistics of a :func:`crevice.trajectory.analyze_trajectory` result.

    Every per-frame quantity of the analyses that were run (the numeric
    columns of ``*_frames.csv``: ``min_radius_A``, ``mean_radius_A``,
    ``max_radius_A`` and ``bottleneck_axis_position_A`` for profiles;
    ``residue_contact_count``; ``cavity_count``; ``tunnel_count``;
    ``network_node_count`` and ``network_edge_count``) is described with
    :func:`crevice.cavity_trajectory.describe_series`: frame quantiles, the
    batch-bootstrap ``confidence_interval`` and the ``sampling_uncertainty``
    record (effective sample size, blocking standard error and moving-block
    bootstrap interval with the block sizes used). Frames without a profile are
    missing values, so profile quantities of such a trajectory get coverage
    only and no interval.

    Parameters
    ----------
    analysis : TrajectoryAnalysis
        Result of :func:`crevice.trajectory.analyze_trajectory`.
    confidence : float, default 0.95
        Interval level; 0 or ``None`` reports descriptive statistics,
        effective sample sizes and blocking only (no intervals).
    block_length : int, optional
        Bootstrap block length in frames (both bootstraps); automatic by default.
    replicates : int, default 2000
        Bootstrap replicates.
    seed : int, default 20260913
        Bootstrap seed.

    Returns
    -------
    dict
        ``quantities`` (name to statistics record), ``frame_count``,
        ``time_ps`` (first, last), ``regular_sampling``, ``analyses``,
        ``definition`` and ``methods``. Flatten with
        :func:`crevice.io.summary_statistics_rows`.
    """
    import numpy as np
    from .cavity_trajectory import describe_series
    frames = list(analysis.frames)
    n = len(frames)
    times = [f.time for f in frames]
    known = [t for t in times if t is not None]
    regular = len(known) == n and (n < 3 or bool(np.allclose(np.diff(known), np.diff(known)[0], rtol=1e-7, atol=1e-7)))
    regular = regular and (n < 2 or bool(np.all(np.diff(known) > 0)))

    def series(get):
        return np.asarray([np.nan if (value := get(f)) is None or value == "" else float(value) for f in frames])

    analyses = set(analysis.analyses)
    wanted = {}
    if "profile" in analyses:
        wanted.update(min_radius_A=lambda f: f.profile.min_radius if f.profile else None,
                      mean_radius_A=lambda f: f.profile.mean_radius if f.profile else None,
                      max_radius_A=lambda f: f.profile.max_radius if f.profile else None,
                      bottleneck_axis_position_A=lambda f: f.profile.bottleneck.t if f.profile else None)
    if "residues" in analyses:
        wanted["residue_contact_count"] = lambda f: len(f.residue_contacts)
    if "cavities" in analyses:
        wanted["cavity_count"] = lambda f: len(f.cavities)
    if "tunnels" in analyses:
        wanted["tunnel_count"] = lambda f: len(f.tunnels)
    if "network" in analyses:
        wanted["network_node_count"] = lambda f: len(f.network.nodes) if f.network else (f.metadata or {}).get("network_node_count")
        wanted["network_edge_count"] = lambda f: len(f.network.edges) if f.network else (f.metadata or {}).get("network_edge_count")
    options = dict(confidence=confidence or 0.95, block_length=block_length, replicates=replicates, seed=seed,
                   regular=regular and bool(confidence))
    quantities = {name: describe_series(series(get), **options) for name, get in wanted.items()}
    return {"quantities": quantities, "frame_count": n,
            "time_ps": [known[0], known[-1]] if known else None, "regular_sampling": regular,
            "analyses": sorted(analyses), "confidence": confidence or None,
            "definition": "statistics over analysed frames in time order; missing values are frames without the "
                          "quantity (for example an unresolved profile)",
            "methods": SAMPLING_METHODS,
            "limitations": "sampling uncertainty of the mean for a stationary trajectory only; not force-field, "
                           "model or geometric-definition error (see --sensitivity for definitions)"}
