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
