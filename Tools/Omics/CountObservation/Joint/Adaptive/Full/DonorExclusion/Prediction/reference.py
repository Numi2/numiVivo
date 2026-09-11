"""Individual-cell SciPy reference, reused from Adaptive/verify.py; no native calls."""
import numpy as np
from scipy import optimize,special

def likelihood(y, library, phi):
    e = np.asarray(library, float) / 1000000.0
    y = np.asarray(y, float)
    total = float(y.sum())
    exposure = float(e.sum())
    if total == 0:
        mode = 0.0
    elif phi == 0:
        mode = total / exposure
    else:

        def score(log_rate):
            q = special.expit(np.log(phi * e) + log_rate)
            return float(np.sum(y * (1 - q) - q / phi))
        center = np.log(total / exposure)
        mode = float(np.exp(optimize.brentq(score, center - 50, center + 50, xtol=1e-13)))
    cache = {}

    def ell(rate):
        rate = float(rate)
        if rate in cache:
            return cache[rate]
        if rate == 0:
            value = 0.0 if total == 0 else -np.inf
        elif total == 0:
            value = -exposure * rate if phi == 0 else float(-np.log1p(phi * e * rate).sum() / phi)
        elif phi == 0:
            value = total * np.log(rate / mode) - exposure * (rate - mode)
        else:
            value = total * np.log(rate / mode) - float(np.dot(y + 1 / phi, np.log1p(phi * e * (rate - mode) / (1 + phi * e * mode))))
        cache[rate] = value
        return value

    def slope(rate, scale):
        if rate == 0:
            return -scale * exposure if total == 0 else np.inf
        if phi == 0:
            score = total - exposure * rate
        else:
            a = phi * e * rate
            score = float(np.sum(y / (1 + a) - a / (1 + a) / phi))
        return score * (1 + scale / rate)
    return (mode, ell, slope)

def moments(cr, tr, w):
    c = np.repeat(cr, len(tr))
    t = np.tile(tr, len(cr))
    mc = float(w @ c)
    mt = float(w @ t)
    lc = float(w @ np.log1p(c))
    lt = float(w @ np.log1p(t))
    return {'controlMeanCPM': mc, 'treatedMeanCPM': mt, 'controlVarianceCPM2': float(w @ (c - mc) ** 2), 'treatedVarianceCPM2': float(w @ (t - mt) ** 2), 'covarianceCPM2': float(w @ ((c - mc) * (t - mt))), 'responseMeanCPM': mt - mc, 'responseVarianceCPM2': float(w @ (t - c - (mt - mc)) ** 2), 'controlMeanLog1pCPM': lc, 'treatedMeanLog1pCPM': lt, 'log1pResponseMean': lt - lc, 'log1pResponseVariance': float(w @ (np.log1p(t) - np.log1p(c) - (lt - lc)) ** 2)}

def sampling(m, phi, libraries):
    e = np.array(libraries, float) / 1000000.0
    mean = float(e.sum() * m['treatedMeanCPM'])
    over = float(phi * (e @ e) * (m['treatedVarianceCPM2'] + m['treatedMeanCPM'] ** 2))
    latent = float(e.sum() ** 2 * m['treatedVarianceCPM2'])
    return dict(plannedCells=len(libraries), plannedLibraryCounts=sum(libraries), meanGeneCounts=mean, conditionalPoissonVariance=mean, conditionalCellOverdispersionVariance=over, latentRateVariance=latent, totalGeneCountVariance=mean + over + latent)
