# vim: set ft=python tw=88 nu ai et ts=4 sw=4:
# ==============================================================================
#  Two-sample F-type bootstrap test
#  Python equivalent of Fstattwosample.R (lecturer's helper)
# ==============================================================================

import numpy as np
from scipy.stats import f as f_dist


def f_stat_two_sample(x, y, t_seq, alpha=0.05, method=1, replications=100,
                      random_state=None):
    """
    Two-sample F-type test for functional data.

    Parameters
    ----------
    x : np.ndarray, shape (n_timepoints, n_x)
        Evaluated curves for group X.
    y : np.ndarray, shape (n_timepoints, n_y)
        Evaluated curves for group Y.
    t_seq : np.ndarray, shape (n_timepoints,)
        Time points at which curves are evaluated.
    alpha : float
        Significance level (default 0.05).
    method : int
        1 = T-statistic (naive), 2 = Bootstrap.
    replications : int
        Number of bootstrap replications (used when method=2).
    random_state : int or None
        Random seed for reproducibility.

    Returns
    -------
    dict with keys:
        'statistics' : float  — observed F statistic
        'pvalue'     : float  — p-value
        'params'     : dict   — df1/df2 (method=1) or bootstrap dist (method=2)
    """
    if random_state is not None:
        np.random.seed(random_state)

    n = x.shape[1]
    m = y.shape[1]
    k = len(t_seq)

    cn = (n * m) / (n + m)

    mu_x = x.mean(axis=1, keepdims=True)   # (k, 1)
    mu_y = y.mean(axis=1, keepdims=True)

    delta_t = mu_x - mu_y                  # (k, 1)

    z_x = x - mu_x                         # centered X
    z_y = y - mu_y                         # centered Y
    z_t = np.hstack([z_x, z_y])            # (k, n+m)

    if n > k or m > k:
        Sigma = (z_t.T @ z_t) / (n - 2)   # (n+m, n+m)
    else:
        Sigma = (z_t @ z_t.T) / (n - 2)   # (k, k)

    A = np.trace(Sigma)
    B = np.trace(Sigma @ Sigma)

    Fstat = float((cn * (delta_t.T @ delta_t).squeeze()) / A)

    if method == 1:
        kappa  = A**2 / B
        pvalue = 1 - f_dist.cdf(Fstat, kappa, (n - 2) * kappa)
        params = {"df1": kappa, "df2": (n - 2) * kappa}

    elif method == 2:
        bt_fstat = np.zeros(replications)

        for i in range(replications):
            rep1       = np.random.choice(n, n, replace=True)
            rep2       = np.random.choice(m, m, replace=True)
            x_star     = x[:, rep1]
            y_star     = y[:, rep2]

            mu_x_star  = x_star.mean(axis=1, keepdims=True)
            mu_y_star  = y_star.mean(axis=1, keepdims=True)
            delta_star = mu_x_star - mu_y_star

            bz_x = x_star - mu_x_star
            bz_y = y_star - mu_y_star
            z_star = np.hstack([bz_x, bz_y])

            if n > k or m > k:
                bt_sigma = (z_star.T @ z_star) / (n - 2)
            else:
                bt_sigma = (z_star @ z_star.T) / (n - 2)

            bt_mu       = (delta_star - delta_t)
            bt_fstat[i] = float(((cn * bt_mu.T @ bt_mu) / np.trace(bt_sigma)).squeeze())

        pvalue = float(np.mean(bt_fstat >= Fstat))
        params = {"bt_fstat": bt_fstat}

    else:
        raise ValueError(f"method must be 1 or 2, got {method}")

    return {"statistics": Fstat, "pvalue": pvalue, "params": params}
