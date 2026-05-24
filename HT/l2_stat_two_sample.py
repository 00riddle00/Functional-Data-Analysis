# vim: set ft=python tw=88 nu ai et ts=4 sw=4:
# ==============================================================================
#  Two-sample L2 norm-based bootstrap test
#  Python equivalent of L2stattwosample.R (lecturer's helper)
# ==============================================================================

import numpy as np
from scipy.stats import chi2


def l2_stat_two_sample(x, y, t_seq, alpha=0.05, method=1, replications=100,
                       random_state=None):
    """
    Two-sample L2 norm-based test for functional data.

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
        1 = chi-squared approximation (naive), 2 = Bootstrap.
    replications : int
        Number of bootstrap replications (used when method=2).
    random_state : int or None
        Random seed for reproducibility.

    Returns
    -------
    dict with keys:
        'statistics' : float  — observed L2 statistic
        'pvalue'     : float  — p-value
        'params'     : dict   — alpha/df (method=1) or bootstrap dist (method=2)
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

    z_x = x - mu_x
    z_y = y - mu_y
    z_t = np.hstack([z_x, z_y])            # (k, n+m)

    if n > k or m > k:
        Sigma = (z_t.T @ z_t) / (n - 2)
    else:
        Sigma = (z_t @ z_t.T) / (n - 2)

    A = np.trace(Sigma)
    B = np.trace(Sigma @ Sigma)

    L2stat = float(cn * (delta_t.T @ delta_t).squeeze())

    if method == 1:
        alp    = B / A
        df     = A**2 / B
        pvalue = 1 - chi2.cdf(L2stat / alp, df)
        params = {"alpha": alp, "df": df}

    elif method == 2:
        bt_l2stat = np.zeros(replications)

        for i in range(replications):
            rep1   = np.random.choice(n, n, replace=True)
            rep2   = np.random.choice(m, m, replace=True)
            x_star = x[:, rep1]
            y_star = y[:, rep2]

            mu_x_star  = x_star.mean(axis=1, keepdims=True)
            mu_y_star  = y_star.mean(axis=1, keepdims=True)
            delta_star = mu_x_star - mu_y_star

            bt_mu        = delta_star - delta_t
            bt_l2stat[i] = float((cn * bt_mu.T @ bt_mu).squeeze())

        pvalue = float(np.mean(bt_l2stat >= L2stat))
        params = {"boot_stat": bt_l2stat}

    else:
        raise ValueError(f"method must be 1 or 2, got {method}")

    return {"statistics": L2stat, "pvalue": pvalue, "params": params}
