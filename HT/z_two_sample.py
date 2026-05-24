# vim: set ft=python tw=88 nu ai et ts=4 sw=4:
# ==============================================================================
#  Two-sample pointwise t-test
#  Python equivalent of Ztwosample.R (lecturer's helper)
# ==============================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t as t_dist


def z_two_sample(x, y, t_seq, alpha=0.05, plot=True, ax=None):
    """
    Two-sample pointwise t-test for functional data.

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
    plot : bool
        Whether to produce a plot (default True).
    ax : matplotlib Axes or None
        Axes to plot on. If None, a new figure is created.

    Returns
    -------
    dict with keys:
        'statistics_pointwise' : np.ndarray — pointwise Z statistics
        'params'               : dict       — critical value
    """
    n = x.shape[1]
    m = y.shape[1]
    k = len(t_seq)

    mu_x    = x.mean(axis=1, keepdims=True)   # (k, 1)
    mu_y    = y.mean(axis=1, keepdims=True)

    delta_t = (mu_x - mu_y).squeeze()          # (k,)

    z_x = x - mu_x
    z_y = y - mu_y
    z_t = np.hstack([z_x, z_y])               # (k, n+m)

    if n > k:
        Sigma = (z_t.T @ z_t) / (n - 2)
    else:
        Sigma = (z_t @ z_t.T) / (n - 2)

    gamma_t      = np.diag(Sigma)              # (k,)
    z_pointwise  = np.sqrt((n * m) / (n + m)) * delta_t / np.sqrt(gamma_t)

    crit     = t_dist.ppf(1 - alpha / 2, df=n - 2)
    crit_val = np.full(k, crit)

    if plot:
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 5))
        ylim_max = max(z_pointwise.max(),  crit) + 0.5
        ylim_min = min(z_pointwise.min(), -crit) - 0.5
        ax.plot(t_seq, z_pointwise, color='black', linewidth=1.5,
                label='Z statistic')
        ax.plot(t_seq,  crit_val, linestyle='--', linewidth=2,
                color='blue', label=f'±critical value (α={alpha})')
        ax.plot(t_seq, -crit_val, linestyle='--', linewidth=2, color='blue')
        ax.set_ylim(ylim_min, ylim_max)
        ax.set_xlabel('Time')
        ax.set_ylabel('Z statistics')
        ax.set_title('Two samples t-test')
        ax.legend()

    return {
        "statistics_pointwise": z_pointwise,
        "params": {"critical_value": crit}
    }
