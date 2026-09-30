"""Shared ERP display statistic."""
import numpy as np

def mean_ci(data):
    """Pointwise 95% t intervals across participant means, not across trials."""
    from scipy.stats import t
    mean = np.mean(data, axis=0)
    half_width = t.ppf(.975, len(data)-1) * np.std(data, axis=0, ddof=1) / np.sqrt(len(data))
    return np.array([mean-half_width, mean+half_width])

