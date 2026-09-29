import numpy as np
from scipy import stats


def jeffreys(k, n, level=0.95):
    """Jeffreys interval for a binomial proportion."""
    if n == 0:
        return (np.nan, np.nan)
    a = (1 - level) / 2
    lo = 0.0 if k == 0 else stats.beta.ppf(a, k + 0.5, n - k + 0.5)
    hi = 1.0 if k == n else stats.beta.ppf(1 - a, k + 0.5, n - k + 0.5)
    return (float(lo), float(hi))


def bootstrap_mean_diff(a, b, draws=4000, seed=0, level=0.95):
    """Percentile bootstrap interval on mean(b) - mean(a)."""
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = np.empty(draws)
    for i in range(draws):
        d[i] = rng.choice(b, len(b)).mean() - rng.choice(a, len(a)).mean()
    q = (1 - level) / 2
    return float(b.mean() - a.mean()), float(np.quantile(d, q)), float(np.quantile(d, 1 - q))


def auc_ci(y, score, draws=10000, seed=0, level=0.95):
    """Percentile bootstrap interval on the area under the curve."""
    y, score = np.asarray(y, int), np.asarray(score, float)
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(draws):
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(y[i])) > 1:
            vals.append(auc(y[i], score[i])[0])
    q = (1 - level) / 2
    return float(np.quantile(vals, q)), float(np.quantile(vals, 1 - q))


def cluster_concentration(labels, clusters, draws=2000, seed=0):
    """Clusters touched by the annotated entries, against random reassignment."""
    labels = np.asarray(labels, bool)
    clusters = np.asarray(clusters)
    uniq, idx = np.unique(clusters, return_inverse=True)
    observed = len(np.unique(idx[labels]))
    rng = np.random.default_rng(seed)
    k = int(labels.sum())
    null = np.empty(draws, int)
    for i in range(draws):
        null[i] = len(np.unique(rng.choice(idx, k, replace=False)))
    return dict(n_clusters=len(uniq), n_annotated=k, observed=observed,
                expected=float(null.mean()),
                lo=float(np.percentile(null, 2.5)), hi=float(np.percentile(null, 97.5)),
                p=float((null <= observed).mean()))


def sign_test(deltas):
    """Two-sided binomial test that shifts are as often positive as negative."""
    d = np.asarray(deltas, float)
    pos, n = int((d > 0).sum()), int((d != 0).sum())
    return pos, n, float(stats.binomtest(pos, n).pvalue)


def auc(labels, scores):
    """Area under the ROC curve, with the equivalent Mann-Whitney p-value."""
    labels = np.asarray(labels, int)
    scores = np.asarray(scores, float)
    pos, neg = scores[labels == 1], scores[labels == 0]
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    return float(u.statistic / (len(pos) * len(neg))), float(u.pvalue), len(labels)


def auc_gap_permutation(labels, scores, group, draws=20000, seed=0):
    """Null for the difference in AUC between two groups, shuffling group membership."""
    labels, scores = np.asarray(labels, int), np.asarray(scores, float)
    group = np.asarray(group, bool)
    observed = auc(labels[group], scores[group])[0] - auc(labels[~group], scores[~group])[0]
    rng = np.random.default_rng(seed)
    null = np.empty(draws)
    for i in range(draws):
        s = rng.permutation(group)
        null[i] = auc(labels[s], scores[s])[0] - auc(labels[~s], scores[~s])[0]
    return float(observed), float((np.abs(null) >= abs(observed)).mean())
