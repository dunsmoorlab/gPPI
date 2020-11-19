from fg_config import *
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
def plot_full_skipped_corr(x,y,title,xlab,ylab):
    from pingouin.utils import _is_sklearn_installed
    _is_sklearn_installed(raise_error=True)
    from scipy.stats import chi2
    from sklearn.covariance import MinCovDet
    X = np.column_stack((x, y))
    nrows, ncols = X.shape
    gval = np.sqrt(chi2.ppf(0.975, 2))
    # Compute center and distance to center
    center = MinCovDet(random_state=42).fit(X).location_
    B = X - center
    B2 = B**2
    bot = B2.sum(axis=1)
    # Loop over rows
    dis = np.zeros(shape=(nrows, nrows))
    for i in np.arange(nrows):
        if bot[i] != 0:  # Avoid division by zero error
            dis[i, :] = np.linalg.norm(B * B2[i, :] / bot[i], axis=1)
    def idealf(x):
        """Compute the ideal fourths IQR (Wilcox 2012).
        """
        n = len(x)
        j = int(np.floor(n / 4 + 5 / 12))
        y = np.sort(x)
        g = (n / 4) - j + (5 / 12)
        low = (1 - g) * y[j - 1] + g * y[j]
        k = n - j + 1
        up = (1 - g) * y[k - 1] + g * y[k - 2]
        return up - low

    # One can either use the MAD or the IQR (see Wilcox 2012)
    # MAD = mad(dis, axis=1)
    iqr = np.apply_along_axis(idealf, 1, dis)
    thresh = (np.median(dis, axis=1) + gval * iqr)
    outliers = np.apply_along_axis(np.greater, 0, dis, thresh).any(axis=0)

    cloud = X[~outliers]
    R = np.random.RandomState(42)
    rs = np.zeros(10000)
    for i in range(10000):
        # _samp = np.random.choice(range(len(cloud)),size=len(cloud))
        _samp = R.choice(range(len(cloud)),size=len(cloud))
        rs[i] = pearsonr(cloud[_samp,0],cloud[_samp,1])[0]
    if rs.mean() > 0:
        p = (1 - np.mean(rs > 0)) * 2
    else:
        p = (1 - np.mean(rs < 0)) * 2

    r_pearson, _ = pearsonr(x[~outliers], y[~outliers])
    ci_l, ci_u = np.percentile(rs,[2.5,97.5])
    
    fig, (ax1, ax3) = plt.subplots(2, figsize=(6, 10))
    # plt.subplots_adjust(wspace=0.3)
    sns.despine()

    # Scatter plot and regression lines
    sns.regplot(x[~outliers], y[~outliers], ax=ax1, color='darkcyan')
    ax1.scatter(x[outliers], y[outliers], color='indianred', label='outliers')
    ax1.scatter(x[~outliers], y[~outliers], color='seagreen', label='good')

    sns.distplot(rs, kde=True, ax=ax3, color='steelblue')
    for i in [ci_l,ci_u]:
        ax3.axvline(x=i, color='coral', lw=2)
    ax3.axvline(x=0, color='k', ls='--', lw=1.5)
    ax3.set_xlabel('Correlation coefficient')
    ax3.set_title(
        'Skipped Pearson r = {}\n95% CI = [{}, {}], P = {}'.format(r_pearson.round(2),
                                                           ci_l.round(2),
                                                           ci_u.round(2),
                                                           p.round(4)),
        y=1.05)
    ax1.set_xlim([i*1.2 for i in ax1.get_xlim()])
    ax1.set_title(title)
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)
    # Optimize layout
    plt.tight_layout()

def skipped_corr(x, y, vis=False, ax=None, color='blue', return_dist=False):

    from pingouin.utils import _is_sklearn_installed
    _is_sklearn_installed(raise_error=True)
    from scipy.stats import chi2
    from sklearn.covariance import MinCovDet
    X = np.column_stack((x, y))
    nrows, ncols = X.shape
    gval = np.sqrt(chi2.ppf(0.975, 2))
    # Compute center and distance to center
    center = MinCovDet(random_state=42).fit(X).location_
    B = X - center
    B2 = B**2
    bot = B2.sum(axis=1)
    # Loop over rows
    dis = np.zeros(shape=(nrows, nrows))
    for i in np.arange(nrows):
        if bot[i] != 0:  # Avoid division by zero error
            dis[i, :] = np.linalg.norm(B * B2[i, :] / bot[i], axis=1)
    def idealf(x):
        """Compute the ideal fourths IQR (Wilcox 2012).
        """
        n = len(x)
        j = int(np.floor(n / 4 + 5 / 12))
        y = np.sort(x)
        g = (n / 4) - j + (5 / 12)
        low = (1 - g) * y[j - 1] + g * y[j]
        k = n - j + 1
        up = (1 - g) * y[k - 1] + g * y[k - 2]
        return up - low

    # One can either use the MAD or the IQR (see Wilcox 2012)
    # MAD = mad(dis, axis=1)
    iqr = np.apply_along_axis(idealf, 1, dis)
    thresh = (np.median(dis, axis=1) + gval * iqr)
    outliers = np.apply_along_axis(np.greater, 0, dis, thresh).any(axis=0)

    cloud = X[~outliers]

    R = np.random.RandomState(42)
    rs = np.zeros(10000)
    for i in range(10000):
        # _samp = np.random.choice(range(len(cloud)),size=len(cloud))
        _samp = R.choice(range(len(cloud)),size=len(cloud))
        rs[i] = pearsonr(cloud[_samp,0],cloud[_samp,1])[0]
    if rs.mean() > 0:
        p = (1 - np.mean(rs > 0)) * 2
    else:
        p = (1 - np.mean(rs < 0)) * 2
    r_pearson, _ = pearsonr(x[~outliers], y[~outliers])
    ci_l, ci_u = np.percentile(rs,[2.5,97.5])
    
    # Scatter plot and regression lines
    if vis and ax == None:
        fig, ax = plt.subplots()
    if vis:
        sns.regplot(x[~outliers], y[~outliers], ax=ax, color=color, scatter=False)
        ax.scatter(x, y, color=color, edgecolor=color)
    print(
            'Skipped Pearson r = {}\n95% CI = [{}, {}], P = {}'.format(r_pearson.round(2),
                                                           ci_l.round(2),
                                                           ci_u.round(2),
                                                           p.round(4)))
    if return_dist:
        return rs