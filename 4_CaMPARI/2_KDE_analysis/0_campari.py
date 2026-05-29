

import numpy as np
from sklearn.neighbors import NearestNeighbors


def p_threshold(p, p_cover=0.95):
    sp = np.sort(p)
    for i in range(1, sp.shape[0]):
        if sp[-i:].sum() > p_cover:
            return sp[-i-1]


def calc_xcorr_overlap(s1, s2, grid_offset=10.0, grid_spacing=5.0, p_cover=0.95):
    xyz1, w1, kde1 = s1
    xyz2, w2, kde2 = s2

    minima1 = xyz1.min(axis=0)
    maxima1 = xyz1.max(axis=0)

    minima2 = xyz2.min(axis=0)
    maxima2 = xyz2.max(axis=0)

    bounds = [[min(minima1[i], minima2[i]), max(maxima1[i], maxima2[i])]
              for i in range(len(minima1))]

    space = [np.arange(mini-grid_offset, maxi+grid_offset, grid_spacing)
             for mini, maxi in bounds]

    grid = np.meshgrid(*space)

    coords = np.vstack(list(map(np.ravel, grid))).T

    logp1 = kde1.score_samples(coords)
    logp2 = kde2.score_samples(coords)

    # regression coeff R
    p1 = np.exp(logp1)
    p2 = np.exp(logp2)

    R = np.corrcoef(p1, p2)[0, 1]

    # masked R
    dV = grid_spacing**3

    np1 = p1*dV
    np2 = p2*dV

    th1 = p_threshold(np1, p_cover)
    th2 = p_threshold(np2, p_cover)
    th = np.min([th1, th2])

    mask = (np1 > th) | (np2 > th)
    R_masked = np.corrcoef(np1[mask], np2[mask])[0, 1]

    # density overlap
    min_logp = np.minimum(logp1, logp2)
    min_p = np.exp(min_logp)

    S_overlap = min_p.sum()*dV

    return R, R_masked, S_overlap


def chamfer_like_distance(x, y, metric='l2', direction='bi', averaging='mean'):
    """Chamfer distance between two point clouds

    Parameters
    ----------
    x: numpy array [n_points_x, n_dims]
        first point cloud
    y: numpy array [n_points_y, n_dims]
        second point cloud
    metric: string or callable, default ‘l2’
        metric to use for distance computation. Any metric from scikit-learn or scipy.spatial.distance can be used.
    direction: str
        direction of Chamfer distance.
            'y_to_x':  computes average minimal distance from every point in y to x
            'x_to_y':  computes average minimal distance from every point in x to y
            'bi': compute both
    Returns
    -------
    chamfer_dist: float
        computed bidirectional Chamfer distance:
            sum_{x_i \in x}{\min_{y_j \in y}{||x_i-y_j||**2}} + sum_{y_j \in y}{\min_{x_i \in x}{||x_i-y_j||**2}}
    """

    if direction == 'y_to_x':
        x_nn = NearestNeighbors(
            n_neighbors=1, leaf_size=1, algorithm='kd_tree', metric=metric).fit(x)
        min_y_to_x = x_nn.kneighbors(y)[0]
        chamfer_dist = np.mean(min_y_to_x)
    elif direction == 'x_to_y':
        y_nn = NearestNeighbors(
            n_neighbors=1, leaf_size=1, algorithm='kd_tree', metric=metric).fit(y)
        min_x_to_y = y_nn.kneighbors(x)[0]
        chamfer_dist = np.mean(min_x_to_y)
    elif direction == 'bi':
        x_nn = NearestNeighbors(
            n_neighbors=1, leaf_size=1, algorithm='kd_tree', metric=metric).fit(x)
        min_y_to_x = x_nn.kneighbors(y)[0]
        min_y_to_x = x_nn.kneighbors(y)[0]
        y_nn = NearestNeighbors(
            n_neighbors=1, leaf_size=1, algorithm='kd_tree', metric=metric).fit(y)
        min_x_to_y = y_nn.kneighbors(x)[0]
        if averaging == 'mean':
            chamfer_dist = np.mean(min_y_to_x) + np.mean(min_x_to_y)
        elif averaging == 'median':
            chamfer_dist = np.median(min_y_to_x) + np.median(min_x_to_y)
        else:
            raise ValueError(
                "Invalid averaging type. Supported types: mean or median")
    else:
        raise ValueError(
            "Invalid direction type. Supported types: \'y_x\', \'x_y\', \'bi\'")

    return chamfer_dist
