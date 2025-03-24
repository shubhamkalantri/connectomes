import numpy as np


def corshrink(data: np.ndarray, lamb: float) -> np.ndarray:
    n = data.shape[0]
    p = data.shape[1]

    data = data - np.mean(data, axis=0)
    data = data / np.std(data, axis=0, ddof=1)

    r = np.corrcoef(data, rowvar=False)

    vr = np.zeros((p, p))
    for i in range(p):
        for j in range(p):
            if i != j:
                vr[i, j] = (1 / (n - 1)) * ((1 - r[i, j] ** 2) ** 2)

    Rhat = (1 - lamb) * r
    np.fill_diagonal(Rhat, 1)

    return Rhat
