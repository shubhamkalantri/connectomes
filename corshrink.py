import numpy as np


def corshrink(data: np.ndarray, lamb: float) -> np.ndarray:
    n = data.shape[0]
    p = data.shape[1]

    w = np.ones((n, 1)) * 1 / n

    h1 = 1 / (1 - w.T @ w)

    r0 = h1 * (data * np.sqrt(w)).T @ data * np.sqrt(w)

    power = (1 - lamb) * r0
    np.fill_diagonal(power, 1)

    return power
