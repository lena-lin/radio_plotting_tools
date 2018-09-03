import numpy as np


def noise_level(clean_map, window_width=100, n_sigma=5):

    y_max = clean_map.shape[0]
    x_max = clean_map.shape[1]
    x_1 = 0 + window_width
    x_2 = x_max - window_width
    y_1 = 0 + window_width
    y_2 = y_max - window_width

    rms = np.array(
           [clean_map[0:y_max, 0:x_1].std(),
            clean_map[0:y_max, x_2:x_max].std(),
            clean_map[0:y_1, 0:x_max].std(),
            clean_map[y_2:y_max, 0:x_max].std()]
        ).mean()

    noise = n_sigma * rms

    return noise
