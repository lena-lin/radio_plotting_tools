from skimage.feature import blob_log


def find_blobs(clean_map, max_sigma_=50, num_sigma_=1, threshold_=0.001):
    blobs = blob_log(
                    clean_map/abs(clean_map).max(),
                    max_sigma=max_sigma_,
                    num_sigma=num_sigma_,
                    threshold=threshold_,
                    )

    return blobs
