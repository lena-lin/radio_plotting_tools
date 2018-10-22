from skimage.feature import blob_log


def find_blobs(clean_map, n_blobs_max, min_sigma=1, max_sigma=20, num_sigma=10, threshold=0.001, overlap=0.7):
    blobs = blob_log(
                    clean_map/abs(clean_map).max(),
                    min_sigma,
                    max_sigma,
                    num_sigma,
                    threshold,
                    overlap
                    )

    if len(blobs) <= n_blobs_max:
        return blobs

    else:
        return find_blobs(clean_map,n_blobs_max, min_sigma=min_sigma+2, max_sigma=max_sigma+10, num_sigma=5, threshold=threshold+0.001)
