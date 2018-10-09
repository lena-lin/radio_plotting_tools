from skimage.feature import blob_log


def find_blobs(clean_map, min_sigma=1, max_sigma_=50, num_sigma_=5, threshold_=0.001):
    blobs = blob_log(
                    clean_map/abs(clean_map).max(),
                    max_sigma=max_sigma_,
                    num_sigma=num_sigma_,
                    threshold=threshold_,
                    )

    if len(blobs) <= 5:
        return blobs

    else:
        print(len(blobs))
        return find_blobs(clean_map, min_sigma=min_sigma+0.5, max_sigma_=50, num_sigma_=5, threshold_=threshold_+0.001)
