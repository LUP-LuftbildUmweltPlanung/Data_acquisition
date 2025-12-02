import time

from encode_to_lmdb_parquet import create_or_open_lmdb
from safetensors.numpy import load
import download_by_shape_functions as func
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import itertools



def compute_mean_std(lmdb_path):
    """
    Calculate the mean and std in the given data (training data!!!)
    :param data: Dataset (for training)
    :return: Mean and std
    """
    db = create_or_open_lmdb(lmdb_path)

    mean = 0.0
    std = 0.0
    total_samples = 0

    counter = 0

    with db.begin() as txn:
        cursor = txn.cursor()
        start = time.time()
        for _, value in cursor:
            data_dict = load(value)

            curr_dict = data_dict["4"].astype(np.float32) / 255.0
            mean += curr_dict.sum()


            total_samples += data_dict["4"].size

            if counter % 1000 == 0:
                print(counter, time.time() - start)
                start = time.time()
            counter += 1

    mean /= total_samples

    counter = 0

    with db.begin() as txn:
        cursor = txn.cursor()
        start = time.time()
        for _, value in cursor:
            data_dict = load(value)

            curr_dict = data_dict["4"].astype(np.float32) / 255.0

            std += np.sum((curr_dict - mean)**2)
            if counter % 1000 == 0:
                print(counter, time.time() - start)
                start = time.time()
            counter += 1

    std /= total_samples
    return mean, np.sqrt(std)

def extract_sum_count(value):
    data_dict = load(value)
    array = data_dict["4"].astype(np.float32) / 255.0
    return array.sum(), array.size

def extract_squared_diff(value, mean):
    data_dict = load(value)
    array = data_dict["4"].astype(np.float32) / 255.0
    return np.sum((array - mean) ** 2)

def compute_mean_std_2pass(lmdb_path, num_workers=16):
    db = create_or_open_lmdb(lmdb_path)

    # === FIRST PASS: Compute mean ===
    total_sum = 0.0
    total_count = 0
    std_sum = 0.0
    counter = 0

    with db.begin() as txn:
        cursor = txn.cursor()
        to_process = []

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            start = time.time()
            for _, value in cursor:

                to_process.append(value)
                if len(to_process) >= 1000:

                    for s, c in list(executor.map(extract_sum_count, to_process)):
                        total_sum += s
                        total_count += c
                    to_process = []
                    print(counter , time.time() - start)
                    start = time.time()
                    counter += 1

            if to_process:
                for s, c in list(executor.map(extract_sum_count, to_process)):
                    total_sum += s
                    total_count += c

    mean = total_sum / total_count

    print(f"mean: {mean}")
    counter = 0
    # === SECOND PASS: Compute std ===
    with db.begin() as txn:
        cursor = txn.cursor()
        to_process = []

        with ProcessPoolExecutor(max_workers=num_workers) as executor:

            for _, value in cursor:
                start = time.time()
                to_process.append(value)
                if len(to_process) >= 1000:

                    for s in list(executor.map(extract_squared_diff, to_process, itertools.repeat(mean))):
                        std_sum += s
                    to_process = []
                    print(counter , time.time() - start)
                    counter += 1

            if to_process:
                for s in list(executor.map(extract_squared_diff, to_process, itertools.repeat(mean))):
                    std_sum += s

    std = np.sqrt(std_sum / total_count)

    print(f"std: {std}")

    return mean, std



###### Example to calculate mean and std from an lmdb file: ######
# final_mean, final_std = compute_mean_std("train.lmdb")

###### Comparing different methods: ######
# log_file = r"norm_calculation.txt"
# log = func.config_logger("info", log_file)
# overall_start = time.time()
# final_mean, final_std = compute_mean_std("train.lmdb")
# log.info(f"mean: {final_mean}, std: {final_std}")
# print(time.time() - overall_start)
# print(final_mean)
# print(final_std)
#
# middle = time.time()
# final_mean2, final_std2 = compute_mean_std_2pass("train.lmdb")
#
# print(time.time() - middle)
# print(final_mean2)
# print(final_std2)
# log.info(f"mean2: {final_mean2}, std2: {final_std2}")
