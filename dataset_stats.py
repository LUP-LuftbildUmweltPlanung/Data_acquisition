import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import encode_to_lmdb_parquet as lmdb_parquet
import matplotlib.patches as mpatches
import calendar


def invalid_dates(parquet_df):
    invalid_dates = parquet_df[~parquet_df['acquisition'].astype(str).str.match(r'^\d{8}$')]
    print("############ invalid dates1  ##########")
    print(invalid_dates)

def extract_acquisition_column(df):
    df_copy = df.copy()
    df_copy['acquisition_str'] = df_copy['acquisition'].astype(str).str.replace(r"\.0$","",regex=True)
    df_copy = df_copy[df_copy['acquisition_str'].str.match(r'^\d{8}$')]


    # Convert "acquisition" to datetime object
    df_copy['acquisition_date'] = pd.to_datetime(df_copy['acquisition_str'].astype(str), format='%Y%m%d', errors="coerce")

    return df_copy

def merge_subsets(list_paths):
    dfs = []
    for elem in list_paths:
        parquet_df_temp = lmdb_parquet.read_all_from_parquet(elem)

        # Print invalid dates of each dataset:
        #invalid_dates(parquet_df)

        dfs.append(parquet_df_temp)

    parquet_df = pd.concat(dfs, ignore_index=False)
    return parquet_df


def preprocess(df,p='month'):
    if p == 'year':
        df['year'] = df['acquisition_date'].dt.year
    elif p == 'month':
        df['month'] = df['acquisition_date'].dt.month
    counts = df[p].value_counts().sort_index()
    percent = 100 * counts / len(df)
    return percent, counts



def plot_year_distribution(year_percent, year_counts, labels=None, key="month"):
    """ Plots the absolute numbers of acquisition dates grouped by key=year or key=month"""

    group_width = 0.8

    all_years = sorted(set().union(*[s.index for s in year_percent]))

    # unify index
    perc_series = [s.reindex(all_years, fill_value=0.0) for s in year_percent]
    cnt_series = [s.reindex(all_years, fill_value=0.0) for s in year_counts]

    n = len(perc_series)
    x = np.arange(len(all_years))
    bar_w = group_width / max(n, 1)

    fig, ax = plt.subplots(figsize=(10, 5))
    offsets = (np.arange(n) - (n - 1) / 2.0) * bar_w

    max_y = 0.0

    for i, (p, c) in enumerate(zip(cnt_series, perc_series)):
        bars = ax.bar(x + offsets[i], p.values, width=bar_w-0.05, label=labels[i], edgecolor="black")
        for bar, count, year in zip(bars, c.values, all_years):
            if 2014 <= int(year) <= 2019:
                bar.set_hatch("//")

            # write values over bars:
            # h = bar.get_height()
            # if h > 0:
            #     ax.text(bar.get_x() + bar.get_width() / 2,
            #              h + 0.5,
            #              f"{int(count)}",
            #              ha="center", va="bottom", fontsize=18)
        max_y = max(max_y, p.values.max())

    if key == "year":
        historical_patch = mpatches.Patch(facecolor="white", edgecolor="black", hatch="//", label="Historical imagery (2014-2019)")

        handles, labels_ = ax.get_legend_handles_labels()
        handles.append(historical_patch)
        labels_.append("Historical imagery (2014-2019)")


    ax.set_xticks(x)
    ax.set_xticklabels(all_years, fontsize=20)
    if key=="year":
        ax.set_xlabel("Year", fontsize=25)
    elif key=="month":
        ax.set_xlabel("Month", fontsize=25)

    ax.set_ylabel("Number of samples", fontsize=25)
    ax.set_title(f"Distribution of acquisition dates per {key}", fontsize=30)
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    if key == "year":
        ax.legend(handles, labels_, fontsize=30)
    else:
        ax.legend(fontsize=30)
    ax.set_ylim(0, max(5, max_y * 1.15))
    ax.tick_params(labelsize=20)

    #plt.tight_layout()
    plt.show()

def plot_distribution_percent(year_percent, year_counts, labels=None, key="month"):
    """ Plots the percentages of acquisition dates per dataset grouped by key=year or key=month"""

    group_width = 0.8

    all_years = sorted(set().union(*[s.index for s in year_percent]))

    # unify index
    perc_series = [s.reindex(all_years, fill_value=0.0) for s in year_percent]
    cnt_series = [s.reindex(all_years, fill_value=0.0) for s in year_counts]

    n = len(perc_series)
    x = np.arange(len(all_years))
    bar_w = group_width / max(n, 1)

    fig, ax = plt.subplots(figsize=(10, 5))
    offsets = (np.arange(n) - (n - 1) / 2.0) * bar_w

    max_y = 0.0


    for i, (p, c) in enumerate(zip(perc_series, cnt_series)):
        bars = ax.bar(x + offsets[i], p.values, width=bar_w-0.05, label=labels[i], edgecolor="black")  # , color=colors[i])
        for bar, count, year in zip(bars, c.values, all_years):
            if 2014 <= int(year) <= 2019:
                bar.set_hatch("//")

            # write values over bars:
            # h = bar.get_height()
            #         if h > 0:
            #             ax.text(bar.get_x() + bar.get_width() / 2,
            #                      h + 0.5,
            #                      f"{int(count)}",
            #                      ha="center", va="bottom", fontsize=10)
        max_y = max(max_y, p.values.max())

    if key == "year":
        historical_patch = mpatches.Patch(facecolor="white", edgecolor="black", hatch="//", label="Historical imagery (2014-2019)")

        handles, labels_ = ax.get_legend_handles_labels()
        handles.append(historical_patch)
        labels_.append("Historical imagery (2014-2019)")



    ax.set_xticks(x)


    if key=="year":
        ax.set_xlabel("Year", fontsize=30)
        ax.set_xticklabels(all_years, fontsize=30)
    elif key=="month":
        ax.set_xlabel("Month", fontsize=30)
        month_labels = [calendar.month_abbr[m] for m in all_years]
        ax.set_xticks(np.arange(len(all_years)))
        ax.set_xticklabels(month_labels, fontsize=30)

    ax.set_ylabel("Percent (%) per dataset", fontsize=30)
    ax.set_title(f"Distribution of acquisition dates per {key}", fontsize=30)
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    if key == "year":
        ax.legend(handles, labels_, fontsize=30)
    else:
        ax.legend(fontsize=30)
    ax.set_ylim(0, max(5, max_y * 1.15))
    ax.tick_params(labelsize=30)

    #plt.tight_layout()
    plt.show()



######### Example to plot dataset stats: ############
# # PATH to parquet with train metadata
# path_to_parquet_train = [r"train_meta.parquet"]
# # PATH to parquets with vali metadata
# path_to_parquet_vali = [r"vali_spati_temp_ind_meta.parquet",
#                         r"vali_spatially_ind_meta.parquet",
#                         r"vali_temp_ind_meta.parquet"]
# # PATH to parquets with test metadata
# path_to_parquet_test = [r"test_spati_temp_ind_meta.parquet",
#                         r"test_spatially_ind_meta.parquet",
#                         r"test_temp_ind_meta.parquet"]
#
# train_df = extract_acquisition_column(merge_subsets(path_to_parquet_train))
# vali_df = extract_acquisition_column(merge_subsets(path_to_parquet_vali))
# test_df = extract_acquisition_column(merge_subsets(path_to_parquet_test))
#
#
# # plot distribution per year (comment line 206)
# dist_key = "year"
#
# # plot distribution per month
# # dist_key = "month"
#
# labels=["Training", "Validation", "Testing"]
#
# train_percent, train_counts = preprocess(train_df, p=dist_key)
# vali_percent, vali_counts = preprocess(vali_df, p=dist_key)
# test_percent, test_counts = preprocess(test_df, p=dist_key)
#
# plot_distribution_percent([train_percent,vali_percent, test_percent],
#                           [train_counts,vali_counts,test_counts],
#                           labels=labels,
#                           key=dist_key)
#
# plot_year_distribution([train_percent,vali_percent, test_percent],
#                        [train_counts,vali_counts,test_counts],
#                        labels=labels,
#                        key=dist_key)