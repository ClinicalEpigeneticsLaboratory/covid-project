import math
import typing as t
from collections import deque

import numpy as np
import pandas as pd
from numba import jit
from tqdm import tqdm
import scipy.stats as s
from pandas import DataFrame, Series
from statsmodels.stats.multitest import multipletests


class StatsAnalysis:

    def __init__(self, df_target: pd.DataFrame, df_control: pd.DataFrame, epic: pd.DataFrame, max_error: float = 0.05,
                 alpha: float = 0.05):
        self.df_target = df_target
        self.df_control = df_control
        self.max_error = max_error
        self.alpha = alpha
        self.epic = epic[["UCSC_RefGene_Group"]].dropna()

    def extract_probes(self, regions_re: str) -> None:
        """
        Analyse only probes in specific gene regions.
        Use name of region e.g. 'body' or to select more than one use re 
        for example TSS200|TSS1500 to select probes in TSS200 or TSS1500.
		"""
        probes = set(self.epic[self.epic["UCSC_RefGene_Group"].str.contains(regions_re)].index)
        probes = set.intersection(probes, set(self.df_target.columns), set(self.df_control.columns))

        self.df_target = self.df_target[probes]
        self.df_control = self.df_control[probes]

    @staticmethod
    @jit(nopython=True)
    def __calculate_error(sample: np.array) -> float:
        error = np.std(sample) / (math.sqrt(len(sample)))
        return error

    def __calculate_min_sample_size(self, sample: pd.Series) -> int:
        if len(sample) > 30:
            a = s.norm.ppf(q=1 - self.alpha / 2)
        else:
            a = s.t.ppf(q=1 - self.alpha / 2, df=len(sample) - 1)

        sample_size = (a ** 2 * np.var(sample)) / (self.max_error ** 2)

        return math.ceil(sample_size)

    def __perform_test(self, sample_a: Series, sample_b: Series) -> t.Tuple[float, str]:
        _, norm_status_a = s.shapiro(sample_a)
        _, norm_status_b = s.shapiro(sample_b)

        if norm_status_a > self.alpha and norm_status_b > self.alpha:
            _, var_test = s.bartlett(sample_a, sample_b)

            if var_test > 0.05:
                dist = "parameric"
                _, pval = s.f_oneway(sample_a, sample_b)

            else:
                dist = "non-parametric"
                _, pval = s.kruskal(sample_a, sample_b)

        else:
            dist = "non-parametric"
            _, pval = s.kruskal(sample_a, sample_b)

        return pval, dist

    @staticmethod
    def extract(df, threshold: float, alpha: float) -> DataFrame:
        return df[(df["Delta mean"].abs() > threshold) & (df["q-value"] <= alpha)]

    @staticmethod
    @jit(nopython=True)
    def __diff(target_sample: np.array, control_sample: np.array) -> t.Tuple[float, float, float, str]:
        target_mean = np.mean(target_sample)
        control_mean = np.mean(control_sample)
        diff = target_mean - control_mean

        if diff > 0:
            status = "Hypermethylated"
        else:
            status = "Hypomethylated"

        return target_mean, control_mean, diff, status

    def run(self) -> DataFrame:

        results = deque()
        markers = set.intersection(set(self.df_target.columns), set(self.df_control.columns))

        for cpg in tqdm(markers):

            target_sample = self.df_target[cpg].values
            control_sample = self.df_control[cpg].values

            min_control_sample_size = self.__calculate_min_sample_size(control_sample)
            min_target_sample_size = self.__calculate_min_sample_size(target_sample)

            if min_target_sample_size > len(target_sample) or min_control_sample_size > len(control_sample):
                continue

            pval, dist = self.__perform_test(target_sample, control_sample)
            control_mean_error, target_mean_error = self.__calculate_error(control_sample), self.__calculate_error(
                target_sample)
            target_mean, control_mean, delta, status = self.__diff(target_sample, control_sample)

            record = {"CpG": cpg, "p-value": pval,
                      "Control CpG mean": control_mean, "Control CpG error": control_mean_error,
                      "Target CpG mean": target_mean, "Target CpG error": target_mean_error,
                      "Min target sample size": min_target_sample_size,
                      "Min control sample size": min_control_sample_size,
                      "Status": dist, "Delta mean": delta, "Status": status}

            results.append(record)

        results = pd.DataFrame(results).set_index("CpG")
        _, results["q-value"], _, _ = multipletests(results["p-value"], method="fdr_bh")

        return results
