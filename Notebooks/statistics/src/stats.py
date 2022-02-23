import math
import typing as t
from collections import deque

import numpy as np
import pandas as pd
import scipy.stats as s
from numba import jit
from pandas import DataFrame, Series
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm


class StatsAnalysis:

    def __init__(self, df_target: pd.DataFrame, df_control: pd.DataFrame, epic: pd.DataFrame, alpha: float = 0.05):
        self.df_target = df_target
        self.df_control = df_control
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
        probes = list(probes)

        self.df_target = self.df_target.loc[:, probes]
        self.df_control = self.df_control.loc[:, probes]

    @staticmethod
    @jit(nopython=True)
    def __calculate_error(sample: np.array) -> float:
        error = np.std(sample) / (math.sqrt(len(sample)))
        return error

    def __perform_test(self, sample_a: Series, sample_b: Series) -> t.Tuple[float, str]:
        _, norm_status_a = s.shapiro(sample_a)
        _, norm_status_b = s.shapiro(sample_b)

        if norm_status_a > self.alpha and norm_status_b > self.alpha:
            _, var_test = s.bartlett(sample_a, sample_b)

            if var_test > self.alpha:
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
        return df[(df["Delta mean"].abs() > threshold) & (df["Adj. p-value"] <= alpha)]

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

    def run(self, name: str = "Run") -> DataFrame:

        results = deque()
        markers = set.intersection(set(self.df_target.columns), set(self.df_control.columns))

        for cpg in tqdm(markers, desc=name):

            target_sample = self.df_target[cpg].values
            control_sample = self.df_control[cpg].values

            pval, dist = self.__perform_test(target_sample, control_sample)
            control_mean_error, target_mean_error = self.__calculate_error(control_sample), self.__calculate_error(target_sample)
            target_mean, control_mean, delta, status = self.__diff(target_sample, control_sample)

            record = {"CpG": cpg, "p-value": pval,
                      "Control CpG mean": control_mean, "Control CpG error": control_mean_error,
                      "Target CpG mean": target_mean, "Target CpG error": target_mean_error,
                      "Status": dist, "Delta mean": delta, "Status": status}

            results.append(record)

        results = pd.DataFrame(results).set_index("CpG")
        _, results["Adj. p-value"], _, _ = multipletests(results["p-value"], method="fdr_bh")

        return results
