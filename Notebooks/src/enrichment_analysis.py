import pickle
import pathlib
from itertools import chain
import typing as t

import numpy as np
import pandas as pd
import plotly.express as px
from IPython.display import display
from scipy.stats import chisquare
from tqdm import tqdm

path_type = t.Union[str, pathlib.Path]


class EnrichmentAnalysis:
    def __init__(self, report: list, mynorm: pd.DataFrame, manifest_path: path_type,
                 target=["UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island"]):
        self.target = target
        self.report = report
        self.mynorm = set(mynorm.index)
        
        self.manifest = pd.read_csv(manifest_path, low_memory=False, index_col=0)
        self.manifest = self.manifest[(self.manifest["CHR"] != "X") & (self.manifest["CHR"] != "Y")][self.target]
        self.manifest["Relation_to_UCSC_CpG_Island"] = self.manifest["Relation_to_UCSC_CpG_Island"].fillna("OpenSea")
        self.manifest["UCSC_RefGene_Group"] = self.manifest["UCSC_RefGene_Group"].fillna("Unknown")

        self.input_probes = {}
        self.results = {}
        self.bg = {}

    @staticmethod
    def extract(data: str) -> t.List[str]:

        if ";" in data:
            data = data.split(";")
            data = list(set(data))
            return data

        else:
            return [data]

    def prepare_bg(self) -> None:
        bg_probes = self.mynorm

        for target in tqdm(self.target):
            bg = self.manifest.loc[bg_probes, target].map(self.extract).tolist()
            bg = pd.Series((chain(*bg)))
            bg = bg.value_counts(normalize=True)
            bg.name = f"BG_{target}"

            self.bg.update({target: bg})

    def calculate_frequency(self) -> None:

        probes = self.report

        for target in tqdm(self.target):
            input_ = self.manifest.loc[probes, target].map(self.extract).tolist()
            input_ = pd.Series((chain(*input_)))
            input_ = input_.value_counts(normalize=True)
            input_.name = f"Input_{target}"

            self.input_probes.update({target: input_})

    def estimate_fc(self) -> None:

        for target in tqdm(self.target):

            bg = self.bg.get(target)
            input_probes = self.input_probes.get(target)

            data = pd.concat((bg, input_probes), axis=1)
            data["FC"] = data.iloc[:, 1] / data.iloc[:, 0]

            for record in data.itertuples():
                index, bg_freq, target_freq, _ = record

                obs = np.round(np.array([target_freq, 1 - target_freq]) * 100, 0)
                exp = np.array([bg_freq, 1 - bg_freq]) * 100

                _, p = chisquare(f_obs=obs, f_exp=exp)
                data.loc[index, "p-value"] = p

            self.results.update({target: data})

    def vis(self, path = False) -> None:

        for target in tqdm(self.target):

            data = self.results.get(target)
            display(data)

            fig = px.bar(data, x=data.index, y=data.columns[:2])

            fig.update_layout(barmode='group')
            fig.update_xaxes(title="")
            fig.update_yaxes(title="Frequency")

            fig.show()
            
            if path:
                fig.write_image(f"{path}/{target}.jpg")
                data.to_csv(f"{path}/df_{target}.csv")
