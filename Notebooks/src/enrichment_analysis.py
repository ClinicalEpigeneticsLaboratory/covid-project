import os
import pickle
import pathlib
import typing as t

import numpy as np
import pandas as pd
import plotly.express as px
from IPython.display import display
from scipy.stats import chisquare
from tqdm import tqdm

path_type = t.Union[str, pathlib.Path]


class EnrichmentAnalysis:
    def __init__(self, report: list, mynorm: pd.DataFrame, manifest_path: path_type = os.environ.get("POETRY_EPIC"),
                 target: list=["UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island"]):
        
        self.target = target
        self.report = report
        self.mynorm = mynorm.index
        
        self.manifest = pd.read_parquet(manifest_path)
        self.manifest = self.manifest[(self.manifest["CHR"] != "X") & (self.manifest["CHR"] != "Y")][self.target]
        self.manifest["Relation_to_UCSC_CpG_Island"] = self.manifest["Relation_to_UCSC_CpG_Island"].fillna("OpenSea")
        self.manifest["UCSC_RefGene_Group"] = self.manifest["UCSC_RefGene_Group"].fillna("Unknown")

        self.input_probes = {}
        self.results = {}
        self.bg = {}

    def prepare_bg(self) -> None:
        bg_probes = self.mynorm

        for target in tqdm(self.target):
            bg = self.manifest.loc[bg_probes, target].str.split(";").explode()
            bg = bg.value_counts(normalize=True)
            bg.name = f"BG_{target}"
            self.bg.update({target: bg})

    def calculate_frequency(self) -> None:

        probes = self.report

        for target in tqdm(self.target):
            input_ = self.manifest.loc[probes, target].str.split(";").explode()
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
