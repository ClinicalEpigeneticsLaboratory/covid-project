from math import log
import itertools

import pandas as pd

class EPIC:
	def __init__(self, epic_path: str):
		self.epic = pd.read_csv(epic_path, index_col=0, low_memory=False)
		self.genes = self.epic["UCSC_RefGene_Name"].dropna()

	def create_bed_file(self, cpgs: list, path: str) -> None:
		data = self.epic.loc[set.intersection(set(cpgs), set(self.epic.index)), ["CHR", "MAPINFO"]]
		data["CHR"]=data["CHR"].map(lambda x: "chr" + str(x))
		data["Start"] = data["MAPINFO"] - 1

		data["End"] = data["MAPINFO"].astype(int)
		data["Start"] = data["Start"].astype(int)

		data.reset_index(inplace=True, drop=False)
		data = data[["CHR", "Start", "End", "IlmnID"]]
		print(data)
		data.to_csv(path, sep="\t", header=None, index=False)
	
	@staticmethod
	def extract_genes(series_of_genes: pd.DataFrame, columnt_to_extract: int = 0) -> pd.DataFrame:
		
		if isinstance(series_of_genes, pd.Series):
			series_of_genes = series_of_genes.to_frame()
		
		df = series_of_genes.iloc[:, columnt_to_extract].str.split(";").tolist()
		df = list(set(itertools.chain(*df)))
		df = pd.DataFrame(df, columns=["Genes"])
    
		return df
