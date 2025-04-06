import pandas as pd
import yaml

# Load the configuration file
with open("config/config.yaml", "r") as config_file:
    config = yaml.safe_load(config_file)

# Access holdout_factors from the config
holdout_factors = config["holdout_factors"]

metadata = pd.read_csv("resources/metadata.csv")

cond1_query = metadata["TrackType"].eq("H3K27ac")
cond1_target = metadata["TrackType"].isin(holdout_factors)
cond2 = metadata["Title"].str.contains("LEIO-|MYO-", na=False)


df_query = metadata.loc[cond1_query & cond2].copy()
df_query["title"] = df_query["Title"].str.split(": ").str[1].str.split(";").str[0]
df_query["sample_id"] = df_query["title"].str.rsplit("-", n=1).str[0]

df_target = metadata.loc[cond1_target & cond2].copy()
df_target["title"] = df_target["Title"].str.split(": ").str[1].str.split(";").str[0]
df_target["sample_id"] = df_target["title"].str.rsplit("-", n=1).str[0]

df = pd.merge(
    df_query.filter(["ID", "TrackType", "CellType", "sample_id"]),
    df_target.filter(["ID", "TrackType", "sample_id"]),
    on="sample_id",
    suffixes=("_query", "_target"),
)
df.to_csv("resources/benchmark_set.csv", index=False)
