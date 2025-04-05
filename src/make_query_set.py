import pandas as pd

df = pd.read_csv("resources/metadata.csv")

cond1 = df["TrackType"].eq("H3K27ac")
cond2 = df["Title"].str.contains("LEIO-|MYO-", na=False)
query_set = df.loc[cond1 & cond2, "ID"].tolist()

# Save the query set to a txt file
with open("resources/query_set.txt", "w") as f:
    for item in query_set:
        f.write(f"{item}\n")
