import pandas as pd

inputs = snakemake.input
path_out = snakemake.output[0]

def load_pixy_data(inputs):
    ds_list = []
    for f in inputs:
        ds_list.append(pd.read_table(f))
    ds = pd.concat(ds_list).assign(pop=lambda x: x["pop"])
    return ds


def calc_pi(ds):
    out = (
        ds.loc[:, ["pop", "count_diffs", "count_comparisons"]]
        .groupby(["pop"])
        .sum()
        .assign(pi=lambda x: x.count_diffs / x.count_comparisons)
        .reset_index()
        .loc[:, ["pop", "pi"]]
    )
    return out


ds_in = load_pixy_data(inputs)
ds_out = calc_pi(ds_in)

ds_out.to_csv(path_out, sep="\t", index=False)
