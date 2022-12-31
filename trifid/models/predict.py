#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" trifid/models/predict.py
https://github.com/fpozoc/trifid/blob/master/trifid/models/predict.py

Predicting a data set of isoforms with TRIFID model

Usage: python -m trifid.model.predict  \
    --config config/config.yaml \
    --features config/features.yaml \
    --model models/trifid.v_1_0_0.pkl
__
--help      | -h    Display documentation.
--config    | -c    Configuration filepath.
--features  | -f    Features selected description filepath.
--model     | -m    Load a pretrained model to perform the predictions.
"""

from __future__ import absolute_import, division, print_function

import argparse
import os
import pickle
import warnings

import numpy as np
import pandas as pd

from ..utils.utils import generate_metrics, parse_yaml


def make_predictions(features: list, ids: list, config: dict, model: object, assembly: str, release: str):
    """[summary]

    Args:
        features (list, ids): [description]
    """
    data_dir = os.path.join("data", "genomes", config["annotation"]["genome_version"], config["annotation"]["db"])
    df = pd.read_csv(os.path.join(data_dir, "trifid_db.tsv.gz"), sep="\t", compression="gzip")
    df_ids = df[ids]
    df_features = df[features]

    if assembly == "GRCh38" and release.startswith("g"):
        pass

    elif assembly == "GRCh37" and release.startswith("g"):
        df[df.columns[df.isna().any()].tolist()] = df[df.columns[df.isna().any()].tolist()].fillna(-1)

    elif assembly == "GRCh38" and release.startswith("r"):
        features_no_rs = list(df_features[df_features["refseq"] == "n"]["feature"].values)
        df.loc[df["ccdsid"].str.contains("CCDS"), "CCDS"] = 1
        df.loc[df["ccdsid"].str.contains("-"), "CCDS"] = 0
        df_predictions = df[ids]
        df[features_no_rs] = df[features_no_rs].fillna(-1)
        df[["corsair_alt", "norm_corsair_alt"]] = df[["corsair_alt", "norm_corsair_alt"]].fillna(0)

    elif assembly == "GRCh37" and release.startswith("r"):
        df_features = df[features]
        features_missing = list(df_features[df_features["refseq"] == "n"]["feature"].values) + [
            "RNA2sj_cds",
            "norm_RNA2sj_cds",
        ]
        df.loc[df["ccdsid"].str.contains("CCDS"), "CCDS"] = 1
        df.loc[df["ccdsid"].str.contains("-"), "CCDS"] = 0
        df_predictions = df[ids]
        df[features_missing] = df[features_missing].fillna(-1)
        df[["corsair_alt", "norm_corsair_alt"]] = df[["corsair_alt", "norm_corsair_alt"]].fillna(0)

    elif assembly == "GRCm39":
        features_missing = ["RNA2sj_cds", "norm_RNA2sj_cds"]
        df[features_missing] = df[features_missing].fillna(-1)

    elif assembly == "GRCm38":
        features_missing = ["RNA2sj_cds", "norm_RNA2sj_cds"]
        df[features_missing] = df[features_missing].fillna(-1)

    elif assembly == "Rnor_6.0":  # review ACDS
        df = df.drop("CCDS", axis=1)
        df_acds = pd.read_csv(
            os.path.join("data", "external", "ACDS", "Rnor_6.0", "e104", "acds.tsv.gz"), sep="\t", compression="gzip"
        )
        df_acds = df_acds[["transcript_id", "CCDS"]]
        df = pd.merge(df, df_acds, how="left", on="transcript_id")
        df["CCDS"] = df["CCDS"].fillna(0)

        features_missing = ["RNA2sj_cds", "norm_RNA2sj_cds"]
        df[features_missing] = df[features_missing].fillna(-1)

    elif assembly == "GRCz11":
        df = df.drop("CCDS", axis=1)
        df_acds = pd.read_csv(
            os.path.join("data", "external", "ACDS", "Rnor_6.0", "e104", "acds.tsv.gz"), sep="\t", compression="gzip"
        )
        df_acds = df_acds[["transcript_id", "CCDS"]]
        df = pd.merge(df, df_acds, how="left", on="transcript_id")
        df["CCDS"] = df["CCDS"].fillna(0)

        features_missing = ["RNA2sj_cds", "norm_RNA2sj_cds"]
        df[features_missing] = df[features_missing].fillna(-1)

    elif assembly == "Sscrofa11.1":
        df = df.drop("CCDS", axis=1)
        df_acds = pd.read_csv(
            os.path.join("data", "external", "ACDS", "Rnor_6.0", "e104", "acds.tsv.gz"), sep="\t", compression="gzip"
        )
        df_acds = df_acds[["transcript_id", "CCDS"]]
        df = pd.merge(df, df_acds, how="left", on="transcript_id")
        df["CCDS"] = df["CCDS"].fillna(0)

        features_missing = ["RNA2sj_cds", "norm_RNA2sj_cds"]
        df[features_missing] = df[features_missing].fillna(-1)

    elif assembly == "Pan_tro_3.0":
        df = df.drop("CCDS", axis=1)
        df_acds = pd.read_csv(
            os.path.join("data", "external", "ACDS", "Rnor_6.0", "e104", "acds.tsv.gz"), sep="\t", compression="gzip"
        )
        df_acds = df_acds[["transcript_id", "CCDS"]]
        df = pd.merge(df, df_acds, how="left", on="transcript_id")
        df["CCDS"] = df["CCDS"].fillna(0)

        features_missing = ["RNA2sj_cds", "norm_RNA2sj_cds"]
        df[features_missing] = df[features_missing].fillna(-1)

    elif assembly == "GRCg6a":
        df[df.columns[df.isna().any()].tolist()] = df[df.columns[df.isna().any()].tolist()].fillna(-1)

    elif assembly == "ARS-UCD1.2":
        df[df.columns[df.isna().any()].tolist()] = df[df.columns[df.isna().any()].tolist()].fillna(-1)

    elif assembly == "BDGP6":
        df[df.columns[df.isna().any()].tolist()] = df[df.columns[df.isna().any()].tolist()].fillna(-1)

    elif assembly == "WBcel235":
        df[df.columns[df.isna().any()].tolist()] = df[df.columns[df.isna().any()].tolist()].fillna(-1)
        df = df.fillna(0)

    df_predictions = generate_metrics(df_ids, df_features, model)
    labels = [
        "gene_id",
        "gene_name",
        "transcript_id",
        "translation_id",
        "flags",
        "ccdsid",
        "appris",
        "ann_type",
        "length",
        "trifid_score",
        "norm_trifid_score",
    ]
    df_predictions[labels].to_csv(
        os.path.join(data_dir, "trifid_predictions.tsv.gz"), index=None, sep="\t", compression="gzip"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Command-line arguments parser", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-a", "--assembly", default="GRCh38", help="Genome assembly version.")
    parser.add_argument("-c", "--config", default="config/config.yaml", help="Configuration filepath.")
    parser.add_argument(
        "-f", "--features", default="config/features.yaml", help="Features selected description filepath."
    )
    parser.add_argument(
        "-m", "--model", default="models/trifid.v_1_0_4.pkl", help="Features selected description filepath."
    )
    parser.add_argument("-r", "--release", default="g27", help="Genome release version.")
    args = parser.parse_args()

    warnings.filterwarnings("ignore")

    config = parse_yaml(args.config)
    df_features = pd.DataFrame(parse_yaml(args.features))
    model = pickle.load(open(args.model), "rb")
    features = df_features[~df_features["category"].str.contains("Identifier")]["feature"].values
    ids = df_features[df_features["category"].str.contains("Identifier")]["feature"].values

    make_predictions(features, ids, config, model, assembly=args.assembly, release=args.release)


if __name__ == "__main__":
    main()
