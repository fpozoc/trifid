#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" trifid/models/train_model.py

Training the TRIFID model.

Usage: python -m trifid.model.train \
    --features config/features.yaml \
    --pretrained
__
--help              |-h     Display documentation.
--custom            |-c     Train with a customized model.
--features          |-f     Features selected description yaml filepath.
--model_selection   |-m     Performs the model selection protocol.
--pretrained        |-p     Train TRIFID with a previously trained model.
"""

from __future__ import absolute_import, division, print_function

import argparse, os, pickle, warnings

import numpy as np
import pandas as pd
import random
from sklearn.ensemble import RandomForestClassifier

from .select import Classifier, ModelSelection
from ..utils import utils


def main():
    parser = argparse.ArgumentParser(
        description='Command-line arguments parser', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-c", "--custom",  action='store_true', default=False, 
        help="Training and saving a customized model.")
    parser.add_argument('-f', '--features', default='config/features.yaml', 
        help="Features selected description filepath.")
    parser.add_argument(
        "-m", "--model_selection",  action='store_true', default=False, 
        help="Perform a nested cv model selection, training and saving the best model.")
    parser.add_argument(
        "-p", "--pretrained",  action='store_true', default=False, 
        help="Train TRIFID with a previously trained model.")
    parser.add_argument(
        "-s", "--seed",  default=123,
        help="Perform a nested cv model selection, training and saving the best model.")
    args = parser.parse_args()

    warnings.filterwarnings("ignore")

    df_features = pd.DataFrame(utils.parse_yaml(args.features))
    features = df_features[~df_features['category'].str.contains('Identifier')]['feature'].values

    df_features = pd.read_csv(os.path.join('data', 'genomes', 'GRCh38', 'g27', 'trifid_db.tsv.gz'), sep='\t', compression='gzip')

    df_training_set = pd.read_csv(os.path.join('data', 'model', 'training_set_initial.g27.tsv.gz'), sep='\t')
    df_training_set.loc[df_training_set['state'].str.contains('F'), 'label'] = 1
    df_training_set.loc[df_training_set['state'].str.contains('U'), 'label'] = 0
    df_training_set = df_training_set.loc[~df_training_set['label'].isnull()]
    df_training_set = df_training_set.loc[df_training_set['added'].str.contains('v1|r|v3')].drop(['added', 'state', 'comment'], axis=1).reset_index(drop=True)
    df_training_set.to_csv(os.path.join('data', 'model', 'training_set_final.g27.tsv.gz'), sep='\t', compression='gzip', index=None)

    if args.model_selection:
        ms = ModelSelection(df_training_set, 
            features_col=df_training_set[features],
            target_col='label',
            random_state=args.seed)
        model = ms.get_best_model(outdir='models')

    elif args.custom:
        custom_model = RandomForestClassifier(
            n_estimators=400, 
            class_weight = None, 
            max_features = 7, 
            min_samples_leaf = 7,
            random_state=args.seed)
        model = Classifier(model = custom_model,
            df = df_training_set, 
            features_col=df_training_set[features].columns,
            target_col='label',
            random_state=args.seed)
        model.save_model(outdir='models')

    elif args.pretrained:
        pretrained_model = pickle.load(open(os.path.join('models', 'selected_model.pkl'), 'rb'))
        model = Classifier(model = pretrained_model,
            df = df_training_set, 
            features_col=df_training_set[features].columns,
            target_col='label',
            random_state=args.seed)


if __name__ == "__main__":
    main()