#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" trifid/models/interpret.py
https://github.com/fpozoc/trifid/blob/master/trifid/models/interpret.py

Interpreting and explaining the TRIFID predictions (https://shap.readthedocs.io/en/latest/)

Classes and functions:
    * TreeInterpretation
"""

from __future__ import absolute_import, division, print_function

import os
import re

import numpy as np
import pandas as pd
import shap
from eli5.sklearn import PermutationImportance
from mlxtend.evaluate import feature_importance_permutation
from rfpimp import cv_importances, dropcol_importances, oob_dropcol_importances
from sklearn.feature_selection import mutual_info_classif
from sklearn.metrics import make_scorer, matthews_corrcoef
from sklearn.pipeline import Pipeline

from ..models.select import Splitter, StratifiedKFold
from ..utils.utils import get_id_patterns, merge_dataframes


class TreeInterpretation(Splitter):
    def __init__(
        self,
        model: object,
        df: list,
        features_col: list,
        target_col: str,
        random_state: int = 123,
        test_size: float = 0.25,
        preprocessing: object = None,
    ):
        """Inherits train/test/split functionality from Splitter class.

        Args:
            model (object): sklearn instance model.
            df (list): pandas DataFrame with training set.
            features_col (list): feature names to use as independent variables.
            target_col (str): target name to use as dependent variable.
            random_state (int, optional): Seed. Defaults to 123.
            test_size (float, optional): Size to split the test. Defaults to 0.25.
            preprocessing (object, optional): Preprocessing step to be added to the model pipeline. Defaults to None.
        """
        super().__init__(df, features_col, target_col, random_state, test_size)
        self.model = model
        self.fit()
        self.features_col = features_col

    def fit(self):
        self.model.fit(self.train_features, self.train_target)
        self.predictions = self.model.predict(self.test_features)
        self.probabilities = self.model.predict_proba(self.test_features)[:, 1]

    @property
    def cv_importances(self):
        """
        https://github.com/parrt/random-forest-importances/blob/master/src/rfpimp.py
        """
        df = (
            cv_importances(self.model, self.train_features, self.train_target, k=5)
            .reset_index()
            .rename(columns={"Feature": "feature", "Importance": "cv_feature_importances"})
            .sort_values(by="cv_feature_importances", ascending=False)
        )
        return df

    @property
    def dropcol_importances(self):
        """
        https://explained.ai/rf-importance/
        """
        df = (
            oob_dropcol_importances(self.model, self.train_features, self.train_target)
            .reset_index()
            .rename(columns={"Feature": "feature", "Importance": "dropcol_importances"})
        )
        return df

    @property
    def feature_importances(self):
        """
        https://scikit-learn.org/stable/auto_examples/ensemble/plot_forest_importances.html
        """
        df = (
            pd.DataFrame(self.model.feature_importances_, index=self.train_features.columns)
            .reset_index()
            .rename(columns={"index": "feature", 0: "feature_importances_sklearn"})
            .sort_values(by="feature_importances_sklearn", ascending=False)
        )
        return df

    @property
    def feature_importance_permutation(self):
        """
        http://rasbt.github.io/mlxtend/user_guide/evaluate/feature_importance_permutation/
        """
        imp_vals, imp_all = feature_importance_permutation(
            predict_method=self.model.predict,
            X=self.train_features.values,
            y=self.train_target.values,
            metric="accuracy",
            num_rounds=10,
            seed=self.random_state,
        )

        imp_mean = np.zeros(imp_all.shape[0])
        imp_std = np.zeros(imp_all.shape[0])

        for n, i in enumerate(imp_all):
            imp_mean[n] = imp_all[n].mean()
            imp_std[n] = imp_all[n].std()

        df = (
            pd.DataFrame([imp_mean, imp_std], columns=self.train_features.columns)
            .T.reset_index()
            .rename(columns={"index": "feature", 0: "feature_importance_permutation", 1: "std"})
        )
        df = df[["feature", "feature_importance_permutation"]].sort_values(
            by="feature_importance_permutation", ascending=False
        )
        return df

    @property
    def mutual_information(self):
        """
        Estimate mutual information for a discrete target variable.

        Mutual information (MI) [1] between two random variables is a non-negative
        value, which measures the dependency between the variables. It is equal to
        zero if and only if two random variables are independent, and higher
        values mean higher dependency.

        https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.mutual_info_classif.html
        """
        df = (
            pd.DataFrame(
                mutual_info_classif(self.train_features, self.train_target, random_state=self.random_state),
                index=self.train_features.columns,
            )
            .reset_index()
            .rename(columns={"index": "feature", 0: "mutual_information"})
            .sort_values(by="mutual_information", ascending=False)
        )
        return df

    @property
    def oob_dropcol_importances(self):
        """
        https://explained.ai/rf-importance/
        """
        df = (
            oob_dropcol_importances(self.model, self.train_features, self.train_target)
            .reset_index()
            .rename(columns={"Feature": "feature", "Importance": "oob_dropcol_importances"})
        )
        return df

    @property
    def permutation_importances(self):
        """
        https://eli5.readthedocs.io/en/latest/
        """
        permutation_importance = PermutationImportance(
            self.model,
            random_state=self.random_state,
            scoring=make_scorer(matthews_corrcoef),
            n_iter=10,
            cv=StratifiedKFold(n_splits=10, shuffle=True, random_state=self.random_state),
        ).fit(self.train_features.values, self.train_target.values)
        df = (
            pd.DataFrame(
                [permutation_importance.feature_importances_, permutation_importance.feature_importances_std_],
                columns=self.train_features.columns,
            )
            .T.reset_index()
            .rename(columns={"index": "feature", 0: "permutation_importance", 1: "std"})
        )
        df = df[["feature", "permutation_importance"]].sort_values(by="permutation_importance", ascending=False)
        return df

    @property
    def shap(self):
        """
        https://shap.readthedocs.io/en/latest/
        """
        explainer = shap.TreeExplainer(self.model)
        shap_values = explainer.shap_values(self.train_features)
        vals = np.abs(shap_values).mean(0)
        std_vals = np.abs(shap_values).std(0)
        df = pd.DataFrame(
            list(zip(self.train_features.columns, vals, std_vals)), columns=["feature", "values_mean", "values_std"]
        )
        df["importance"] = df["values_mean"].mean()
        df["std"] = df["values_std"].mean()
        df = df[["feature", "importance", "std"]].reset_index(drop=True).rename(columns={"importance": "shap"})
        df = df[["feature", "shap"]].sort_values(by="shap", ascending=False)
        return df

    @property
    def target_permutation(self):
        """
        Original author Olivier Grellier
        https://www.kaggle.com/ogrellier/feature-selection-target-permutations
        """
        n_splits = 5
        n_runs = 10
        imp_df = np.zeros((len(self.train_features.columns), n_splits * n_runs))
        np.random.seed(self.random_state)
        idx = np.arange(len(self.train_target))
        for run in range(n_runs):
            np.random.shuffle(idx)
            perm_target = self.train_target.iloc[idx]
            folds = StratifiedKFold(n_splits, True, None)
            oof = np.empty(len(self.train_features))
            for fold_, (trn_idx, val_idx) in enumerate(folds.split(perm_target, perm_target)):
                trn_dat, trn_tgt = self.train_features.iloc[trn_idx], perm_target.iloc[trn_idx]
                val_dat, val_tgt = self.train_features.iloc[val_idx], perm_target.iloc[val_idx]
                self.model.fit(trn_dat, trn_tgt)
                imp_df[:, n_splits * run + fold_] = self.model.feature_importances_
                oof[val_idx] = self.model.predict_proba(val_dat)[:, 1]
        bench_imp_df = np.zeros((len(self.train_features.columns), n_splits * n_runs))
        for run in range(n_runs):
            np.random.shuffle(idx)
            perm_target = self.train_target.iloc[idx]
            perm_data = self.train_features.iloc[idx]

            folds = StratifiedKFold(n_splits, True, None)
            oof = np.empty(len(self.train_features))

            for fold_, (trn_idx, val_idx) in enumerate(folds.split(perm_target, perm_target)):
                trn_dat, trn_tgt = perm_data.iloc[trn_idx], perm_target.iloc[trn_idx]
                val_dat, val_tgt = perm_data.iloc[val_idx], perm_target.iloc[val_idx]
                self.model.fit(trn_dat, trn_tgt)
                bench_imp_df[:, n_splits * run + fold_] = self.model.feature_importances_
                oof[val_idx] = self.model.predict_proba(val_dat)[:, 1]

        bench_mean = bench_imp_df.mean(axis=1)
        perm_mean = imp_df.mean(axis=1)

        df = pd.DataFrame(
            list(zip(self.train_features.columns, bench_mean, perm_mean)), columns=["feature", "mdi", "mda_target"]
        )
        df["ratio_mdi-mda"] = df["mdi"] / df["mda_target"]
        df["ratio_mdi-mda"].fillna(0, inplace=True)
        df = df[["feature", "ratio_mdi-mda"]].sort_values(by="ratio_mdi-mda", ascending=False)
        return df

    @property
    def merge_feature_importances(self):
        df = merge_dataframes(
            self.cv_importances,
            self.dropcol_importances,
            self.feature_importance_permutation,
            self.feature_importances,
            self.mutual_information,
            self.oob_dropcol_importances,
            self.permutation_importances,
            self.target_permutation,
            on_type="feature",
        )
        return df

    def local_explanation(self, df_features, sample: str, waterfall: bool = False) -> object:
        """https://github.com/slundberg/shap

        Args:
            df_features ([type]): Features of the sample to interpret.
            sample (str): Ensembl Transcript identifier or Gene Name.
            waterfall (bool, optional): If True, return a waterfall plot. Default is False.

        Returns:
            object: pandas DataFrame with local feature values and shap values.
        """
        if sample.startswith(get_id_patterns()):
            idx = "transcript_id"
        elif re.match("[\s\S]*\.\d[a-z]\.\d", sample):
            idx = "transcript_id"
        else:
            idx = "gene_name"
        df_features = df_features[["gene_name", "transcript_id"] + list(self.features_col)]
        df_idx = df_features.set_index(["gene_name", "transcript_id"])
        df_sample = df_idx.iloc[df_idx.index.get_level_values(idx) == sample]
        explainer = shap.TreeExplainer(self.model)
        if idx == "transcript_id":
            shap_values = explainer.shap_values(df_sample)
            df = pd.DataFrame(
                list(zip(np.abs(shap_values).mean(0)[0], df_sample.values[0])),
                columns=["shap", "feature"],
                index=df_sample.columns,
            ).sort_values("shap", ascending=False)
            if waterfall:
                # shap.plots.waterfall(shap_values[0])
                # shap_value_waterfall = shap_values[0].values[:,0]
                print(df_sample)
                base_value = explainer.expected_value
                shap.plots._waterfall.waterfall_legacy(base_value[0], shap_values[0])

        elif idx == "gene_name":
            explain_gene = {}
            for i in range(0, len(df_sample)):
                explain_gene[df_sample.index[i][1]] = np.abs(explainer.shap_values(df_sample.iloc[i])).mean(0)
            df = pd.DataFrame(explain_gene).T
            df.columns = self.features_col
            df.insert(0, "gene_name", sample)
            df.index.set_names(["transcript_id"])
            df = df.reset_index().rename(columns={"index": "transcript_id"})
            df = df.set_index(["gene_name", "transcript_id"]).T
            df["std"] = df.T.std()
            df["sum"] = df.T.sum()
            df = df.sort_values(by="sum", ascending=False).T.round(3)
        return df.round(3)

    def waterfall_plot(self, model, df_features, sample: str):
        """ """
        if sample.startswith(get_id_patterns()):
            idx = "transcript_id"
        elif re.match("[\s\S]*\.\d[a-z]\.\d", sample):
            idx = "transcript_id"
        else:
            idx = "gene_name"
        df_features = df_features[["gene_name", "transcript_id"] + list(self.features_col)]
        df_idx = df_features.set_index(["gene_name", "transcript_id"])
        df_sample = df_idx.iloc[df_idx.index.get_level_values(idx) == sample]
        explainer = shap.TreeExplainer(model)
