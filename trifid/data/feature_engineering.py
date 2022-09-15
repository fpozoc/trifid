#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" trifid/data/feature_engineering.py
https://github.com/fpozoc/trifid/blob/master/trifid/data/feature_engineering.py

Feature Engineering Module

It contains functions to let the user to create, deploy and improve the existing
features from our pandas DataFrames.

Classes and functions:
    * build_appris_features
    * build_categorical_features
    * build_corsair_alt_features
    * build_phylocsf_features
    * build_qpfam_featues
    * build_qsplice_features
    * build_features
    * load_data
"""

from __future__ import absolute_import, division, print_function

import os, sys

import numpy as np
import pandas as pd
from loguru import logger 

from . import loaders
from ..utils import utils


def build_appris_features(df:pd.DataFrame)->pd.DataFrame:
    """Prepare the APPRIS DataFrame for the feature engineering.

    Args:
        df (pd.DataFrame): Initial DataFrame with APPRIS features.

    Returns:
        pd.DataFrame: Processed DataFrame with APPRIS features.
    """
    features = ['crash_p', 'crash_m', 'firestar', 'matador3d', 'spade', 'thump']
    df = utils.group_normalization(df, features=features)
    df = utils.group_normalization(df, features=['corsair'], nmax=4) 
    return df


def build_categorical_features(df:pd.DataFrame)->pd.DataFrame:
    """Prepare the one hot encoded categorical features.

    Args:
        df (pd.DataFrame): Initial DataFrame with categorical features.

    Returns:
        pd.DataFrame: Processed DataFrame with categorical features.
    """
    if 'tsl_1' not in df.columns :
        if df['tsl'].shape[0] != 0:
            df = utils.one_hot_encoding(df, ['tsl'])
        else: 
            for i in range(0,7):
                df.insert(len(df.columns), f'tsl_{i}', np.nan)
    if 'level_1' not in df.columns :
        if df['level'].shape[0] != 0:
            df = utils.one_hot_encoding(df, ['level'])
        else: 
            for i in range(0,4):
                df.insert(len(df.columns), f'level_{i}', np.nan)
    else: 
        pass
    return df


def build_corsair_alt_features(df:pd.DataFrame)->pd.DataFrame:
    """Prepare the ALT-corsair DataFrame for the feature engineering.

    Args:
        df (pd.DataFrame): Initial DataFrame with ALT-corsair features.

    Returns:
        pd.DataFrame: Processed DataFrame with ALT-corsair features.
    """
    features = ['corsair_alt']
    if df[features].shape[0] != 0:
        df = utils.group_normalization(df, features=['corsair_alt'], nmax=0.25)
    else:
        for feature in features:
            df.insert(len(df.columns), f'norm_{feature}', np.nan)        
    return df


def build_phylocsf_features(df:pd.DataFrame)->pd.DataFrame:
    """Prepare the PhyloCSF DataFrame for the feature engineering.

    Args:
        df (pd.DataFrame): Initial DataFrame with PhyloCSF features.

    Returns:
        pd.DataFrame: Processed DataFrame with PhyloCSF features.
    """
    features = ['ScorePerCodon', 'RelBranchLength', 'PhyloCSF_Psi']
    if df[features].shape[0] != 0:
        for feature in features:
            df[feature] = df[feature].astype(float)
            df[feature] = df.groupby('gene_id')[feature].transform(lambda x: x.fillna(x.min()))
        df = utils.impute(df, itype='percentile', percentile=3, features=features)
        df = utils.group_normalization(df, features=features)
        df = utils.fragments_correction(df, features=features)
    else:
        for feature in features:
            df.insert(len(df.columns), f'norm_{feature}', np.nan)        
    return df


def build_qpfam_featues(df:pd.DataFrame)->pd.DataFrame:
    """Prepare the QPfam DataFrame for the feature engineering.

    Args:
        df (pd.DataFrame): Initial DataFrame with QPfam features.
    
    Returns:
        pd.DataFrame: Processed DataFrame with QPfam features.
    """
    features = ['Lost_residues_pfam', 'Gain_residues_pfam']
    df = utils.group_normalization(df, features=features)
    df.loc[df['Lost_residues_pfam'] < 10, 'norm_Lost_residues_pfam'] = 0
    df = df.drop(features, axis=1) # keeping only normalized values
    df['spade_loss'] = df.groupby('gene_name')['spade'].transform(lambda x: max(x) - x)
    df.loc[df['spade_loss']>50, 'spade_loss'] = 50
    df = utils.group_normalization(df, features=['spade_loss'])
    df['norm_spade_loss'] = 1 - df['norm_spade_loss']    
    return df


def build_qsplice_features(df:pd.DataFrame)->pd.DataFrame:
    """Prepare the QSplice DataFrame for the feature engineering.

    Args:
        df (pd.DataFrame): Initial DataFrame with QSplice features.

    Returns:
        pd.DataFrame: Processed DataFrame with QSplice features.
    """
    features = ['RNA2sj_cds', 'RNA2sj'] 
    if df[features].shape[0] != 0:
        df['RNA2sj'] = df['RNA2sj'].fillna(0)
        df['RNA2sj_cds'] = df['RNA2sj_cds'].fillna(0)
        df = utils.group_normalization(df, features=features)
        df = utils.fragments_correction(df, features=features)
    else:
        for feature in features:
            df.insert(len(df.columns), f'norm_{feature}', np.nan)
    return df


def build_features(df:pd.DataFrame)->pd.DataFrame:
    """Prepare the DataFrame for the feature engineering.

    Args:
        df (pd.DataFrame): Initial DataFrame with features.
    
    Returns:
        pd.DataFrame: Processed DataFrame with features.
    """
    df = utils.delta_score(df, ['length'])
    df = build_categorical_features(df)
    df = build_appris_features(df)
    df = build_corsair_alt_features(df)
    df = build_phylocsf_features(df)
    df = build_qpfam_featues(df)
    df = build_qsplice_features(df)

    df = utils.unity_ranger(df, features=[col for col in df if col.startswith('norm')])
    df = utils.reorder_cols(df)
    return df


def load_data(config:dict, assembly:str, release:str)->pd.DataFrame:
    """Load the data from the databases.

    Args:
        config (dict): Configuration dictionary.
        assembly (str): Assembly name.
        release (str): Release name.

    Returns:
        pd.DataFrame: DataFrame with the data.
    """

    annotation_path = config['genomes'][assembly][release]['annotation']
    appris_path = config['genomes'][assembly][release]['appris_data']
    corsair_alt_path = config['genomes'][assembly][release]['corsair_alt']
    qpfam_path = config['genomes'][assembly][release]['qpfam']
    qsplice_path = config['genomes'][assembly][release]['qsplice']
    phylocsf_path = config['genomes'][assembly][release]['phylocsf']
    reference_path = config['genomes'][assembly][release]['reference']
    sequences_path = config['genomes'][assembly][release]['sequences']

    try:
        if 'gtf' in annotation_path:
            df_annotation = loaders.load_annotation(
                annotation_path, db=release[0])
        elif '-' in annotation_path:
            df_annotation = pd.DataFrame(
                columns=['transcript_id', 'basic', 'CCDS', 'level_1', 'level_2', 'level_3',
                        'nonsense_mediated_decay', 'readthrough', 'StartEnd_NF', 'tsl_1',
                        'tsl_2', 'tsl_3', 'tsl_4', 'tsl_5', 'tsl_6'])
        logger.info(utils.get_df_info(df_annotation))
    except NameError:
        print('Please, provide a right genome annotation filepath to serve as reference in the config file.') 

    try:
        if 'appris_data' in appris_path:
            df_appris = loaders.load_appris(appris_path)
        logger.info(utils.get_df_info(df_appris))
    except NameError:
        print('Please, provide a right APPRIS filepath to serve as reference in the config file.') 

    try:
        if 'corsair_alt' in corsair_alt_path:
            df_corsair_alt = loaders.load_corsair_alt(corsair_alt_path)
        elif '-' in corsair_alt_path:
            df_corsair_alt = pd.DataFrame(
                columns=['transcript_id', 'corsair_alt'])
        logger.info(utils.get_df_info(df_corsair_alt))
    except NameError:
        print('Please, provide a right ALT-corsair filepath in the config file.') 

    try:
        if 'qpfam' in qpfam_path:
            df_qpfam = loaders.load_qpfam(qpfam_path)
        elif '-' in qpfam_path:
            df_qpfam = pd.DataFrame(
                columns=['transcript_id', 'pfam_score', 'pfam_domains_impact_score',
                'perc_Damaged_State', 'perc_Lost_State', 'Lost_residues_pfam',
                 'Gain_residues_pfam'])
        logger.info(utils.get_df_info(df_qpfam))
    except NameError:
        print('Please, provide a right QPfam filepath in the config file. Visit https://gitlab.com/fpozoc/qpfam for more info.') 

    try:
        if 'qsplice' in qsplice_path:
            df_qsplice = loaders.load_qsplice(qsplice_path)
        elif '-' in qsplice_path:
            df_qsplice = pd.DataFrame(
                columns=['transcript_id', 'RNA2sj', 'RNA2sj_cds'])
        logger.info(utils.get_df_info(df_qsplice))
    except:
        print('Please, provide a right QSplice filepath in the config file. Visit https://gitlab.com/fpozoc/qsplice for more info.') 
    
    try:
        if 'PhyloCSF' in phylocsf_path:
            df_phylocsf = loaders.load_phylocsf(phylocsf_path)
        elif '-' in phylocsf_path:
            df_phylocsf = pd.DataFrame(
                columns=['transcript_id', 'ScorePerCodon', 'RelBranchLength', 'PhyloCSF_Psi'])
        logger.info(utils.get_df_info(df_phylocsf))
    except:
        print('Please, provide a right PhyloCSF filepath in the config file.') 

    try:
        if 'qduplications' in reference_path:
            df_reference = loaders.load_reference(reference_path)
        elif '-' in reference_path:
            df_reference = pd.DataFrame(
                columns=['transcript_id', 'ann_type', 'transcript_ref'])
        logger.info(utils.get_df_info(df_reference))
    except:
        print('Please, provide a right annotations reference filepath in the config file.')

    try:
        if 'fa' in sequences_path:
            df_sequences = loaders.load_sequences(sequences_path)
        elif '-' in sequences_path:
            df_sequences = pd.DataFrame(
                columns=['transcript_id', 'sequence'])
        logger.info(utils.get_df_info(df_sequences))
    except NameError:
        print('Please, provide a right FASTA sequences filepath in the config file.') 

    df = utils.merge_dataframes(df_appris, df_sequences, df_annotation, 
                               df_corsair_alt, df_reference, df_qsplice, df_qpfam, df_phylocsf)

    df = df.loc[~df['transcript_id'].str.contains('PAR')] # removing pars 
    return df