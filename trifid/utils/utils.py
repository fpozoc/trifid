#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" trifid/utils/utils.py
https://gitlab.com/bu_cnio/trifid/-/blob/master/trifid/models/utils.py

Helper Functions

It contains some convenient utilities not included with the standard modules

Classes and functions:
    * Statistics
    * balanced_training_set: returns a 1:1 proportion labeled training set.
    * create_dir: returns the path to directory recently created.
    * delta_score - returns a pandas DataFrame with delta Score feature (commonly used for length).
    * fragments_correction - returns a pandas DataFrame with fragments isoforms corrected by feature.
    * generate_training_set - returns a pandas DataFrame with training set for TRIFID ready.
    * generate_trifid_metrics - returns a pandas DataFrame with TRIFID scores calculated.
    * get_df_info: returns pandas pandas DataFrame useful info.
    * get_id_patterns - returns a tuple with a set of transcript identifiers initial patterns.
    * group_normalization - returns a pandas DataFrame with features normalized. 
    * impute - returns a pandas DataFrame with features imputed.
    * merge_dataframes: returns a merged DataFrame on transcript identifier.
    * one_hot_encoding - returns a pandas DataFrame with features encoded.
    * reduce_mem_usage: returns a pandas pandas DataFrame with optimized memory usage.
    * reorder_cols: returns a pandas DataFrame with columns reordered by object and then float columns.
    * round_df_floats: returns pandas DataFrame with float features rounded.
    * parse_yaml: returns a pandas DataFrame with new features.
    * timer: returns info from time consumed in a coding block.
    * unity_ranger - returns a pandas DataFrame with features between 0 and 1. 
"""

from __future__ import absolute_import, division, print_function

from contextlib import contextmanager
import functools, gzip, os, yaml

import numpy as np
import pandas as pd
from loguru import logger 


class Statistics():
    def __init__(self, df:list, nr:bool=False):
        if nr:
            self.df = self._remove_redundants(df)
        else:
            self.df = df
        self.df_pri = self.df[self.df['appris'].str.contains('PRINCIPAL')].reset_index(drop=True)
        self.df_alt = self.df[~self.df['appris'].str.contains('PRINCIPAL')].reset_index(drop=True)
        
    def get_stats(self, cutoff:float=0.5, norm_double_check:bool=False, cutoff_feature:str='trifid_score'):

        if norm_double_check:
            pf = ((self.df_pri['trifid_score']>=cutoff) & (self.df_pri['norm_trifid_score_max']>=cutoff)).sum()
            pnf = ((self.df_pri['trifid_score']<cutoff) & (self.df_pri['norm_trifid_score_max']<cutoff)).sum()
            af = ((self.df_alt['trifid_score']>=cutoff) & (self.df_alt['norm_trifid_score_max']>=cutoff)).sum()
            anf = ((self.df_alt['trifid_score']<cutoff) & (self.df_alt['norm_trifid_score_max']<cutoff)).sum()
        else:
            pf = (self.df_pri[cutoff_feature]>=cutoff).sum()
            pnf = (self.df_pri[cutoff_feature]<cutoff).sum() 
            af = (self.df_alt[cutoff_feature]>=cutoff).sum() 
            anf = (self.df_alt[cutoff_feature]<cutoff).sum() 

        p = self.df_pri.shape[0]
        a = self.df_alt.shape[0]

        stats = {
            "PRINCIPAL": {
                "Functional": pf,
                "Non functional": pnf,
                "Total": p,
                "Percentage of functional": round(((pf)/p)*100, 2),
            },
            "ALTERNATIVE": {
                "Functional": af,
                "Non functional": anf,
                "Total": a,
                "Percentage of functional": round(((af)/a)*100, 2),
            },
            "Total": {
                "Functional": pf+af,
                "Non functional": pnf+anf,
                "Total": p+a,
                "Percentage of functional": round(((pf+af)/(p+a))*100,2),
            }
        }
        df_stats = pd.DataFrame(stats).T
        df_stats[df_stats.columns[0:3]] = df_stats[df_stats.columns[0:3]].astype(int)
        return df_stats
    
    def _remove_redundants(self, df):
        df = df.loc[df['nr'] == 'yes'
        ].loc[df['flags'].str.contains('^protein_coding$')].reset_index(drop=True)
        return df


def balanced_training_set(df:pd.DataFrame, seed:int=1) -> pd.DataFrame:
    """It creates a balanced proportion of the labels.

    Args:
        df (pd.DataFrame): labelled data set.
        seed (int, optional): Random state. Defaults to 1.

    Returns:
        pd.DataFrame: Balanced data set.
    """    
    return pd.concat([
        df[df['label'] == 1], 
        df[df['label'] == 0].sample(df[df['label'] == 1].shape[0], random_state=seed)
    ]).reset_index(drop=True)


def create_dir(dirpath: str) -> str: 
    """mkdir -p python equivalent
    
    Arguments:
        dirpath {str} -- Path to create the new folder
    
    Returns:
        absolute_path {str} -- Absolute path of the new folder
    """    
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    absolute_path = os.path.abspath(dirpath)


def delta_score(df:pd.DataFrame, features:list, mode:str='appris', groupby_feature:str='gene_id') -> pd.DataFrame:
    """Lenght Delta Function. (It also apply for any other feature).

    It creates a feature that substract the largest feature value of isoform to the current one. 
    Afterwards, it creates an score of normalizing by the largest length of isoform per gene.

    Args:
        df (pd.DataFrame): input data set.
        features (list): List (variable length) of features to transform.
        mode (str, compulsory): Mode to select the reference. APPRIS PRINCIPAL or longest isoform.
        groupby_feature (str, optional): Feature to calculate the score differences. Defaults to 'gene_id'.

    Returns:
        pd.DataFrame: data set that returns the input with normalization column updated.
    """    
    for feature in features:
        df['ccdsid_number'] = df['ccdsid'].str.split('S').str[1].astype(float).fillna(0)
        df.loc[df['tsl'] == 1, 'tsl_first'] = 1
        df.loc[~df['tsl'] == 1, 'tsl_first'] = 0
        df.loc[df['tsl'] == 2, 'tsl_second'] = 1
        df.loc[~df['tsl'] == 2, 'tsl_second'] = 0
        df = df.sort_values(
            by=['gene_id', 'score', 'corsair_alt', 'tsl_first', 'tsl_second', 'ccdsid_number', 'length'],
            ascending=[False, False, False, False, False, True, False])
        if mode == 'longest':
            df[f'{feature}_reference'] = df.groupby(
                ['gene_id'])[feature].transform(max)
            df[f'{feature}_delta'] = df[f'{feature}_reference']-df[feature]
            df[f'{feature}_delta_score'] = (df[f'{feature}_reference']-df[f'{feature}_delta'])/df[f'{feature}_reference']
        elif mode == 'appris':
            df[f'{feature}_reference'] = df.groupby(
                ['gene_id'])[feature].transform('first')
            df[f'{feature}_delta_score'] = 1-abs((df[f'{feature}_reference']-df[feature])/df[f'{feature}_reference'])
            df[f'{feature}_delta_score'][df[f'{feature}_delta_score'] > 1] = 1
            df[f'{feature}_delta_score'][df[f'{feature}_delta_score'] < 0] = 0
        df = df.drop([f'{feature}_reference', 'ccdsid_number', 'tsl_second'], 1)
        logger.info(f'{feature}_delta created.')
        logger.info(df[f'{feature}_delta_score'].describe().round(1).to_dict())
    return df


def fragments_correction(df_iso:pd.DataFrame, features:list) -> pd.DataFrame:
    """Fragments corrector Function

    Correcting the value of sequence fragments by the homologous sequence.

    Args:
        df_iso (pd.DataFrame): data set isoforms database.
        features (list): List of features to correct

    Returns:
        pd.DataFrame: data set with feature columns of isoforms fragments corrected by fragments.
    """    
    df_iso['ann_type'] = df_iso['ann_type'].fillna('-')
    df_iso = df_iso.reset_index(drop=True)  

    appris_labels = ['Principal', 'Alternative']
    df_labels = []
    for label in appris_labels:
        if label == 'Principal':
            df = df_iso.loc[df_iso['ann_type'].str.contains(
                f'^Redundant {label}$|^{label}$')]
        elif label == 'Alternative':
            df = df_iso.loc[
                (df_iso['ann_type'] == 'Redundant {}'.format(label)) | 
                (df_iso['transcript_id'].str.contains(
                    '|'.join(list(df_iso['transcript_ref'].unique()[1:]))))]
            df = df.drop('transcript_ref', axis=1)
        df = df[['gene_id', 'transcript_id',  'ann_type', 'StartEnd_NF']+features]
        for feature in features:
            df[f'check_{feature}'] = df.loc[df['ann_type'] == label][feature]
            df[f'check_{feature}'] = df.groupby('gene_id')[f'check_{feature}'].transform('max')
            df.loc[df[feature] < df[f'check_{feature}'], f'check_{feature}'] = df[feature]
            df[feature] = df[f'check_{feature}']
            df = df.drop(f'check_{feature}', axis=1)
        df_labels.append(df)
    df = pd.concat([df_labels[0], df_labels[1]], axis = 0)

    df = pd.merge(df_iso[['transcript_id'] + features], df, on = 'transcript_id', how='left')
    features_x = df[[col for col in df if col.endswith('x')]].columns
    features_y = df[[col for col in df if col.endswith('y')]].columns
    for feature_x, feature_y in zip(features_x, features_y):
        df.loc[df[feature_y].isnull(), feature_y] = df[feature_x]
    for n, _ in enumerate(features_y):
        df = df.rename(columns={features_y[n]: features[n]})
    df = df.drop(features_x, axis=1)
    df = df[['transcript_id'] + features]
    logger.info(f'Features corrected:\n---')
    for featfix in features:
        logger.info(f'{featfix}: {len(np.where(df[featfix] != df_iso[featfix])[0])}')
    for feature in features:
        df_iso[feature] = df[feature]
    return df_iso 


def generate_training_set(df_features:pd.DataFrame, filepath:str) -> pd.DataFrame:
    """Creating a training set from GENCODE 27 isoforms and their proteomics evidence.

    Args:
        df_features (pd.DataFrame): TRIFID database data set from GENCODE27.
        filepath (str): Filepath to load the labeled isoforms.

    Returns:
        pd.DataFrame: data set with the training set.
    """    
    df = pd.read_csv(filepath, sep='\t', compression='gzip')
    df.loc[df['evidence'].str.contains('Functional|Principal'), 'label'] = 1
    df.loc[~df['evidence'].str.contains('Functional|Principal|Neutral'), 'label'] = 0
    df = df[['transcript_id', 'label']]
    df = pd.merge(df_features, df, how='left', on='transcript_id').drop_duplicates('transcript_id').reset_index(drop=True)
    df = df.loc[~df['label'].isnull()].reset_index(drop=True)
    df['duplicate_score'] = (
        df['tsl_1'] *1.5 + 
        df['tsl_2'] * 1.25 + 
        df['tsl_3'] * 1 + 
        df['tsl_4'] * 0.75 + 
        df['tsl_5'] * 0.5 + 
        df['level_1'] * 1.25 + 
        df['level_2'] * 0.75 + 
        df['level_3'] * -1 + 
        df['norm_ScorePerCodon'] * 1.5 + 
        df['norm_PhyloCSF_Psi'] *1.5 + 
        df['norm_RelBranchLength']*1.5 + 
        df['nonsense_mediated_decay']*0.01
    )
    df = df.sort_values(by=['gene_name', 'duplicate_score', 'sequence'], ascending=[True, False, True])
    df = df.drop('duplicate_score', axis=1).reset_index(drop=True)
    df = df.drop_duplicates(subset=['sequence'], keep='first').reset_index(drop=True)
    df['label'] = df['label'].astype(int)
    return df


def generate_trifid_metrics(df:pd.DataFrame, features:pd.DataFrame, model:object, nmax_norm_median:bool=False) -> pd.DataFrame:
    """Generate TRIFID scores to evaluate the predictions of the method.

    Args:
        df (pd.DataFrame): data set to updaate the predictions.
        features (pd.DataFrame): data set with features to make the predictions.
        model (object): scikit-learn model to make the predictions.
        nmax_norm_median (float): number of min maximum thresohld in  the normalization procedure.

    Returns:
        pd.DataFrame: data set with metrics updated.
    """    
    if hasattr(model, 'make_prediction'):
        df['trifid_score'] = model.make_prediction(features, probability=True)
    else:
        df['trifid_score'] = model.predict_proba(features)[:, 1]
    if nmax_norm_median:
        nmax_norm_median = df['trifid_score'].median()
    else:
        nmax_norm_median = 0.5
    df['norm_trifid_score'] = df.groupby('gene_id')['trifid_score'].transform(
                lambda x: 0 if (x == 0).all() else (
                    1 if ((len(set(x)) == 1) & (x >= 0.5).all()) else (
                        x) / (max(nmax_norm_median, x.max()))))
    df = df.round({
        'trifid_score': 4, 
        'norm_trifid_score': 4
        })
    df.insert(8, 'length', features['length'])
    return df


def get_df_info(df:pd.DataFrame) -> str:
    """Getting pandas DataFrame info for logging.

    Args:
        df (pd.DataFrame): data set as input.

    Returns:
        str: Info message.
    """    
    ncols = df.shape[1]
    nrows = df.shape[0]
    colstr = ';'.join(df.columns)
    msg = f"{ncols} features ({colstr}) with {nrows} instances loaded."
    return msg


def get_id_patterns() -> tuple:
    """Returns a tuple with a set of transcript identifiers initial patterns.

    Args:

    Returns:
        tuple: Regex patterns for species selected.
    """    
    tids_reg = {
            "ensembl":
            {
                "homo_sapiens": "ENST0", 
                "mouse": "ENSMUST",
                "danio_rerio": "ENSDART",
                "rattus_norvegicus": "ENSRNOT",
                "sus_scrofa": "ENSSSCT",
                "pan_troglodytes": "ENSPTRT",
                "gallus_gallus": "ENSGALT",
                "bos_taurus": "ENSBTAT",
                "drosophila_melangaster": "FBtr",
                # "caenorhabditis_elegans": "[\s\S]*\.\d[a-z]\.\d",
            },
            "refseq":
            {
                "homo_sapiens": ("NM", "XM", "YP")
            }
        }
    return tuple(tids_reg['ensembl'].values())+tuple(tids_reg['refseq'].values())[0]

def get_nr_transcripts(df:pd.DataFrame) -> pd.DataFrame:
    """Getting non redundant transcripts from genome.

    Args:
        df (pd.DataFrame): data set as input.

    Returns:
        pd.DataFrame: Non redundant transcript lists.
    """
    df = df.loc[df['ann_type'].str.contains('^Alternative$|^Alternative Duplication$|^Principal$|^Principal Duplication')]
    nr_transcripts = df.groupby(['gene_id', 'sequence'])['trifid_score', 'transcript_id'].max().reset_index().transcript_id.values
    return nr_transcripts


def group_normalization(df:pd.DataFrame, features:list, nmax:int=0, nmin:int=0, groupby_feature:str='gene_id')->pd.DataFrame:
    """Normalizer by group.

    Args:
        df (pd.DataFrame): input data set.
        features (list): list of features to normalize.
        nmax (int, optional): maximum number of score value. Defaults to 0.
        nmin (int, optional): minimum number of score value. Defaults to 0.
        groupby_feature (str, optional): feature to group. Defaults to 'gene_id'.

    Returns:
        pd.DataFrame: data set with normalized features added.
    """    

    phigh = 0.9 # Percentile 90
    plow = 0.1 # Percentile 10

    for feature in features:
        if feature in ['firestar', 'matador3d', 'spade']:
            df[f'norm_{feature}'] = df.groupby(groupby_feature)[feature].transform(
                lambda x: 0 if (x == 0).all() else (
                    1 if (
                        len(set(x)) == 1) else (
                            x-min(nmin, x.min())) / (max(nmax, x.max())-min(nmin, x.min())))
                )
        else:
            df[f'{feature}'] = df[f'{feature}'].astype(float)
            df[f'norm_{feature}'] = df.groupby(groupby_feature)[feature].transform(
                lambda x: 0 if (x == 0).all() else (
                    (x-df[feature].quantile(plow))/(df[feature].quantile(phigh)- df[feature].quantile(plow)) if (
                        len(set(x)) == 1) else (
                            x-min(nmin, x.min())) / (max(nmax, x.max())-min(nmin, x.min())))
                )
        logger.info(f'{feature} normalized.')
    return df


def impute(df:pd.DataFrame, features:list, n:int=None, itype:str='class', column:str=None, condition:str=None, percentile:float=None) -> pd.DataFrame:
    """Imputator

    It contains some functions to impute vectors in different ways.

    Args:
        df (pd.DataFrame): input data set.
        features (list): List of features to normalize.
        n (int, optional): Maximum number of score value. Defaults to None.
        itype (str, optional): Type of imputation to perform. Defaults to 'class'.
        column (str, optional): Column to select in condition. Defaults to None.
        condition(str, optional): Value to select in condition. Defaults to None.
        percentile(float, optional): Percentile to impute with. Defaults to None.

    Returns:
        pd.DataFrame: data set with imputed columns.
    """    
    for feature in features:
        if itype == 'class':
            df.loc[df[feature].isnull(), feature] = n
        elif itype == 'conditional':
            df.loc[df[column].str.contains(condition, na=False), feature] = n
        elif itype == 'percentile':
            df.loc[df[feature].isnull(), feature] = df[feature].quantile(percentile/100)  
        elif itype == 'same_as_norm':
            df.loc[df[feature].isnull(), feature] = df[feature.replace('norm_','')]         
    return df


def merge_dataframes(*args, on_type:str='transcript_id', how_type:str='left', pivot_on:int=0, nimpute:int=None) -> pd.DataFrame:
    """Pandas DataFrame merger Function.
    
    This function merges several DataFrame to create an unique database which contains same isoforms 
    as reference database. It uses pandas merge method.
    
    https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.merge.html


    Args:
        on_type (str, optional): Merge on feature selected. Defaults to 
    'transcript_id'.
        how_type (str, optional): Merge method. Defaults to 'left'.
        pivot_on (int, optional): Represent the first DataFrame to merge when 
    how_type='left'. Defaults to 0.
        nimpute (int, optional): If user wants to impute some feature, nimpute
    is the number to fill na. Defaults to None.

    Raises:
        ValueError: Merge function has to receive more than one DataFrame inside 
    list argument

    Returns:
        pd.DataFrame: pandas DataFrame mergefd
    """    
    args = list(args)
    if len(args) <= 1:
        raise ValueError(
            'Merge function has to receive more than one DataFrame inside list argument.')
    if how_type != 'left':
        pivot_on == None
    args.insert(0, args.pop(pivot_on))
    df = functools.reduce(lambda left, right: pd.merge(
        left, right, on=on_type, how=how_type),args).drop_duplicates(subset=on_type).reset_index(drop=True)
    if nimpute != None:
        df = df.fillna(nimpute)
    return df


def one_hot_encoding(df:pd.DataFrame, features:list) -> pd.DataFrame:
    """One Hot Encoder

    It encodes features selected as "one hot" mood and removing initial feature.
    https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.OneHotEncoder.html

    Args:
        df (pd.DataFrame): input data set.
        features (list): Feature list to perform the method

    Returns:
        pd.DataFrame: data set with encoded features.
    """    
    for feature in features:
        df[feature] = df[feature].astype(int)
        one_hot = pd.get_dummies(df[feature], prefix=feature)
        df = df.drop(feature, axis=1)
        df = df.join(one_hot)
    return df


def open_files(filepath:str) -> object:
    """Openning both compressed and non-compressed files.

    Args:
        filepath (str): File path to the file.

    Returns:
        str: open file object.
    """    
    if filepath.endswith('.gz'):
        open_file = gzip.open(filepath, 'rt')
    else: 
        open_file = open(filepath, 'r')
    return open_file


def parse_yaml(yaml_file:str) -> str:
    '''YAML parsing function

    This function parses a configuration file in yaml format (http://zetcode.com/python/yaml/).

    Parameters
    ----------
    yaml_file: str
        Config file path in yaml format.

    Returns
    -------
    config: dict
        Dictionary with configuration data structure.
    '''
    with open(yaml_file, 'r') as config:
        try:
            config = yaml.safe_load(config)
        except yaml.YAMLError as exc:
            logger.info(exc)
    return config


def reduce_mem_usage(df:pd.DataFrame, verbose:bool=False, round_float:int=False) -> pd.DataFrame:
    """Memory reducer Function

    It reduces memory usage of pandas DataFrame. 
    Inspired in https://www.kaggle.com/artgor/artgor-utils

    Args:
        df (pd.DataFrame): input data set.
        verbose (bool, optional): Verbosity control. Defaults to False.
        round_float (int, optional): Round to 4 dec. Defaults to False.

    Returns:
        pd.DataFrame: data set with reduced memory usage.
    """
    start_mem_usg = df.memory_usage().sum() / 1024**2
    logger.info("Memory usage of properties dataframe is :", start_mem_usg, " MB")
    NAlist = []
    for col in df.columns:
        if df[col].dtype != object:
            if verbose:
                logger.info(
                    "******************************\nColumn: {}\ndtype before: {}\n".format(col, df[col].dtype))
            IsInt = False
            mx = df[col].max()
            mn = df[col].min()
            if not np.isfinite(df[col]).all():
                NAlist.append(col)
                df[col].fillna(mn-1, inplace=True)
            asint = df[col].fillna(0).astype(np.int64)
            result = (df[col] - asint)
            result = result.sum()
            if result > -0.01 and result < 0.01:
                IsInt = True
            if IsInt:
                if mn >= 0:
                    if mx < 255:
                        df[col] = df[col].astype(np.uint8)
                    elif mx < 65535:
                        df[col] = df[col].astype(np.uint16)
                    elif mx < 4294967295:
                        df[col] = df[col].astype(np.uint32)
                    else:
                        df[col] = df[col].astype(np.uint64)
                else:
                    if mn > np.iinfo(np.int8).min and mx < np.iinfo(np.int8).max:
                        df[col] = df[col].astype(np.int8)
                    elif mn > np.iinfo(np.int16).min and mx < np.iinfo(np.int16).max:
                        df[col] = df[col].astype(np.int16)
                    elif mn > np.iinfo(np.int32).min and mx < np.iinfo(np.int32).max:
                        df[col] = df[col].astype(np.int32)
                    elif mn > np.iinfo(np.int64).min and mx < np.iinfo(np.int64).max:
                        df[col] = df[col].astype(np.int64)
            else:
                df[col] = df[col].astype(np.float32)
                if round_float:
                    df[col] = df[col].round(round_float)
            if verbose:
                logger.info("dtype after: ", df[col].dtype)
    logger.info("___MEMORY USAGE AFTER COMPLETION:___")
    mem_usg = df.memory_usage().sum() / 1024**2
    logger.info("Memory usage is: ", mem_usg, " MB")
    logger.info("This is ", 100*mem_usg/start_mem_usg, "% of the initial size")
    return df, NAlist


def reorder_cols(df:pd.DataFrame) -> pd.DataFrame:
    """Reorder columns Function

    Ordering columns of DataFrame: Strings at the start of the df

    Args:
        df (pd.DataFrame): pandas DataFrame

    Returns:
        pd.DataFrame: pandas DataFrame
    """    
    return df[list(df.select_dtypes(include='object').columns)+list(df.select_dtypes(exclude='object').columns)]


def round_df_floats(df:pd.DataFrame, n:int=4) -> pd.DataFrame:
    """Rounder floats Function

    It rounds to 4 (default) all the floated columns of the DataFrame

    Args:
        df (pd.DataFrame): input data set.
        n (int, optional): round integer. Defaults to 4.

    Returns:
        pd.DataFrame: data set with floated columns rounded
    """    
    df[df.select_dtypes(include='float', exclude=None).columns] = df.select_dtypes(include='float', exclude=None).round(n)
    return df


@contextmanager
def timer(title:str):
    """
    https://docs.python.org/3/library/contextlib.html

    """
    import time
    t0 = time.time()
    yield
    logger.info("{} - done in {:.0f}m".format(title, round((time.time()-t0)/60), 2))


def unity_ranger(df:pd.DataFrame, features:list) -> pd.DataFrame:
    """Unity ranger Function

    It truncates to 1 values higher than 1 and to 0 values lower than 0.

    Args:
        df (pd.DataFrame): input data set.
        features (list): Feature list to perform the method

    Returns:
        pd.DataFrame: data set with corrected features.
    """    
    for feature in features:
        df.loc[df[feature] > 1, feature] = 1
        df.loc[df[feature] < 0, feature] = 0
    return df
