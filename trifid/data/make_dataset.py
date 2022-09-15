#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" trifid/data/make_dataset.py
https://github.com/fpozoc/trifid/blob/master/trifid/data/make_dataset.py

Making the initial TRIFID data set

Usage: python -m trifid.make_dataset \
    --config config/config.yaml \
    --features config/features.yaml \
    --assembly GRCh38 \
    --release g27
___
--help      | -h    Display documentation
--assembly  | -a    Genome assembly version
--config    | -c    Configuration filepath
--features  | -f    Features selected description filepath
--release   | -h    Genome release version
"""

from __future__ import absolute_import, division, print_function

import argparse, os, warnings

import pandas as pd
from loguru import logger

from .feature_engineering import build_features, load_data
from ..utils import utils


def main():    
    parser = argparse.ArgumentParser(description='Command-line arguments parser', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--assembly',
                    default='GRCh38', help="Genome assembly version.")
    parser.add_argument('-c', '--config',
                    default='config/config.yaml', help="Configuration filepaths.")
    parser.add_argument('-f', '--features',
                    default='config/features.yaml', help="Features selected description filepath.")
    parser.add_argument('-r', '--release',
                    default='g27', help="Genome release version.")
    args = parser.parse_args()

    warnings.filterwarnings("ignore")
    
    config = utils.parse_yaml(args.config)
    df_features = pd.DataFrame(utils.parse_yaml(args.features))
    data_dir = os.path.join('data', 'genomes', args.assembly, args.release)
    utils.create_dir(data_dir)
    logger.add(os.path.join(data_dir, 'trifid_db.{time}.log'))
    logger.info(f'TRIFID has started and its output will be ready in {data_dir}')

    df = load_data(config=config, assembly=args.assembly, release=args.release)
    df = build_features(df)
    df[df_features.feature.values].to_csv(
        os.path.join(data_dir, 'trifid_db.tsv.gz'), index=None, sep='\t', compression='gzip')

if __name__ == "__main__":
    main()