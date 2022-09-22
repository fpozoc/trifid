#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""trifid/preprocessing/label_fragments.py

Labelling genome isoforms in duplications, fragments, RT or NMD

Usage: python -m trifid.preprocessing.label_fragments  \
    --gtf data/genome_annotation/GRCh38/g27/gencode.v27.annotation.gtf.gz \
    --seqs data/external/appris/GRCh38/g27/appris_data.transl.fa.gz \
    --principals data/external/appris/GRCh38/g27/appris_data.principal.txt \
    --outdir data/external/label_fragments/GRCh38/g27 \
    --rm
___
--help       | -h    Display documentation
--gtf        | -g    Gtf annotation file.
--outdir     | -o    Output directory.
--principals | -p    APPRIS principals file.
--seqs       | -s    Protein reference file with sequences in fasta format (.gz files allowed).
--rm         | -r    If user wants to remove intermediate files.
--trifid     | -t    TRIFID (optional | previous version) file to sort APPRIS file before label.

Classes and functions:
    * main
    * generate_annotations - returns a custom pandas DataFrame with gtf annotations.
    * generate_sequences - returns a custom pandas DataFrame with protein sequences.
    * get_NR_list - returns the perl command line order to call the script. 
    * get_seqlen - returns the perl command line order to call the script.
"""

from __future__ import absolute_import, division, print_function

import argparse, os, warnings

import pandas as pd
from gtfparse import read_gtf

from ..utils.utils import create_dir, get_id_patterns
from ..data.loaders import Fasta


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--gtf', '-g', 
                        help='Custom or not GENCODE/Ensembl gtf file.', 
                        type=str)
    parser.add_argument('--outdir', '-o',
                        help='Output directory.', 
                        type=str)
    parser.add_argument('--principals', '-p', 
                        help='APPRIS principal list filepath.', 
                        type=str)
    parser.add_argument('--seqs', '-s', 
                        help='APPRIS/GENCODE FASTA sequences.', 
                        type=str)
    parser.add_argument('--rm', '-r',
                        help='If user wants to remove intermediate files.',
                        action='store_true', default=False)
    parser.add_argument('--trifid', '-t',
                        help='TRIFID (optional | previous version) file to sort APPRIS file before label.',
                        type=str, default=None)
    args = parser.parse_args()

    warnings.filterwarnings('ignore')

    create_dir(args.outdir)

    df_fasta = generate_sequences(fasta_path=args.seqs)
    sequences_file = os.path.join(args.outdir, 'appris.pc_sequences.tsv')
    df_fasta.to_csv(sequences_file, index=None, sep='\t', header=None)

    df_gtf = generate_annotations(gtf_path=args.gtf)
    df_gtf.to_csv(os.path.join(args.outdir, 'gencode.pc_annotations.tsv'), index=None, sep='\t', header=None)

    os.system(f"grep 'PRINCIPAL' {args.principals} > {os.path.join(args.outdir, 'appris.principals.tsv')}")
    os.system(get_seqlen(outdir=args.outdir))

    if args.trifid:
        df_trifid = pd.read_csv(args.trifid, sep='\t', compression='gzip')
        df_pc_annotations = pd.read_csv(os.path.join(args.outdir, 'gencode.pc_annotations.out.tsv'), sep='\t', 
                                        names=['gene_id', 'transcript_id', 'sequence', 'length', 'ann_label'])
        df_pc_annotations = pd.merge(df_pc_annotations, df_trifid[['transcript_id', 'appris', 'trifid_score']], on='transcript_id', how='left')
        df_pc_annotations = df_pc_annotations.sort_values(by=['gene_id', 'appris', 'trifid_score'], ascending=[True, True, False])
        df_pc_annotations.to_csv(os.path.join(args.outdir, 'gencode.pc_annotations.out.tsv'), sep='\t', header=None, index=None)

    os.system(get_NR_list(outdir=args.outdir))
    os.system(f"gzip {os.path.join(args.outdir, 'gencode.qduplications.tsv')}")

    if args.rm:
        os.system(f"rm {os.path.join(args.outdir, 'appris.pc_sequences.tsv')}")
        os.system(f"rm {os.path.join(args.outdir, 'appris.principals.tsv')}")
        os.system(f"rm {os.path.join(args.outdir, 'gencode.pc_annotations.out.tsv')}")
        os.system(f"rm {os.path.join(args.outdir, 'gencode.pc_annotations.tsv')}")


def generate_annotations(gtf_path:str) -> pd.DataFrame:
    """It creates a custom pandas DataFrame with gtf annotations.

    Args:
        gtf_path (str): .gtf path.

    Returns:
        pd.DataFrame: gtf annotations output.
    """    
    df = read_gtf(gtf_path)
    df = df.loc[df['feature'] == 'transcript']
    if df['transcript_id'].values[0].startswith(get_id_patterns()):
        df['transcript_id'] = df['transcript_id'].str.rsplit(".", 1).str[0]
        df['gene_id'] = df['gene_id'].str.rsplit(".", 1).str[0]
    df.loc[df.tag.str.contains('readthrough_transcript'), 'readthrough'] = 'readthrough'
    df.loc[df.tag.str.contains('start_NF'), 'NF'] = 'Start NF'
    df.loc[df.tag.str.contains('end_NF'), 'NF'] = 'End NF'
    df = df[df.gene_type == 'protein_coding']
    df = df.loc[df['transcript_type'].str.contains(
        'protein_coding|nonsense_mediated_decay|non_stop_decay|polymorphic_pseudogene|IG|TR')]
    df = df[['gene_name', 'gene_id', 'transcript_id', 'gene_type', 'transcript_type', 
             'readthrough', 'NF', 'transcript_support_level']].reset_index(drop=True)
    return df


def generate_sequences(fasta_path:str) -> pd.DataFrame:
    """It creates a custom pandas DataFrame with protein fasta sequences.

    Args:
        fasta_path (str): .fasta path.

    Returns:
        pd.DataFrame: fasta protein annotations output.
    """    
    df = Fasta(fasta_path).load
    df['gene_id'] = df['id'].str.split('|').str[2]
    df['transcript_id'] = df['id'].str.split('|').str[1]
    if df['transcript_id'].values[0].startswith(get_id_patterns()):
        df['gene_id'] = df['gene_id'].str.rsplit(".", 1).str[0]
        df['transcript_id'] = df['id'].str.split('|').str[1]
        df['transcript_id'] = df['transcript_id'].str.rsplit(".", 1).str[0]
    df['length'] = df['id'].str.split('|').str[-1]
    df = df[['gene_id', 'transcript_id', 'sequence']]
    return df


def get_NR_list(outdir:str) -> str:
    """Perl script caller 

    Args:
        outdir (str): Directory to store the output.

    Returns:
        str: command line order to run the program.
    """    
    perl_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../utils', 
                            'get_NR_list.pl')
    cmd = f'perl {perl_script} {outdir}'
    return cmd


def get_seqlen(outdir:str)->str:
    """Perl script caller 

    Args:
        outdir (str): Directory to store the output.

    Returns:
        str: command line order to run the program.
    """    
    perl_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../utils', 
                            'get_seqlen.pl')
    cmd = f'perl {perl_script} {outdir}'
    return cmd


if __name__ == "__main__":
    main()