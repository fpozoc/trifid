#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""trifid/preprocessing/qsplice.py

Quantifying splice junctions coverage from data released by STAR and mapping it to genome positions.

Usage: python -m trifid.preprocessing.qsplice \
    --gff   ~/hdd1/data/genome_annotation/GRCh38/g27/gencode.v27.annotation.gff3.gz \
    --outdir data/external/qsplice/GRCh38/g27 \
    --samples /home/fpozoc/hdd2/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g27 \
    --version g
___
--help         | -h    Display documentation
--gff          | -g    Gff annotation file.
--outdir       | -o    Genome annotation version. GENCODE: `g` + `nversion`.
--samples      | -s    Customized SJ file.
--rm           | -r    If user wants to remove intermediate files.
--version      | -v    Directory which contains the files to be globbed and concatenated.

Classes and functions:
    * annotate_introns - returns introns annotated with position and CDS coverage.
    * concat_samples - returns a processed set RNA-seq SJ.out.tab samples.
    * generate_introns - returns an annotation file with introns positions generated.
    * load_annotations - returns a DataFrame with annotations processed.
    * map_junctions_positions - returns the coverage of splice junction genome positions.
    * score_per_junction - returns a score per splice junction.
    * score_per_transcript - returns a score per transcript.
"""

from __future__ import absolute_import, division, print_function

import argparse, glob, os, re, subprocess

from loguru import logger
import pandas as pd
import numpy as np

from ..data.loaders import GFF
from ..utils.utils import create_dir, get_id_patterns


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--experiment', '-e',
                        help='Output directory.', 
                        type=str, default='emtab2836')
    parser.add_argument('--file', '-f', 
                        help='Custom splice junctions file (in gzip)',
                        type=str)
    parser.add_argument('--gff', '-g', 
                        help='Custom or not GENCODE/Ensembl gff file.', 
                        type=str)
    parser.add_argument('--outdir', '-o',
                        help='Output directory.', 
                        type=str)
    parser.add_argument('--rm', '-r',
                        help='If user wants to remove intermediate files.',
                        action='store_true', default=False)
    parser.add_argument('--samples', '-s', 
                        help='Directory which contains the files to be globbed and concatenated', 
                        type=str)
    parser.add_argument('--tissue_ids', '-t', 
                        help='Directory which contains the files to be globbed and concatenated', 
                        type=str)
    parser.add_argument('--version', '-v',
                        help='Genome annotation version. GENCODE: `g` + `nversion`',
                        type=str)
    args = parser.parse_args()
    
    logger.info(f"Program has been launched succesfully.")

    create_dir(args.outdir)

    gff_path = generate_introns(args.gff)
    logger.info(f'Introns generated.')
    df_annotations = load_annotations(gff_path, db=args.version[0])
    logger.info(f'Annotations generated.')
    df_introns_annotated = annotate_introns(df_annotations)
    logger.info(f'CDS coverage and complete database generated.')
    
    df_introns = df_introns_annotated[df_introns_annotated['type'] == 'intron'].drop(
        ['exon_id','exon_number','start_cds','end_cds'], axis=1)
    logger.info(f'Introns annotation generated: {df_introns.shape[0]}')
    logger.info(f'Introns coding coverage: {df_introns.cds_coverage.value_counts().to_dict()}')

    exons_annotation_path = gff_path.replace('gz', 'exons_cds.tsv.gz')
    df_exons = df_introns_annotated[~df_introns_annotated['type'].str.contains('intron')].drop(
        ['intron_number'], axis=1)
    df_exons.to_csv(exons_annotation_path, index=None, sep='\t', compression='gzip')
    logger.info(f'Exons annotation generated: {df_introns.shape[0]}')
    logger.info(f'Exons coding coverage {df_exons.cds_coverage.value_counts().to_dict()}')

    if args.file:
        # cat STAR/g27/ERR*/SJ.out.tab > SJ.out.tab.concat && gzip SJ.out.concat 
        df_sj = pd.read_csv(args.file, compression='gzip', sep='\t', 
                            names=['seqname', 'start', 'end', 'nstrand', 'unique_reads', 'tissue'])
        df_sj = df_sj.drop('nstrand', axis=1)
    else:
        df_sj = concat_samples(indir=f'{args.samples}/*/SJ.out.tab', tissue_ids=args.tissue_ids) # concatenating several SJ.out after being annotated and parsed in same pandas DataFrame
        df_sj['tissue'] = df_sj['tissue'].str.split('_').str[0]  # Getting max values per position and per tissue group sample
    df_sj_max_position, df_sj_max_tissue =  map_junctions_positions(df_sj)

    if args.rm == False:
        df_sj_max_position.to_csv(os.path.join(args.outdir, f'sj_maxp.{args.experiment}.tsv.gz'), index=None, sep='\t', compression='gzip')
        df_sj_max_tissue.to_csv(os.path.join(args.outdir, f'sj_maxt.{args.experiment}.tsv.gz'), index=None, sep='\t', compression='gzip')

    df_junction_score = score_per_junction(df_introns, df_sj_max_position)
    df_junction_score.to_csv(os.path.join(args.outdir, f'sj_maxp.emtab2836.mapped.tsv.gz'), sep='\t', index=None, compression='gzip')

    df_qsplice = score_per_transcript(df_junction_score)
    df_qsplice.to_csv(os.path.join(args.outdir, f'qsplice.{args.experiment}.tsv.gz'), sep='\t', index=None, compression='gzip')


def annotate_introns(df: pd.DataFrame) -> pd.DataFrame:
    """Annotating introns cds positions merging exons and cds info.

    Arguments:
        df (pd.DataFrame) -- Input from whole introns annotation.

    Returns:
        pd.DataFrame -- output DataFrame with introns annotations.
    """
    df_whole = _add_intron_number(df)
    df_exons = _get_exons_cds(df_whole)
    df_introns = df_whole[df_whole['type'] == 'intron'].reset_index(drop=True)
    df_annot = pd.concat([df_exons, df_introns], axis=0).reset_index(drop=True)
    del df_introns
    del df_exons
    del df_whole
    # df_annot = df_annot[df_annot['gene_type'] == 'protein_coding']
    df_annot['nexons'] = df_annot.groupby('transcript_id')['type'].transform(lambda x: x.str.contains('exon').sum())
    df_annot['ncds'] = df_annot.groupby('transcript_id')['cds_coverage'].transform(lambda x: x.str.contains('full|partial').sum())
    df_annot = df_annot.sort_values(by=['gene_id', 'transcript_id', 'start', 'end'], ascending=[False, True, True, True]).reset_index(drop=True)
    df_annot.loc[df_annot['cds_coverage'].isnull(), 'l_intron_coverage'] = df_annot.groupby('transcript_id')['cds_coverage'].shift(1)
    df_annot.loc[df_annot['cds_coverage'].isnull(), 'r_intron_coverage'] = df_annot.groupby('transcript_id')['cds_coverage'].shift(-1)
    df_annot.loc[(df_annot['cds_coverage'].isnull()) & (df_annot['l_intron_coverage'].str.contains('none') | df_annot['r_intron_coverage'].str.contains('none')), 'cds_coverage'] = 'none'
    df_annot.loc[(df_annot['cds_coverage'].isnull()) & (df_annot['l_intron_coverage'].str.contains('full|partial') | df_annot['r_intron_coverage'].str.contains('full|partial')), 'cds_coverage'] = 'full'
    df_annot = df_annot.drop(['l_intron_coverage', 'r_intron_coverage'], axis=1)
    return df_annot


def concat_samples(indir: str, tissue_ids:str) -> pd.DataFrame:
    """Concatenate several SJ.out.tab files already processed in same pandas DataFrame.

    Args:
        indir (str): annotation file directory.
        sample_ids (str): samples identifiers.

    Returns:
        pd.DataFrame: SJ.out.tab concatenated with tissue annotations.
    """    
    globbed_dir = glob.glob(f'{indir}')
    annotation_dict = _parse_emtab(tissue_ids)
    df_samples = pd.concat([_process_sj(filepath, annotation_dict) for filepath in globbed_dir]).reset_index(drop=True)
    logger.info(f'{len(globbed_dir)} RNA-seq samples loaded and concatenated.')
    return df_samples


def generate_introns(inpath: str) -> str:
    """Generating introns file from gff*.gz file path.
    
    tidy option described here http://genometools.org/pipermail/gt-users/2015-August/000794.html
    
    Arguments:
        inpath (str) -- The name of the annotation file to process.

    Returns:
        str -- Path where file has been stored.
    """
    outpath = re.sub(r'gff', r'introns.gff', inpath)
    # Ensembl hack
    # sed -i "s/Parent=transcript:/Parent=/g" {inpath}' 
    # sed -i "s/ID=transcript:/ID=/g" {inpath}'
    cmd_gt = f"zcat {inpath} | gt gff3 -tidy -retainids -addintrons | gzip > {outpath}" # 
    subprocess.call(cmd_gt, shell=True)
    return outpath


def load_annotations(filepath:str, db:str, features:list=['CDS', 'exon', 'intron']) -> pd.DataFrame:
    """Loading features from genome annotation file.

    Args:
        filepath (str): .gff annotation file path.
        db (str): Annotation reference (GENCODE, RefSeq or UniProt allowed).
        features (list): Feature annotations selected to filter.

    Returns:
        pd.DataFrame: Features from gff.
    """    
    df_gff = GFF(filepath, db=db).load
    df_gff.loc[df_gff['type'] == 'intron', 'transcript_id'] = df_gff['Parent']
    df_gff = df_gff.loc[df_gff['type'].str.contains('|'.join(features))]
    if db.lower().startswith('r'): # RefSeq annotations
        df_gff['exon_id'] = df_gff['ID'].str.extract(r'^(?:exon|CDS)-(.*?)$')
        df_gff['exon_number'] = df_gff['exon_id'].str.split('-').str[1]
        df_gff.loc[df_gff['type'] == 'intron', 'transcript_id'] = df_gff['Parent'].str.split('rna-').str[1]
    df_gff = df_gff[['seqname', 'type', 'start', 'end', 'strand', 'gene_id', 
    'gene_name', 'transcript_id', 'exon_id', 'exon_number']]
    logger.info(f"{df_gff.shape[0]} annotations loaded and stored with this columns:\n{';'.join(df_gff.columns)}")
    return df_gff


def map_junctions_positions(df:pd.DataFrame) -> pd.DataFrame:
    """This function use a pandas DataFrame with concatenated RNA-seq (SJ.out.tab)
    concatenated to extract the highest coverage per position and per tissue and 
    position.

    Arguments:
        df {list} -- pandas DataFrame splice junctions and reads concatenated.

    Returns:
        df_sj_max_position {list} -- pandas DataFrame with maximum coverage per junction position.
        df_sj_max_tissue {list} -- pandas DataFrame with maximum coverage per junction position and tissue.
    """    
    df_sj_max_position = df.sort_values(
        by=['start', 'end', 'unique_reads'], 
        ascending=[True, True, False]).drop_duplicates(
            subset=['start', 'end'], keep='first').reset_index(drop=True)
    logger.info(f"Mean unique reads per position: {round(df_sj_max_position['unique_reads'].mean(), 3)}")

    df_sj_max_tissue = df.sort_values(
        by=['start', 'end', 'unique_reads'], 
        ascending=[True, True, False]).drop_duplicates(
            subset=['start', 'end', 'tissue'], keep='first').reset_index(drop=True)
    logger.info(f"Mean unique reads per tissue and position: {round(df_sj_max_tissue['unique_reads'].mean(), 3)}")
    return df_sj_max_position, df_sj_max_tissue


def score_per_junction(df_introns:list, df_sj_max_position:list)->list:
    """Calculating means and score per junction.

    Arguments:
        df_introns {list} -- pandas DataFrame with introns.
        df_sj_max_postion {list} -- pandas DataFrame with junction read positions.

    Returns:
        list -- pandas DataFrame with score per exon, gene and transcript.
    """    
    df = pd.merge(df_introns, df_sj_max_position, how='left', on=['seqname', 'start', 'end'])
    df.loc[df['unique_reads'].isnull(), 'unique_reads'] = 0
    df.loc[df['unique_reads'].isnull(), 'tissue'] = '-'
    df['gene_mean'] = df.groupby('gene_id')['unique_reads'].transform(lambda x:x.mean())
    df['gene_mean'] = df['gene_mean'].fillna(df.groupby('gene_id')['gene_mean'].transform('mean'))
    df['gene_mean_cds'] = df.loc[df['cds_coverage'] == 'full'].groupby('gene_id')['unique_reads'].transform(lambda x:x.mean())
    df['gene_mean_cds'] = df['gene_mean_cds'].fillna(df.groupby('gene_id')['gene_mean_cds'].transform('mean'))
    df['RNA2sj'] = df['unique_reads']/df['gene_mean']
    df['RNA2sj_cds'] = df['unique_reads']/df['gene_mean_cds']
    df[df.columns[-6:]] = round(df[df.columns[-6:]], 4)
    if df['transcript_id'].values[0].startswith(get_id_patterns()):
        df['gene_id'] = df['gene_id'].str.rsplit(".", 1).str[0]
        df['transcript_id'] = df['transcript_id'].str.rsplit(".", 1).str[0]
    df = df.sort_values(by=['seqname', 'gene_id', 'transcript_id', 'start', 'end'])
    logger.info(f'{df.shape[0]} junctions from protein-coding genes loaded and quantified.')
    return df


def score_per_transcript(df:pd.DataFrame) -> pd.DataFrame:
    """Calculating mins and exporting qsplice

    Arguments:
        df (pd.DataFrame) -- Score per junction.

    Returns:
        pd.DataFrame -- Scores per transcript (qsplice output).
    """    
    df['cds_coverage'] = df.groupby('transcript_id')['cds_coverage'].transform(
        lambda x: 'null' if (x.str.contains('none')).all() else x)
    df.loc[df['cds_coverage'] == 'none', 'unique_reads'] = df['unique_reads'].max()
    df = df.loc[df.groupby('transcript_id')['unique_reads'].idxmin()].sort_values(
        by=['seqname', 'gene_id', 'transcript_id', 'start', 'end'], 
        ascending=[True, True, True, True, True]).reset_index(drop=True)
    df = df.drop(['type', 'cds_coverage', 'start', 'end', 'strand'], axis=1)
    if df['transcript_id'].values[0].startswith(get_id_patterns()):
        df['gene_id'] = df['gene_id'].str.rsplit(".", 1).str[0]
        df['transcript_id'] = df['transcript_id'].str.rsplit(".", 1).str[0]
    df = df.sort_values(
        by=['seqname', 'gene_id', 'transcript_id', 'RNA2sj_cds'], 
        ascending=[True, True, True, False])
    logger.info(f'{df.shape[0]} transcripts from protein-coding genes loaded and quantified.')
    return df


def _add_intron_number(df: pd.DataFrame) -> pd.DataFrame:
    """Adding intron number and correct positions for CDS. This function sorts 
    the positions in order to get the proper exon-CDS order in both 5' and 3', 
    and below it adds the intron number and correct CDS positions correcting by 
    strand.

    Arguments:
        df (pd.DataFrame) -- Initial (and previous) introns annotation.

    Returns:
        pd.DataFrame -- Introns context added.
    """
    df = df.sort_values(
        by=['seqname', 'transcript_id', 'start', 'end'], 
        ascending=[True, True, True, False]) # sorting the positions in order to get the proper exon-CDS order in both 5' and 3'.
    df['intron_number'] = df['exon_number']
    intronf_list = ['gene_id','transcript_id']
    # intronf_list = ['gene_id','gene_type','transcript_id']
    for feature in intronf_list:
        df.loc[df['type'].str.contains('exon|intron'), feature] = df.fillna(method='ffill')
    df.loc[df['type'].str.contains('exon|intron'), 'intron_number'] = df.fillna(method='ffill')
    df.loc[~df['type'].str.contains('intron'), 'intron_number'] = '-'
    df.loc[df['type'].str.contains('exon|CDS'), 'exon_number'] = df.fillna(method='ffill')
    df.loc[df['type'].str.contains('exon|CDS'), 'exon_id'] = df.fillna(method='ffill')
    return df


def _get_exons_cds(df: pd.DataFrame) -> pd.DataFrame:
    """Getting exons with cds annotations.

    Arguments:
        df (pd.DataFrame) -- Initial (and previous) exons and introns annotation.

    Returns:
        pd.DataFrame -- Exons annotation.
    """

    exon = df[df['type'] == 'exon'].drop(['intron_number', 'type'], axis=1).reset_index(drop=True)
    cds = df[df['type'] == 'CDS'].drop(['intron_number', 'type'], axis=1).reset_index(drop=True)
    df = pd.merge(exon, cds, how='left',
                     on=list(df.drop(['start', 'end', 'type', 'intron_number'], axis=1).columns),
                     indicator='cds_coverage')
    df['cds_coverage'] = df['cds_coverage'].replace('left_only', 'none')
    df['cds_coverage'] = df['cds_coverage'].replace('both', 'full')
    df.insert(1, 'type', 'exon')
    df = df.rename(columns={'start_x':'start',
                            'end_x': 'end',
                            'start_y': 'start_cds',
                            'end_y': 'end_cds'})

    df.loc[df['start_cds'].isnull(), 'start_cds'] = df['start']
    df.loc[df['end_cds'].isnull(), 'end_cds'] = df['end']
    df['start_cds'] = df['start_cds'].astype(int)
    df['end_cds'] = df['end_cds'].astype(int)
    df['cds_coverage'] = df['cds_coverage'].astype(str)
    df.loc[(df['start_cds'] != np.nan) & ((df['start_cds'] != df['start']) | (df['end_cds'] != df['end'])), 'cds_coverage'] = 'partial'
    df = df.drop_duplicates(['start', 'end', 'seqname', 'gene_id', 'transcript_id'])
    return df


def _parse_emtab(inpath: str) -> dict:
    '''Parsing experiment data

    source: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2836/samples/
    It takes a table as input and returns a dictionary with identifier as key 
    and tissue as value.

    Args:
        inpath (str): annotation file path.

    Returns:
        samples_ids (dict): experiments identifiers with tissues as values.
    '''
    df = pd.read_csv(inpath, sep ='\t')
    df['Comment[ENA_RUN]'] = df['Comment[ENA_RUN]'].apply(lambda x: x + '.1') # Adds '.1' to sample
    df  = df[['Comment[ENA_RUN]', 'Source Name']].drop_duplicates(
        subset='Comment[ENA_RUN]').sort_values(by='Comment[ENA_RUN]'
        ).reset_index(drop=True) # Dropping duplicates names
    samples_ids = dict(zip(df['Comment[ENA_RUN]'], df['Source Name']))
    return samples_ids


def _process_sj(inpath: str, samples_ids: dict = None) -> list:
    '''SJ.out.tab processor. Removing non mapped regions and unannotated maps.

    STAR manual
    https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf

    Args:
        inpath (str): SJ.out.tab infile.
        samples_ids (dictionary): Select if you want to add a tissue annotation for that SJ file. In this case introduce a dictionary with dirname identifier as keys and tissue names as values.

    Returns:
        df (list): pandas DataFrame with file processed.
    '''
    df = pd.read_csv(
        inpath, sep='\t', 
        names=['seqname', 'start', 'end', 'nstrand', 'intron-motif', 'annotated', 'unique_reads', 'multimapping_reads', 'overhang'])
    df = df[df['annotated'] == 1][['seqname', 'start', 'end', 'nstrand', 'unique_reads']]
    if samples_ids:
        df['tissue'] = samples_ids[os.path.basename(os.path.dirname(inpath))] # adding tissue
    df = df.drop('nstrand', axis=1)
    return df


if __name__ == "__main__":
    main()