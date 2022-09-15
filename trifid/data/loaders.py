#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" trifid/data/loaders.py
https://github.com/fpozoc/trifid/blob/master/trifid/data/loaders.py

Isoform Data Collector

Classes and functions:
    * Fasta
    * GFF
    * GTF
    * load_annotation - returns GTF genome annotations.
    * load_appris - returns APPRIS selected isoforms and scores.
    * load_corsair_alt - returns ALT-Corsair scores.
    * load_corsair_alt_exons - returns ALT-CorsairExons scores.
    * load_phylocsf - returns PhyloCSF scores.
    * load_qpfam - returns QPfam scores.
    * load_qsplice - returns QSplice scores.
    * load_reference - returns annotation type reference categories.
    * load_sequences - returns protein sequences.
    * load_spade - returns APPRIS SPADE scores.
"""

from __future__ import absolute_import, division, print_function

import gzip

import numpy as np
import pandas as pd
from Bio import SeqIO
from gtfparse import read_gtf

from ..utils.utils import get_id_patterns, open_files


class Fasta:
    """
    Management FASTA sequences class
    
    Usage
    >> f = Fasta(path=fasta_path, db=db_name)
    >> df_fasta = f.load
    """    
    def __init__(self, path:str=None, db:str=None):
        self.path = path
        self.db = db

    @property
    def load(self)->list:
        """Loading FASTA sequences file in a pandas DataFrame
        
        Database name has to be inserted in the class declaration. In case of 
        not insert anything, the name of the sequence will be the first word of 
        the identifier.

        Returns:
            list -- pandas DataFrame
        """        
        df = self._to_df()
        return df

    def _load_seqs_dict(self)->dict:
        seqs = list()
        with open_files(self.path) as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                seqs.append({
                    'id': record.description,
                    'sequence': str(record.seq),
                })
        return seqs

    def _to_df(self)->list:
        return pd.DataFrame(self._load_seqs_dict())
    
    def _output(self, fid:str, fseq:str)->str:
        return f">{fid}\n{fseq}\n"


class GFF:
    def __init__(self, path:str, db:str=None):
        self.path = path
        self.db = db
        self.headers = ['seqname', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'tag']

    @property
    def load(self, complete:bool=False)->list:
        df = self._read(complete)
        return df

    def _read(self, extra_attributes:bool=False)->list:
        df = pd.read_csv(self.path, sep='\t', comment='#', compression='gzip', names = self.headers)
        if self.db:
            if self.db in ['gencode', 'g', 'gn']:
                df = self._parse_gencode(df, extra_attributes)
            elif self.db in ['UniProt', 'uniprot', 'u', 'up']:
                pass
            elif self.db in ['RefSeq', 'refseq', 'NCBI', 'r', 'rs']:
                df = self._parse_refseq(df, extra_attributes)
        else:
            pass
        return df

    def _parse_gencode(self, df:list, extra_attributes:bool=False):
        attributes = {
            'main':['ID', 'Parent', 'gene_id', 'transcript_id', 'gene_type', 'gene_name', 
                    'transcript_type','transcript_name', 'exon_id', 'exon_number'],
            'optional':['ont', 'hgnc_id', 'havana_gene', 'havana_transcript', 
                        'transcript_support_level', 'level', 'tag']
                        }
        if extra_attributes:
            attributes = attributes['main'] + attributes['optional']
        else:
            attributes = attributes['main']
        for att in attributes:
            df[att] = df['tag'].str.extract(r'(?:^|;){}=(.*?)(?:$|;)'.format(att))
        return df

    def _parse_refseq(self, df:list, extra_attributes:bool=False):
        attributes = {
            'main': ['ID', 'Parent', 'Dbxref', 'gbkey', 'gene', 'transcript_id', 'gene_biotype', 'protein_id'],
            'optional':['product', 'description', 'exception', 'gene_synonym', 
                        'pseudo', 'model_evidence', 'Name', 'Note', 'inference', 
                        'transl_except', 'function', 'experiment', 'old_locus_tag', 
                        'mobile_element']
                        }
        if extra_attributes:
            attributes = attributes['main'] + attributes['optional']
        else:
            attributes = attributes['main']
        for att in attributes:
            df[att] = df['tag'].str.extract(r'(?:^|;){}=(.*?)(?:$|;)'.format(att))
        df = df.drop('tag', axis=1)
        df['gene_type'] = df['gene_biotype'].fillna(method='ffill')
        df = df.drop('gene_biotype', axis=1)
        df['gene_id'] = df['Dbxref'].str.extract(r'GeneID:(.*?),').fillna(method='ffill')
        df['gene_name'] = df['Parent'].str.split('gene-').str[1].fillna(method='ffill')
        df['transcript_id'] = df['Parent'].str.split('rna-').str[1]
        df['exon_id'] = df['ID'].str.extract(r'^(?:exon|CDS)-(.*?)$')
        df['exon_number'] = df['exon_id'].str.split('-').str[1]
        return df


class GTF:
    """
    Management GTF class
    
    Usage
    >> g = GTF(path, db)
    >> df_gtf = f.load
    """    
    def __init__(self, path:str=None, db:str=None):
        self.path = path
        self.db = db
    
    @property
    def load(self)->list:
        df = self._read()
        return df

    def _read(self)->list:
        from gtfparse import read_gtf
        df = read_gtf(self.path)
        if self.db:
            if self.db in ['gencode', 'g', 'gn']:
                df = self._parse_gencode(df)
            elif self.db in ['UniProt', 'uniprot', 'u', 'up']:
                pass # to do
            elif self.db in ['RefSeq', 'refseq', 'NCBI', 'r', 'rs']:
                df = self._parse_refseq(df)
            elif self.db in ['Ensembl', 'e']:
                df = self._parse_ensembl(df)
        else:
            pass
        return df

    def _parse_gencode(self, df:list, feature:str='transcript')->list:
        headers = ['seqname', 'source', 'feature', 'start', 'end', 'strand',
                    'gene_id', 'gene_name', 'gene_type', 
                    'transcript_id', 'transcript_name', 'transcript_type',
                    'protein_id', 'transcript_support_level', 'level', 'ccdsid',
                    'exon_number', 'exon_id', 
                    'tag']
        if 'transcript_support_level' in df.columns:
            pass
        else:
            df['transcript_support_level'] = -1
        df = df[headers]
        df['appris'] = df['tag'].str.split('appris_').str[1].str.split(',').str[0]
        df.loc[df['tag'].str.contains('start_NF'), 'NF'] = 'SNF'
        df.loc[df['tag'].str.contains('end_NF'), 'NF'] = 'ENF'
        if feature == 'transcript':
            df = df.loc[df['feature'] == feature].reset_index(drop=True)
            df = df.rename(columns={'transcript_support_level': 'tsl'})
            df['tsl'] = df['tsl'].replace(np.nan, 6)
            df['tsl'] = df['tsl'].replace('', 6)
            df['tsl'] = df['tsl'].replace('NA', 6)
            df['CCDS'] = df['ccdsid'].str.contains('CCDS', na=False).astype(int)
            df['StartEnd_NF'] = df['tag'].str.contains(
                'cds_start_NF|cds_end_NF|mRNA_start_NF|mRNA_end_NF', 
                na=False).astype(int)
            df['RNA_supported'] = df['tag'].str.contains(
                'nested_454_RNA_Seq_supported|RNA_Seq_supported_only|RNA_Seq_supported_partial', 
                na=False).astype(int)
            for feature in ['basic', 'CCDS', 'NAGNAG', 'readthrough']:
                df[feature] = df['tag'].str.contains(feature, 
                na=False).astype(int)
            for feature in ['nonsense_mediated_decay', 'non_stop_decay', 'pseudogene']:
                df[feature] = df['transcript_type'].str.contains(feature, 
                na=False).astype(int)
        return df

    def _parse_refseq(self, df:list, feature:str='transcript')->list:
        headers = ['seqname', 'source', 'feature', 'gbkey', 'start', 'end', 'strand',
                    'hgnc_id', 'gene_id', 'gene_name', 'gene_type', 
                    'transcript_id', 'readthrough', 
                    'exon_number', 
                    'product']
        df = df.rename(columns={'gene':'gene_name', 'gene_biotype':'gene_type'})
        df['gene_id'] = df['tag'].str.extract(r'^GeneID:(.*?)(?:$|,)"')
        df['hgnc_id'] = df['tag'].str.extract(r',HGNC:HGNC:(.*?)(?:$|,)"')
        df.loc[df['product'].str.contains('readthrough', na=False), 'readthrough'] = True
        df = df[headers]
        # df = df.loc[df['feature'] == feature]
        return df

    def _parse_ensembl(self, df:list, feature:str='transcript')->list:
        headers = ['seqname', 'source', 'feature', 'start', 'end', 'strand',
                    'gene_id', 'gene_name', 'gene_type', 
                    'transcript_id', 'transcript_name', 'transcript_type',
                    'protein_id', 'transcript_support_level',
                    'exon_number', 'exon_id', 
                    'tag']
        if 'transcript_support_level' in df.columns:
            pass
        else:
            df['transcript_support_level'] = -1
        df = df[headers]
        df.loc[df['tag'].str.contains('start_NF'), 'NF'] = 'SNF'
        df.loc[df['tag'].str.contains('end_NF'), 'NF'] = 'ENF'
        if feature == 'transcript':
            df = df.loc[df['feature'] == feature].reset_index(drop=True)
            df = df.rename(columns={'transcript_support_level': 'tsl'})
            df['tsl'] = df['tsl'].replace(np.nan, 6)
            df['tsl'] = df['tsl'].replace('', 6)
            df['tsl'] = df['tsl'].replace('NA', 6)
            df['CCDS'] = df['tag'].str.contains('CCDS', na=False).astype(int)
            df['StartEnd_NF'] = df['tag'].str.contains(
                'cds_start_NF|cds_end_NF|mRNA_start_NF|mRNA_end_NF', 
                na=False).astype(int)
            for feature in ['basic', 'CCDS']:
                df[feature] = df['tag'].str.contains(feature, 
                na=False).astype(int)
            for feature in ['nonsense_mediated_decay', 'non_stop_decay', 'pseudogene']:
                df[feature] = df['transcript_type'].str.contains(feature, 
                na=False).astype(int)
        return df


def load_annotation(filepath:str, db:str='g')->list:
    """Annotations (.gtf) loader

    GTF Annotation Format https://mblab.wustl.edu/GTF22.html.
    GENCODE GTF docs: https://www.gencodegenes.org/pages/data_format.html.
    RefSeq GTF docs: https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/.

    It takes a GTF file from GENCODE annotation and returns for the selected
    features. It uses the module gtfparse 1.2.0 
    (https://pypi.org/project/gtfparse/) and the function read_gtf to parse all 
    rows of GTF file into a Pandas DataFrame. Several steps are performed to 
    create a final DataFrame which intends to get value from features of 
    transcripts: 
        - Removing the point from transcript id and gene id.
        - Creating a new category of TSL (transcript_support_level) classifying 
        the transcripts without info or with Not Available (NA) info (from 1 to 
        5 in the rest of the transcripts).This category will be labelled as 6.
        - Capturing some features to binarize creating new ones as booleans:
         * StartEnd_NF: If isoform has a Start or End Region that could not be 
        confirmed.
         * basic: If isoform appears in the subset of representative isoforms.
         * NMD: If isoform type is Non-Sense Mediated Decay.
         * NAGNAG: If isoform with Alternative acceptor motifs 3bp like tandem 
        acceptor sites or NAGNAG. In-frame type of variation where, at the 
        acceptor site, some variants splice after the first AG and others after 
        the second AG.
         * CCDS: If isoform is labeled as Consensus Coding DNA Sequence, 
        confirming coding regions between Ensembl, UCSC, NCBI and HAVANA.
         * readthrough: If isoform overlaps two or more independent loci but is 
        considered to belong to a third, separate locus.    

    Args:
        filepath (str): .gz compressed text file path.
        db (str, optional): GENCODE or RefSeq are allowed. Defaults to 'g'.

    Returns:
        list: pandas DataFrame which contains GTF genome annotations.
    """    
    if db=='g':
        df = GTF(filepath, db).load
        if df['transcript_id'].values[0].startswith(get_id_patterns()):
            df['transcript_id'] = df['transcript_id'].str.rsplit(".", 1).str[0]            
        df = df[['transcript_id', 'CCDS', 'StartEnd_NF', 
                'RNA_supported', 'basic', 'NAGNAG', 
                'readthrough', 'nonsense_mediated_decay', 'level']]
    elif db=='e':
        df = GTF(filepath, db).load
        if df['transcript_id'].values[0].startswith(get_id_patterns()):
            df['transcript_id'] = df['transcript_id'].str.rsplit(".", 1).str[0]            
        df["RNA_supported"] = ''
        df["NAGNAG"] = ''
        df["readthrough"] = ''
        df["level_1"] = ''
        df["level_2"] = ''
        df["level_3"] = ''
        df = df[['transcript_id', 'CCDS', 'StartEnd_NF', 
                'RNA_supported', 'basic', 'NAGNAG', 
                'readthrough', 'nonsense_mediated_decay', 'level_1', 'level_2', 'level_3']]
    return df


def load_appris(filepath:str)->list:
    """APPRIS loader
    
    It takes transcript annotations and scores from APPRIS Web Server.
    APPRIS Official doc: http://appris.bioinfo.cnio.es/#/help/intro

        - It retrieves only transcripts flags with TRANSLATION. These are the 
        final reference dataset. 
        - Crash score has 2 integers. To separate them, it creates 2 new 
    categories. One score for signal peptides and the other for mitochondrial 
    signal sequences. 
        - It manages corsair score. Alignments with the same species scores just 
    0.5. Alignments with closest species score 1.5. It is not right for the 
    objective of the study to differenciate between these types of isoforms. 
    To deal with it, we truncate corsair scores lower than 1.5 to 0.

    Args:
        filepath (str): text file path downloaded from APPRIS 
        (http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.appris.txt).

    Returns:
        list: pandas DataFrame which contains APPRIS selected isoforms.
    """    
    headers = {
        'gene_id': 'object',
        'gene_name': 'object',
        'transcript_id': 'object',
        'translation_id': 'object',
        'translation': 'category',
        'flags': 'object',
        'start_stop_codons': 'category',
        'ccdsid': 'object',
        'tsl': 'int64',
        'length': 'float64',
        'firestar': 'float64',
        'matador3d': 'float64',
        'corsair': 'float64',
        'spade': 'float64',
        'thump': 'float64',
        'crash': 'object',
        'inertia': 'category',
        'proteo': 'float64',
        'score': 'float64',
        'appris': 'category', 
    }
    df = pd.read_csv(filepath, sep='\t', names=headers)
    df = df[df['translation'] == 'TRANSLATION']
    df['gene_id'] = df['gene_id'].astype(str)
    df[['crash_p', 'crash_m']] = df['crash'].str.split(',', expand=True)
    df['crash_p'] = df['crash_p'].replace(
        '-', df.loc[~df['crash_p'].str.contains('-$')]['crash_p'].median()
        ).astype('float64')
    df['crash_m'] = df['crash_m'].replace(
        '-', df.loc[~df['crash_m'].str.contains('-$')]['crash_m'].median()
        ).astype('float64')
    df['tsl'] = df['tsl'].replace('-', 6)
    df['corsair'].replace('-', 0, inplace=True)
    if df['corsair'].dtype == object:
        df['corsair'] = df['corsair'].str.replace(r'\$-', '0').astype(float)
    df['corsair'].fillna(0, inplace=True)
    df['corsair'] = df.groupby('gene_id')['corsair'].transform(lambda x: 0 if (x<=1.5).all() else x)
    for name, dtype in headers.items():
        if (dtype == 'float64') and (df[name].dtype == 'object'):
            df.loc[df[name].str.contains(r'[^0-9]', na=False), name] = 0
        df[name] = df[name].astype(dtype)
    if df['transcript_id'].values[0].startswith(get_id_patterns()):
        df['transcript_id'] = df['transcript_id'].str.rsplit(".", 1).str[0]
        df['translation_id'] = df['translation_id'].str.rsplit(".", 1).str[0]
    df = df.drop(
        # ['translation', 'crash', 'inertia', 'start_stop_codons', 'proteo', 'score'], axis=1
        ['translation', 'crash', 'inertia', 'start_stop_codons', 'proteo'], axis=1
        ).reset_index(drop=True)
    return df


def load_corsair_alt(filepath:str)->list:
    """ALT-Corsair loader

    Args:
        filepath (str): .gz compressed text file path.

    Returns:
        list: pandas DataFrame which contains ALT-Corsair scores.
    """    
    df = pd.read_csv(filepath, sep='\t', names=['transcript_id', 'corsair_alt'])
    if df['transcript_id'].values[0].startswith(get_id_patterns()):
        df['transcript_id'] = df['transcript_id'].str.rsplit(".", 1).str[0]
    df['corsair_alt'] = df['corsair_alt'].fillna(0)
    return df


def load_corsair_alt_exons(filepath:str)->list:
    """ALT-CorsairExons loader

    File provided to the function has one exon position per row and include in 
    one column all the isoforms which has these exons. It expand the dataframe 
    to generate one isoform per row with the data available and compress it 
    getting only the minimum per isoform.

    Preproccessing:
    df = pd.read_csv(filepath, sep='\t', names=['label', 'minexon_corsair_alt'])
    df['chr'] = df['label'].str.split('-').str[0].str.split('>').str[1]
    df['chr'] = 'chr' + df['chr']
    df['start'] = df['label'].str.split('-').str[1].str.split(':').str[0]
    df['end'] = df['label'].str.split('-').str[1].str.split(':').str[1].str.split(':').str[0]
    df['transcript_id'] = df['label'].str.split('_').str[1]
    df['frame'] = df['label'].str.split('_').str[0].str.split(':').str[2]

    # GENCODE
    df_gencode_me = df.set_index(
        df.loc[:, df.columns != 'transcript_id'].columns.values.tolist()
        )['transcript_id'].str.split('+', expand=True).stack().reset_index(
            name='transcript_id').drop('level_{}'.format(df.shape[1]-1), axis=1)

    # RefSeq
    df_refseq = read_gtf(filepath)
    df_refseq = df_refseq.loc[df_rs['feature']=='exon']
    df_refseq = df_refseq[['transcript_id', 'start', 'end', 'frame']]
    df_refseq_me  = pd.merge(
        df, df_refseq, on = ['start', 'end', 'frame'], 
        how='left').dropna().drop_duplicates().sort_values(by='transcript_id').reset_index(drop=True)

    Args:
        filepath (str): .gz compressed text file path.

    Returns:
        list: pandas DataFrame which contains ALT-CorsairExons scores.
    """    
    df = pd.read_csv(filepath, sep='\t')
    df = df.groupby('transcript_id')['minexon_corsair_alt'].min().reset_index()
    df['minexon_corsair_alt'] = df['minexon_corsair_alt'].fillna(0)
    if df['transcript_id'].values[0].startswith(get_id_patterns()):
        df['transcript_id'] = df['transcript_id'].str.rsplit(".", 1).str[0]
    return df


def load_phylocsf(filepath: str) -> list:
    """PhyloCSF loader

    Reference: https://academic.oup.com/bioinformatics/article/27/13/i275/178183
    Github: https://github.com/mlin/PhyloCSF/wiki
    Source: https://www.dropbox.com/sh/vrxv1o51s56t3r4/AADB7r7DqjO74YgZIK2m-2n-a?dl=0

    File provided gives data about (mini) exons. File has one row per exon 
    position, consequently, exons shared by 2+ exons appear separated by ";". 
    It has been stacked and divided to represent one row per transcript. 
    Once we have it, the function onlyh takes the minimum value.

    Technically, function extract some features  grouped by isoform, 
    getting a punctuation of the evolutionary for each isoform. 
    File does not contain data for Non Sense Mediated Decay (NMD) variants. 
    Therefore, these isoforms will be imputed as lower values. Score perCodon 
    represents a Codon Substitution Frequencies (CSF) score. Psi represents a 
    different way to represent the CSF. RelBranchLength represents a branch 
    length estimation score.

    Args:
        filepath (str): .gz compressed text file path.

    Returns:
        list: pandas DataFrame which contains PhyloCSF scores.
    """
    df = pd.read_csv(filepath, sep='\t', compression='gzip')
    df = df.loc[df['RelBranchLength'] > 0.1]
    df = df.loc[df['NumCodons'] > 3]
    scorecodon = pd.concat([
        pd.Series(row['ScorePerCodon'], row['Transcripts'].split(','))
        for _, row in df.iterrows()]).reset_index(name='ScorePerCodon').rename(
            columns={'index': 'transcript_id'}
            ).groupby(['transcript_id'])['ScorePerCodon'].min().reset_index()
    relbranch = pd.concat([
        pd.Series(row['RelBranchLength'], row['Transcripts'].split(','))
        for _, row in df.iterrows()]).reset_index(name='RelBranchLength').rename(
            columns={'index': 'transcript_id'}
            ).groupby(['transcript_id'])['RelBranchLength'].min().reset_index()
    psi = pd.concat([
        pd.Series(row['PhyloCSF_Psi'], row['Transcripts'].split(','))
        for _, row in df.iterrows()]).reset_index(name='PhyloCSF_Psi').rename(
            columns={'index': 'transcript_id'}
            ).groupby(['transcript_id'])['PhyloCSF_Psi'].min().reset_index()
    df = pd.merge(scorecodon, relbranch).merge(psi)
    df['transcript_id'] = df['transcript_id'].str.rsplit(".", 1).str[0]
    return df


def load_qpfam(filepath:str)->list:
    """QPfam loader

    https://gitlab.com/fpozoc/qpfam

    Quantifying Pfam effects over reference isoform of every protein-coding gene 
    for the entire human genome.

    Args:
        filepath (str): .gz compressed text file path.

    Returns:
        list: pandas DataFrame which contains qpfam scores.
    """    
    df = pd.read_csv(filepath, compression='gzip', sep='\t')
    df = df[['transcript_id', 'pfam_score', 'pfam_domains_impact_score',
             'perc_Damaged_State', 'perc_Lost_State', 'Lost_residues_pfam',
             'Gain_residues_pfam']]
    return df


def load_qsplice(filepath:str)->list:
    """QSplice loader

    Quantifying splice junctions coverage from data released by STAR and mapping it to genome positions.
    https://gitlab.com/fpozoc/qsplice

    Args:
        filepath (str): .gz compressed text file path.

    Returns:
        list: pandas DataFrame which contains QSplice scores.
    """    
    df = pd.read_csv(filepath, compression='gzip', sep='\t')
    df['RNA2sj'] = df['RNA2sj'].replace('-', np.nan).astype(float).fillna(0)
    df['RNA2sj_cds'] = df['RNA2sj_cds'].replace('-', np.nan).astype(float).fillna(0)
    df['RNA2sj'][df['RNA2sj'] > 1] = 1
    df['RNA2sj_cds'][df['RNA2sj_cds'] > 1] = 1
    df = df[['transcript_id', 'RNA2sj','RNA2sj_cds']]
    if df['transcript_id'].values[0].startswith(get_id_patterns()):
        df['transcript_id'] = df['transcript_id'].str.rsplit(".", 1).str[0]
    return df


def load_reference(filepath:str)->list:
    """Annotation type reference loader

    It provides the labels below to serve as reference correcting the isoforms 
    that actually are fragments (and sometimes they could score better than the
    full transcript).

    Label categories: [Alternative', 'Principal', 'Alternative.NMD', 
    'Redundant Principal', 'Principal Duplication', 'Alternative Duplication', 
    'Redundant Alternative', 'Principal.RT', 'Principal.NMD', 'Alternative.RT']          

    Args:
        filepath (str): .gz compressed text file path.

    Returns:
        list: pandas DataFrame which contains annotation type reference 
    categories.
    """    
    df = pd.read_csv(filepath, compression='gzip', sep='\t', 
    names=['gene_id', 'transcript_id', 'ann_type', 'rtag', 'transcript_ref']) 
    df = df[['transcript_id', 'ann_type', 'transcript_ref']]
    return df


def load_sequences(filepath:str)->list:
    """FASTA sequences loader

    It takes a fasta file from GENCODE/RefSeq genome annotation and returns a 
    pandas DataFrame.

    It uses SeqIO module from BioPython (https://biopython.org/wiki/SeqIO) to 
    parse an standard fasta file and retrive id and seq and parse info from the 
    id.

    Args:
        filepath (str): .gz compressed text file path.

    Returns:
        list: pandas DataFrame which contains protein sequences.
    """    
    df = Fasta(filepath).load
    df.insert(0, 'transcript_id', df['id'].str.split('|').str[1])
    df = df.drop('id', axis=1)
    if df['transcript_id'].values[0].startswith(get_id_patterns()):
        df['transcript_id'] = df['transcript_id'].str.rsplit(".", 1).str[0]
    return df


def load_spade(path:str)->list:
    """Loading SPADE (APPRIS) data.
    http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/e90v24/appris_method.spade.gtf

    Args:
        path (str): SPADE path 

    Returns:
        list: pandas DataFrame
    """    
    headers = {
        'seqname': 'object', 
        'source': 'object', 
        'feature': 'object', 
        'start': 'int64', 
        'end': 'int64', 
        'score': 'float64', 
        'strand': 'object',
        'frame': 'int', 
        'transcript_id': 'object', 
        'gene_id': 'object', 
        'hmm_name': 'object', 
        'evalue': 'float64', 
        'pep_start': 'int64',
        'pep_end': 'int64'
            }
    df = read_gtf(path)
    df['hmm_name'] = df['note'].str.extract(r'^hmm_name:(.*?),')
    df['evalue'] = df['note'].str.extract(r'evalue:(.*?),')
    df['pep_start'] = df['note'].str.extract(r'pep_start:(.*?),')
    df['pep_end'] = df['note'].str.extract(r'pep_end:(.*?)$')
    df = df.drop('note', axis=1)
    for col in df.columns:
        df[col] = df[col].astype(headers[col])
    return df
