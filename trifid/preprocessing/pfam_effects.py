#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""trifid/preprocessing/pfam_effects.py

Quantifying Pfam effects over reference isoform of every protein-coding gene for 
the entire genome

Usage: python -m trifid.preprocessing.pfam_effects \
    --appris ~/hdd1/data/appris/GRCh38/g27/appris_data.appris.txt \
    --jobs 10 \
    --seqs ~/hdd1/data/appris/GRCh38/g27/appris_data.transl.fa.gz \
    --spade ~/hdd1/data/appris/GRCh38/g27/appris_method.spade.gtf.gz \
    --outdir data/external/pfam_effects/GRCh38/g27 \
    --rm
___
--help      | -h    Display documentation.
--appris    | -a    APPRIS scores data set.
--jobs      | -j    List of cpus used to run the programs in parallel.
--outdir    | -o    Output directory.
--rm        | -r    If user wants to remove intermediate files.
--seqs      | -s    Protein reference file with sequences in fasta format (.gz files allowed).
--spade     | -p    Spade file to be included in the pipeline.


Classes and functions:
    * main
    * analyse_transcripts - returns the perl command line order to call the script.
    * annotation_reference - returns the reference transcript per gene.
    * load_pfam - returns the `pfam_effects` processed file.
    * make_spade - returns the SPADE DataFrame.
    * mp_msa
    * qpfam_effects - returns the TRIFID scores for this module.
"""

from __future__ import absolute_import, division, print_function

import argparse
import functools
import multiprocessing as mp
import os
import warnings

import pandas as pd
from loguru import logger

from ..data.loaders import Fasta, load_appris, load_spade
from ..utils.utils import create_dir, get_id_patterns, merge_dataframes, open_files


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--appris",
        "-a",
        help="It could be necessary to add an APPRIS path file for some version incompatibilities with GENCODE.",
        type=str,
    )
    parser.add_argument("--jobs", "-j", help="List of cpus used to run the programs in parallel.", type=int)
    parser.add_argument("--outdir", "-o", help="Output directory.", type=str)
    parser.add_argument(
        "--rm", "-r", help="If user wants to remove intermediate files.", action="store_true", default=False
    )
    parser.add_argument(
        "--seqs", "-s", help="Protein reference file with sequences in fasta format (.gz files allowed).", type=str
    )
    parser.add_argument("--spade", "-p", help="Spade file to be included in the pipeline.", type=str)
    args = parser.parse_args()

    warnings.filterwarnings("ignore")

    create_dir(args.outdir)

    logger.info(f"Program has been launched succesfully.")

    mp_msa(appris_path=args.appris, fasta_path=args.seqs, outdir=args.outdir, cpus=args.jobs)
    logger.info("""Muscle MSA has been generated.""")

    df_spade = make_spade(spade_path=args.spade, outdir=args.outdir)
    logger.info("""SPADE files has been generated.""")

    df_appris_msa = annotation_reference(args.appris, args.seqs, args.outdir, save=True)
    pfam_effects_file = os.path.join(args.outdir, "Pfam_effects.tsv")
    os.system(analyse_transcripts(outdir=args.outdir, outfile=pfam_effects_file))
    logger.info(f"""The pfam effects has been generated.""")

    df_qpfam = qpfam_effects(df_reference=df_appris_msa, df_spade=df_spade, pfam_effects_filepath=pfam_effects_file)
    logger.info(f'{df_qpfam["transcript_id"].nunique()} transcripts quantified.')
    df_qpfam = df_qpfam.drop_duplicates("transcript_id")
    if df_qpfam["transcript_id"].values[0].startswith(get_id_patterns()):
        df_qpfam["transcript_id"] = df_qpfam["transcript_id"].str.replace(".", "_", 1)
    else:
        df_qpfam["transcript_id"] = df_qpfam["transcript_id"].str.replace("_", ".")
    df_qpfam.to_csv(os.path.join(args.outdir, "qpfam.tsv.gz"), sep="\t", index=None, compression="gzip")
    os.system(f"rm -rf {os.path.join(args.outdir, 'tmp*')}")

    if args.rm:
        os.system(f"rm -rf {os.path.join(args.outdir, 'muscle')}")
        os.system(f"rm -rf {os.path.join(args.outdir, 'spade3')}")
        os.system(f"rm {os.path.join(args.outdir, 'Pfam_effects.tsv')}")
        os.system(f"rm {os.path.join(args.outdir, 'spade_references.tsv')}")
        os.system(f"rm {os.path.join(args.outdir, 'ListOfSawnOffTranscripts.tsv')}")
        os.system(f"rm {os.path.join(args.outdir, 'pc_translations.tsv')}")


def analyse_transcripts(outdir: str, outfile: str) -> str:
    """Perl script caller to analyse Muscle MSA (ClustalW type) and quantify the
    effect over Pfam domain of every transcript.

    Args:
        outdir (str): Directory to store the output.
        outfile (str): Path of the output file.

    Returns:
        str: command line string to run the program.
    """
    perl_script = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "../utils", "analyse_appris_spade_transcripts_nf.pl"
    )
    cmd = f"perl {perl_script} {outdir} {outfile}"
    return cmd


def annotation_reference(appris_path: str, fasta_path: str, outdir: str = None, save: bool = None):
    """Getting one transcript per gene as reference. Priority will be selected by:

    1. Protein coding label.
    2. Best SPADE (APPRIS score).
    3. TSL 1 label.
    4. CCDS id label.
    5. Highest number of residues.
    6. CCDS lowest number.
    7. APPRIS tag.

    Args:
        appris_path (str): APPRIS path.
        fasta_path (str): FASTA sequences path.
        save (bool, optional): If user wants to save the reference. Defaults to None.
    """
    df_appris = load_appris(appris_path)
    df_fasta = Fasta(fasta_path).load
    df_fasta["transcript_id"] = df_fasta["id"].str.split("|").str[1]
    if df_fasta["transcript_id"].values[0].startswith(get_id_patterns()):
        df_fasta["transcript_id"] = df_fasta["transcript_id"].str.rsplit(".", 1).str[0]
    df_reference = pd.merge(df_appris, df_fasta[["transcript_id", "sequence"]], how="left", on="transcript_id")
    df_reference["transcript_id"] = df_reference["transcript_id"].str.replace(".", "_")
    df_priority = _select_reference_priority(df_reference)
    df_pfam_effects_reference = df_priority.drop_duplicates("gene_id")[
        ["gene_name", "gene_id", "transcript_id", "appris"]
    ].reset_index(drop=True)
    df_reference.loc[df_reference["transcript_id"].isin(df_priority["transcript_id"]), "pfam_effects_msa"] = "Reference"
    df_reference.loc[
        ~df_reference["transcript_id"].isin(df_priority["transcript_id"]), "pfam_effects_msa"
    ] = "Transcript"
    df_appris_msa = df_reference[["gene_id", "transcript_id", "pfam_effects_msa", "appris", "sequence"]]

    if save:
        df_pfam_effects_reference.to_csv(
            os.path.join(outdir, "spade_references.tsv"), sep="\t", index=None, header=None
        )
        df_appris_msa[["gene_id", "transcript_id", "sequence"]].to_csv(
            os.path.join(outdir, "pc_translations.tsv"), sep="\t", index=None, header=None
        )

    return df_appris_msa


def load_pfam(file: str) -> pd.DataFrame:
    """File provided to the function contains data from pfam domains in isoform
    structures.

    A single row represents the effect of AS over pfam domains. Initial file has
    data for pfam domain number, name, start and end. Also it adds the type of
    event, the lost residues over the protein, and the event residues over the
    pfam domain. It also classifies the State of the event (Damaged, Lost or
    Intact), and the event type (C or N terminal swap, Homology, Deletion,
    Insertion, C or N terminal Deletion, NAGNAG, Two Proteins and Substitution)
    calculating a pfam score, which represents the impact over the isoform,
    being the lowest the one with less impact over structure.

    Then, it separates transcripts because more than one transcript can
    share same effect. Moreover, it converts pfam_score to 0-1 scale being the 1
    the less affected by AS and 0 an isoform totally damaged. It represents the
    quantitative impact over residues of domains. Therefore, it calculates
    another score getting the worst or more damaged event per isoform having
    higher punctuation an instance with all residues intacts than the one which
    has lost/damaged its pfam domains.

    Furthermore, it counts impact generated by event type, taking into account
    that State is either lost or damage per isoform. It also counts the State of
    each event summarized by isoform. It gets in an extra column the number of
    events per transcript. Percentages of these counts divided by events per
    transcript were also calculated.

    Finally, the last predictive feature created was pfam_domains_impact_score.
    With data available, it creates the score taking number of domains per
    isoform and adding the ones which are classified as Lost or Damaged and
    substrating the other ones that are classified as Intact. Then, it tries to
    get a value which represents the integrity of the structure, but instead of
    previous score it performs the operation in a different way.

    Args:
        file (str): pfam effects source file.

    Returns:
        list: pandas DataFrame with potentially predictive structure features.
    """
    df = pd.read_csv(
        file,
        sep="\t",
        names=[
            "gene_id",
            "alternative_principal(;sep)",
            "Pfam_domain_number",
            "Pfam_domain_name",
            "Pfam_domain_start",
            "Pfam_domain_end",
            "Event_type",
            "Lost_residues_total",
            "Lost_residues_pfam",
            "State",
            "Identity_percent",
            "Gaps_percent",
            "pfam_score",
        ],
    )
    df["gene_id"] = df["gene_id"].astype(str)
    df = df.loc[~df["Event_type"].isnull()]
    df["transcript_id"] = df["alternative_principal(;sep)"].str.split(" ").str[0]
    df["principal"] = df["alternative_principal(;sep)"].str.split(" ").str[1]
    df["pfam_score"] = 1 - (df["pfam_score"] / 100)
    impacts = _instance_counter(df, "Event_type")
    states = _instance_counter(df, "State")
    nevents = _n_counter(df, "Pfam_domain_number")
    minpfamscore = df.groupby(["transcript_id"])["pfam_score"].min().fillna(1).reset_index()
    lost_residues_total = _unique_counter(df, "Lost_residues_total", "Lost")
    gain_residues_total = _unique_counter(df, "Lost_residues_total", "Gain")
    lost_residues_pfam = _unique_counter(df, "Lost_residues_pfam", "Lost")
    gain_residues_pfam = _unique_counter(df, "Lost_residues_pfam", "Gain")
    lost_gain = merge_dataframes(
        lost_residues_total,
        gain_residues_total,
        lost_residues_pfam,
        gain_residues_pfam,
        on_type="transcript_id",
        how_type="outer",
        pivot_on=0,
        nimpute=0,
    )
    df = merge_dataframes(
        impacts, states, nevents, minpfamscore, lost_gain, on_type="transcript_id", how_type="left", pivot_on=0
    )
    df.loc[df["Lost_State"] >= 1, "pfam_score"] = 0
    return df


def make_spade(spade_path: str, outdir: str) -> pd.DataFrame:
    """Splitting spade transcript in several files.

    Args:
        spade_path (str): APPRIS SPADE file path.
        outdir (str): Directory to store the list of pfam domain files.

    Returns:
        pd.DataFrame: pandas DataFrame with SPADE method preprocessed.
    """
    if spade_path.endswith(".gtf.gz"):
        df_spade = load_spade(spade_path)
    elif spade_path.endswith(".txt"):
        df_spade = pd.read_csv(
            spade_path, names=["transcript_id", "hmm_name", "pep_start", "pep_end", "score"], sep="\t"
        )
        df_spade.insert(4, "feature", "-", allow_duplicates=False)
    if df_spade["transcript_id"].values[0].startswith(get_id_patterns()):
        df_spade["transcript_id"] = df_spade["transcript_id"].str.rsplit(".", 1).str[0]
    df_spade["transcript_id"] = df_spade["transcript_id"].str.replace(".", "_")
    df_spade = df_spade[["transcript_id", "hmm_name", "pep_start", "pep_end", "feature", "score"]]
    spadedir = os.path.join(outdir, "spade3")
    os.system(f"mkdir -p {spadedir}")
    for tid in df_spade["transcript_id"].unique():
        df_spade[df_spade["transcript_id"] == tid].to_csv(
            os.path.join(spadedir, f"{tid}.pfam"), index=None, header=None, sep="\t"
        )
    return df_spade


def mp_msa(appris_path: str, fasta_path: str, outdir: str, cpus: int):
    """Running a MSA (paired) in parallel

    Methods implemented: muscle
    If user want to run the pair alignment with another reference, change `run_msa_reference`

    Args:
        appris_path (str): APPRIS path.
        fasta_path (str): FASTA sequences path.
        outdir (str): Output directory to store the results.
        cpus (int): Number of cpus used to run in parallel this program.
    """
    df_appris_msa = annotation_reference(appris_path, fasta_path, outdir)
    gb = df_appris_msa.groupby("gene_id")
    df_group = [gb.get_group(x) for x in gb.groups]
    with mp.Pool(processes=cpus) as process:
        process.map(functools.partial(_run_msa_reference, outdir=outdir), df_group)


def qpfam_effects(df_reference: list, df_spade: list, pfam_effects_filepath: str) -> list:
    """Quantifying pfam effects and report transcript scores

    Args:
        df_appris_msa (list): pandas DataFrame from reference from the MSA alignment
        df_spade (list): pandas DataFrame from SPADE method.
        pfam_effects_filepath (str): File that summarizes Pfam effects over set
    of non-reference isoforms.

    Returns:
        list: final DataFrame with effects quantified.
    """
    df_pfam = load_pfam(pfam_effects_filepath)
    df_npfam_domains = (
        df_spade["transcript_id"]
        .value_counts()
        .reset_index()
        .rename(columns={"index": "transcript_id", "transcript_id": "n_pfam_domains"})
    )
    df = pd.merge(df_reference, df_npfam_domains, how="left", on="transcript_id")
    df.loc[df["pfam_effects_msa"].str.contains("Reference"), "n_pfam_domains_reference"] = df["n_pfam_domains"]
    df["n_pfam_domains_reference"] = df.groupby("gene_id")["n_pfam_domains_reference"].transform("max")
    df = pd.merge(df, df_pfam, how="left", on="transcript_id")
    df["pfam_score"] = df["pfam_score"].fillna(1)
    events_columns = df.select_dtypes(include="float").drop(["pfam_score"], axis=1).columns
    df[events_columns] = df[events_columns].fillna(0).astype(int)
    df.loc[df["Lost_residues_pfam"] < 10, "norm_Lost_residues_pfam"] = 0
    df = _get_pfam_percentages(df)
    df = df[
        [
            "gene_id",
            "transcript_id",
            "pfam_score",
            "pfam_domains_impact_score",
            "perc_Damaged_State",
            "perc_Lost_State",
            "Lost_residues_pfam",
            "Gain_residues_pfam",
            "pfam_effects_msa",
            "appris",
        ]
    ]
    return df


def _instance_counter(df: pd.DataFrame, column: str) -> pd.DataFrame:
    """Instance counter Function

    It counts the number of difference Instances in column Event type and State.

    Args:
        df (pd.DataFrame): input
        column (str): feature to count

    Returns:
        pd.DataFrame: transcript identifier + events per isoform
    """
    for instance in [x for x in df[column].unique() if str(x) != "nan"]:
        if column == "Event_type":
            df[f"{instance}_{column}"] = (
                (df[column] == instance) & ((df["State"].str.contains("Lost|Damaged")))
            ).astype(int)
        elif column == "State":
            df = df.drop_duplicates(["transcript_id", "Pfam_domain_number", "Pfam_domain_name", "State"])
            df[f"{instance}_{column}"] = df[column].str.contains(instance, na=False).astype(int)
    col_header_impact = list()
    for col_header in df.columns:
        if col_header.endswith(column):
            col_header_impact.append(col_header)
    df = df.groupby(["transcript_id"])[col_header_impact].sum().reset_index()
    return df


def _get_pfam_percentages(df: pd.DataFrame) -> pd.DataFrame:
    """Getting pfam scores from a processed Pfam transcript quantifier.

    Args:
        df (pd.DataFrame): pandas DataFrame with effects quantified over Pfam domains

    Returns:
        pd.DataFrame: pandas DataFrame with scores
    """
    df["pfam_domains_impact_score"] = 1 - (
        (df["Damaged_State"] + df["Lost_State"]) / df["n_pfam_domains_reference"]
    ).fillna(1)
    df["perc_Damaged_State"] = (df["Damaged_State"] / df["n_pfam_domains_reference"]).fillna(0)
    df["perc_Lost_State"] = (df["Lost_State"] / df["n_pfam_domains_reference"]).fillna(0)
    df["Intact_State"] = df["n_pfam_domains_reference"] - (df["Lost_State"] + df["Damaged_State"])
    df["perc_Intact_State"] = (df["Intact_State"] / df["n_pfam_domains_reference"]).fillna(0)
    df.loc[df["perc_Intact_State"] < 0, "perc_Intact_State"] = 0
    return df


def _run_muscle(infile: str, outfile: str) -> str:
    """Multiple Sequence ALigment

    'MUSCLE stands for MUltiple Sequence Comparison by Log- Expectation.
    MUSCLE is claimed to achieve both better average accuracy and better
    speed than ClustalW2 or T-Coffee, depending on the chosen options.'
    https://www.ebi.ac.uk/Tools/msa/muscle/

    Command line reference: https://drive5.com/muscle/manual/options.html

    Arguments:
        infile {str} -- Specified using the -in option. Must be in
            FASTA format. If any gaps are present in the input file,
            they will be discarded before making the alignment.
            If the -in option is not specified, or is specified
            as - (a minus sign), then input is taken from
            standard input.
        outfile {str} -- Specified using the -out option. By default,
            the output is in FASTA format with gaps added to align the
            sequences. If the -out option is not specified, or is
            specified as - (minus sign), then output is written to
            standard output. I use "afa" (aligned fasta) as the
            conventional filename extension.

    Returns:
        str -- Command line argument selected with CLUSTALW format
            as output and running in silent mode
    """
    cmd = f"muscle -in {infile} -out {outfile} -quiet -clw"
    return cmd


def _n_counter(df: pd.DataFrame, column: str) -> pd.DataFrame:
    """It counts the number of Events per Isoform.
    Note that n_events * n_pfam_domains = n_states

    Args:
        df (pd.DataFrame): Input.
        column (str): Feature to count.
    Returns:
        pd.DataFrame: Transcript identifier + events per isoform.
    """
    df["n_events"] = df[column].str.contains("^0$", na=False).astype(int)
    df = df.groupby(["transcript_id"])["n_events"].sum().reset_index()
    df["n_events"] = df["n_events"] + 1
    return df


def _run_msa_reference(df: list, outdir: str, msa_method: str = "muscle"):
    """Running a pair sequence aligment against a transcript reference per gene
    previously selected.

    Arguments:
        df {list} -- DataFrame with annotations and sequences to create the
            dictionary reference.
        outdir {str} -- Directory to store the results.

    Keyword Arguments:
        msa_method {str} -- Method to perform the MSA. At the moment, only
            muscle has been implemented (default: {'muscle'})
    """
    df["gene_id"] = df["gene_id"].astype(str)
    df_ref = df.loc[df["pfam_effects_msa"].str.contains("Reference")]
    if msa_method == "muscle":
        os.system(f"mkdir -p {os.path.join(outdir, 'muscle')}")
        os.system(f"mkdir -p {os.path.join(outdir, 'tmp_muscle')}")
    for gid, tid, seq in zip(df["gene_id"], df["transcript_id"], df["sequence"]):
        tids_ref = df_ref[df_ref["gene_id"] == gid]["transcript_id"].values
        seqs_ref = df_ref[df_ref["gene_id"] == gid]["sequence"].values
        if tid in tids_ref:
            pass
        elif tid not in tids_ref:
            for tid_ref, seq_ref in zip(tids_ref, seqs_ref):
                tid_ref_noext = tid_ref.split(".")[0]
                gid_noext = gid.split(".")[0]
                tid_noext = tid.split(".")[0]
                name = f"{gid_noext}.{tid_ref_noext}.{tid_noext}"
                infile = os.path.join(outdir, "tmp_muscle", name + ".faa")
                with open(infile, "w") as dualfasta_file:
                    dualfasta_file.write(f">{tid_ref_noext}\t{gid_noext}\n{seq_ref}\n>{tid_noext}\t{gid_noext}\n{seq}")
                os.system(_run_muscle(infile, os.path.join(outdir, "muscle", name + ".ali")))


def _select_reference_priority(df_reference: pd.DataFrame, sorted_list=[""]) -> pd.DataFrame:
    """Selecting the transcript reference by 'spade' score or 'appris' tag.
    Args:
        df_reference (pd.DataFrame): pandas DataFrame previously constructed
        feature (str, optional): Defaults to 'spade'.
    Returns:
        pd.DataFrame -- pandas DataFrame with reference selected
    """

    df_reference.loc[df_reference["flags"] == "protein_coding", "protein_coding"] = 1
    df_reference.loc[df_reference["flags"] != "protein_coding", "protein_coding"] = 0
    df_reference.loc[df_reference["ccdsid"].str.contains("CCDS"), "CCDS"] = 1
    df_reference.loc[~df_reference["ccdsid"].str.contains("CCDS"), "CCDS"] = 0
    df_reference.loc[df_reference["tsl"] == 1, "tsl_1"] = 1
    df_reference.loc[df_reference["tsl"] != 1, "tsl_1"] = 0
    df_reference["ccdsid_number"] = df_reference["ccdsid"].str.split("S").str[1].astype(float).fillna(0)
    df_reference.loc[df_reference["appris"].str.contains("PRINCIPAL"), "appris_order"] = 1
    df_reference.loc[~df_reference["appris"].str.contains("PRINCIPAL"), "appris_order"] = 0

    df = df_reference[
        df_reference.groupby(["gene_id"])["protein_coding"].transform(max) == df_reference["protein_coding"]
    ]
    df = df[df.groupby(["gene_id"])["spade"].transform(max) == df["spade"]]
    df = df[df.groupby(["gene_id"])["tsl_1"].transform(max) == df["tsl_1"]]
    df = df[df.groupby(["gene_id"])["CCDS"].transform(max) == df["CCDS"]]
    df = df[df.groupby(["gene_id"])["length"].transform(max) == df["length"]]
    df = df[df.groupby(["gene_id"])["ccdsid_number"].transform(min) == df["ccdsid_number"]]
    df = df[df.groupby(["gene_id"])["appris_order"].transform(min) == df["appris_order"]]

    df = df[(~df[["gene_id", "sequence"]].duplicated()) | (df["sequence"].isnull())]
    return df


def _unique_counter(df: pd.DataFrame, column: str, condition: str) -> pd.DataFrame:
    """It counts the number of residues Inserted or Lost/Swapped in an isoform.

    Args:
        df (pd.DataFrame): Input
        column (str): Feature to count
        condition (str): Gain or Lost depending if you want to account for Lost residues or Gain residues

    Returns:
        pd.DataFrame: Transcript identifier + condition selected unique counts
    """
    if condition == "Gain":
        df = (
            df.loc[df["Event_type"].str.contains("Insert", na=False)]
            .groupby("transcript_id")
            .agg({column: lambda x: x.unique().sum()})
            .reset_index()
        )
        df.columns = df.columns.str.replace("Lost", "Gain")
    elif condition == "Lost":
        df = (
            df.loc[~df["Event_type"].str.contains("Insert", na=False)]
            .groupby("transcript_id")
            .agg({column: lambda x: x.unique().sum()})
            .reset_index()
        )
    df = df.fillna(0)
    return df


if __name__ == "__main__":
    main()
