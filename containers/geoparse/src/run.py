#! /usr/local/bin/python
import os

import click

import typing
from typing import List
from downloader import Downloader
from functional import seq
from typing import *
from pprint import pprint
import pandas as pd


@click.command()
@click.option('--location', default="series", type=click.Path(), help='where to save downloaded files')
@click.option('--filetype', default="sra", type=click.Choice(['sra', 'fastq', 'fasta']), help='file type (if fastq then will try to extract with fastq-dump)')
@click.option('--keep_sra', default=False, type=bool, help='if we should keep sra after downloading')
@click.argument('samples', nargs=-1, required=True)
def download(location: str, filetype: str,  keep_sra: bool, samples):
    #for instance location series GSM1696283 GSM1696284
    fastq_dump_options = {
        'skip-technical': None,
        'clip': None,
        'split-files': None,
        'readids': None,
        'read-filter': 'pass',
        'dumpbase': None,
        'gzip': None
    }
    sra_kwargs = {
        "keep_sra": keep_sra,
        'filetype': filetype,
        "fastq_dump_options": fastq_dump_options
    }
    files = cast(typing.Dict[str, str], download_gsms(list(samples), sra_kwargs, location))
    frame = pd.DataFrame.from_dict(files, orient="index")
    p = os.path.abspath(os.path.join(location, "output.tsv"))
    frame.to_csv(p, sep="\t", header=False)
    return p


def download_gsms(gsms: List[str], sra_kwargs: Dict[str, str], location: str) -> typing.Dict[str, str]:
    from downloader import Downloader
    d = Downloader(location)
    files = seq(gsms)\
        .flat_map(lambda gsm_id: seq(d.download_gsm(gsm_id, sra_kwargs)))\
        .dict()
    return files

if __name__ == '__main__':
    download()

