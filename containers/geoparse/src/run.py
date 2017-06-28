#! /usr/local/bin/python

import click

import typing
from typing import List
from downloader import Downloader
from functional import seq
from typing import *
from pprint import pprint

@click.command()
@click.option('--location', default="series", help='Number of greetings.')
@click.argument('samples', nargs=-1, required=True)
def download(location: str, samples):
    #for instance location series GSM1696283 GSM1696284
    return download_gsms(list(samples), "sra", location)


def download_gsms(gsms: List[str], sra_filetype: str = "sra", location: str = "./") -> typing.Dict[str, str]:
    from downloader import Downloader
    d = Downloader(location)
    #.filter(lambda kv: kv[0].endswith(".sra"))\
    files = seq(gsms)\
        .map(lambda gsm_id: d.download_gsm(gsm_id, sra_filetype))\
        .dict()
    pprint(files)
    return files

if __name__ == '__main__':
    download()

