#! /usr/local/bin/python

import click
from conditions import download

@click.command()
@click.argument('sra', nargs=-1, required=True)
def cli(sra: str):
    say_hello_to(sra)

if __name__ == '__main__':
    cli()

