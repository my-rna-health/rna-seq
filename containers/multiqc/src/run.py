#! /usr/local/bin/python

import click
from hello import say_hello_to

@click.command()
@click.argument('sra', nargs=-1, required=True)
def cli(sra):
    say_hello_to(sra)

if __name__ == '__main__':
    cli()


