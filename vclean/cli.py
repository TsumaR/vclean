import argparse
import sys

from vclean.modules import run_for_cli as run_vclean
from vclean.modules import download_database_for_cli as download_database

def cli():
    parser = argparse.ArgumentParser(
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"""vClean: a tool for assessing the contamination from environmental viral metagenomes

usage: vclean <program> [options]

program:
    run                  run full pipeline
    download_database    download the dependent database""",
    )

    # subparsers = parser.add_subparsers(dest='command')
    subparsers = parser.add_subparsers(help=argparse.SUPPRESS)

    # Command: run
    run_parser = subparsers.add_parser(
        'run',
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"""Run vClean

usage: vclean <input> <output> [options]""",
    )
    run_vclean.fetch_arguments(run_parser)

    # Command: download_database
    download_parser = subparsers.add_parser(
        'download_database',
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Download the latest version of vClean's databasen
        \nusage: vclean download_database <destination>""",
    )
    download_database.fetch_arguments(download_parser)
    
    # Handling no arguments case
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    elif len(sys.argv) == 2:
        # Execute the appropriate function based on the command
        if sys.argv[1] == 'run':
            run_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == 'download_database':
            download_parser.print_help()
            sys.exit(0)
        else:
            parser.print_help()
            sys.exit(0)

    args = vars(parser.parse_args())
    args["func"](args)

if __name__ == '__main__':
    cli()
