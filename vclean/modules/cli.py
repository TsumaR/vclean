import argparse
import sys
import vClean

def cli():
    parser = argparse.ArgumentParser(
        usage='vClean <command> [options]',
        description='vClean: a tool for assessing the quality of viral metagenomes',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    subparsers = parser.add_subparsers(dest='command')

    # Command: run
    run_parser = subparsers.add_parser(
        'run',
        help='run full pipeline for viral genome analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""usage: run.py [-h] [-i INPUT] [-d CHECKV_DIR] [-tmp TMP] [-t THREADS] [-p PROTEIN] [-n GENE] [-o OUTPUT] [--translate_table T_TABLE] [-m MODE]
              [--skip_feature_table SKIP_FEATURE_TABLE] [--skip_lgb_step SKIP_LGB_STEP] [--f_table F_TABLE] [-pr THRESHOLD]"""
    )
    vClean.run.fetch_arguments(run_parser)

    # Command: download_database
    download_parser = subparsers.add_parser(
        'download_database',
        help='download the latest version of vClean database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""usage: vClean download_database [-d DESTINATION]"""
    )
    vClean.download_database.fetch_arguments(download_parser)

    # Handling no arguments case
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # Parse and execute the appropriate function
    args = parser.parse_args()

    # Execute the appropriate function based on the command
    if args.command == 'run':
        vClean.run.execute(args)
    elif args.command == 'download_database':
        vClean.download_database.execute(args)
    else:
        parser.print_help()
        sys.exit(0)

if __name__ == "__main__":
    cli()

