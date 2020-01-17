# Copyright (c)  2020, BooleanTrans (previously CardioTrans)
# Author: Rohit Suratekar
# Email: rsuratekar [at] iimcb.gov.pl
# URL: https://github.com/rohitsuratekar/CardioTRaNs
# Organization: Winata Lab, IIMCB, Warsaw
# Description: This project aims to generate Transcription Regulatory
# Network of the developing zebrafish heart. Analysis is focused on specific
# question, however, all the functions are robust enough to handle any kind
# of similar dataset
#
# Main configuration file for the project


import argparse

from helper import ConfigParser


def test():
    from analysis.deseq2 import run

    run(config)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Argument parsing for the "
                                                 "BooleanTrans project")

    parser.add_argument("-c",
                        "--config",
                        type=str,
                        help="Full path of the configuration file.")

    # Parse the arguments
    args = parser.parse_args()
    # Sanity check
    if args.config is None:
        raise FileNotFoundError("Configuration file not provided. Please "
                                "provide valid config file to start the "
                                "analysis. Example of such file is present "
                                "at the base of this package. \nUsage: "
                                "main.py -c /full/path/config.json")
    # Generate the configuration object
    config = ConfigParser(args.config)

    # Test function
    test()
