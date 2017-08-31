#!/usr/bin/env python

"""
basic module to use as an example.
"""
from clip_analysis import plotting_methods as pm
import argparse

def main():
    """
    Main program.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--wd",
        required=True,
        type=basestring,
    )

    args = parser.parse_args()


if __name__ == "__main__":
    main()
