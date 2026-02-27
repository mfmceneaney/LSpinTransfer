#!/usr/bin/env python3

import argparse
import glob
import os
import sys
import pandas as pd
import numpy as np
from decimal import Decimal, ROUND_HALF_UP


def parse_args():
    parser = argparse.ArgumentParser(
        description="Format CSV columns using per-column minimum digit-place precision."
    )

    parser.add_argument(
        "paths",
        nargs="+",
        help="CSV file paths or glob patterns (e.g. data/*.csv)"
    )

    parser.add_argument(
        "--columns",
        nargs="+",
        required=True,
        help="Existing column names to format"
    )

    parser.add_argument(
        "--new-columns",
        nargs="+",
        required=False,
        help="New column names"
    )

    parser.add_argument(
        "--precision",
        nargs="+",
        required=True,
        help="Minimum digit place per column (e.g. 1e-2, 0.01, 1e1)"
    )

    parser.add_argument(
        "--reset-index",
        nargs=2,
        metavar=("COLUMN", "START"),
        help="Reset a column as a sequential index starting at START"
    )

    parser.add_argument(
        "--times-bounds",
        nargs=2,
        metavar=("LOW", "HIGH"),
        help="Set the LOW and HIGH bounds for using VALUE X10^{exp} LaTeX notation"
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite original CSV files"
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print out information used for reformatting numeric values",
    )

    return parser.parse_args()


def expand_paths(path_args):
    files = set()
    for p in path_args:
        matches = glob.glob(p)
        if matches:
            files.update(matches)
        elif os.path.isfile(p):
            files.add(p)
    return sorted(files)


def precision_exponent(precision):
    """Return exponent p such that precision = 10^p."""
    if precision>=10:
        return int(np.floor(np.log10(abs(precision))))
    d = Decimal(str(precision))
    return d.as_tuple().exponent


def round_to_precision(value, precision):
    d_value = Decimal(str(value))
    d_prec = Decimal(str(precision))
    return (d_value / d_prec).to_integral_value(rounding=ROUND_HALF_UP) * d_prec


def format_value(x, precision, verbose, lower_bound, upper_bound):
    if pd.isna(x):
        return x

    try:
        x = float(x)
    except Exception:
        return x

    precision = Decimal(str(precision))
    p = precision_exponent(precision)  # precision = 10^p

    if np.abs(x)<precision and x!=0.0:
        return format_value(0.0, precision, verbose, lower_bound, upper_bound)

    if x == 0:
        return f"{Decimal(0):.{max(0, -p)}f}"

    abs_x = abs(x)

    # Decide on scientific notation
    if abs_x < lower_bound or abs_x > upper_bound:
        exponent = max(int(np.floor(np.log10(abs_x))), p+1)

        # Normalize
        mantissa = x / (10 ** exponent)

        # Round mantissa to correct precision
        if verbose:
            print("INFO: x               = ",x)
            print("INFO: mantissa        = ",mantissa)
            print("INFO: precision       = ",precision)
            print("INFO: type(mantissa)  = ",type(mantissa))
            print("INFO: type(precision) = ",type(precision))
            print("INFO: exponent        = ",exponent)
            print("INFO: p               = ",p)
        mantissa = round_to_precision(mantissa, precision / Decimal(10 ** exponent) ) #, precision

        # REQUIRED mantissa decimal places
        mantissa_decimals = abs(p - (exponent)) #max(0, -(p - (exponent)))
        mantissa_str = f"{mantissa:.{mantissa_decimals}f}"
        return rf"{mantissa_str}\times 10^{{{exponent}}}"

    else:
        rounded = round_to_precision(x, precision)
        decimals = max(0, -p)
        return f"{rounded:.{decimals}f}"


def reset_index_column(df, column_name, start):
    start = int(start)
    df[column_name] = range(start, start + len(df))
    return df


def process_csv(csv_path, columns, new_columns, precisions, overwrite, reset_index_args, times_bounds_args, verbose):
    df = pd.read_csv(csv_path)

    # Set bounds for using times notation
    lower_bound, upper_bound = (float(el) for el in times_bounds_args) if times_bounds_args else (1e-3, 10)

    # Reset index column if requested
    if reset_index_args:
        idx_col, idx_start = reset_index_args
        df = reset_index_column(df, idx_col, idx_start)

    for col, prec in zip(columns, precisions):
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in {csv_path}")

        df[col] = df[col].apply(lambda x: format_value(x, prec, verbose, lower_bound, upper_bound))

    if new_columns: df = df.rename(columns=dict(zip(columns, new_columns)))

    output_path = csv_path if overwrite else os.path.join(
        os.getcwd(), os.path.basename(csv_path)
    )

    df.to_csv(output_path, index=False)
    print(f"Wrote: {output_path}")


def main():
    args = parse_args()

    if not args.new_columns or len(args.new_columns)==0:
        if not (len(args.columns) == len(args.precision)):
            sys.exit("Error: --columns and --precision must match in length")
    else:
        if not (len(args.columns) == len(args.new_columns) == len(args.precision)):
            sys.exit("Error: --columns, --new-columns, and --precision must match in length")

    csv_files = expand_paths(args.paths)
    if not csv_files:
        sys.exit("Error: No CSV files found")

    for csv in csv_files:
        if args.verbose:
            print("INFO: Loading csv: "+csv)
        process_csv(
            csv,
            args.columns,
            args.new_columns,
            args.precision,
            args.overwrite,
            args.reset_index,
            args.times_bounds,
            args.verbose
        )


if __name__ == "__main__":
    main()
