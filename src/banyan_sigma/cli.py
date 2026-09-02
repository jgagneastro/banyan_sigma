"""Command-line interface for batch BANYAN Sigma classifications."""

import argparse
import math
from pathlib import Path

import pandas as pd

from .core import membership_probability


REQUIRED_COLUMNS = ("RA", "DEC", "PMRA", "PMDEC", "EPMRA", "EPMDEC")
OPTIONAL_COLUMN_PAIRS = (("RV", "ERV"), ("PLX", "EPLX"), ("DIST", "EDIST"))
RESULT_COLUMNS = ("INPUT_ROW", "ASSOCIATION", "PROBABILITY")


def _probability(value):
    """Parse a probability constrained to the closed interval [0, 1]."""
    try:
        probability = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number between 0 and 1") from exc
    if not math.isfinite(probability) or not 0 <= probability <= 1:
        raise argparse.ArgumentTypeError("must be between 0 and 1")
    return probability


def _column_names(targets):
    """Map case-insensitive standard column names to the input CSV headers."""
    lookup = {}
    for column in targets.columns:
        normalized = str(column).upper()
        if normalized in lookup:
            raise ValueError(
                "Input column names must be unique when compared case-insensitively: "
                f"{lookup[normalized]!r} and {column!r}"
            )
        lookup[normalized] = column

    missing = [column for column in REQUIRED_COLUMNS if column not in lookup]
    if missing:
        raise ValueError("Input CSV is missing required columns: " + ", ".join(missing))

    column_names = {column: lookup[column] for column in REQUIRED_COLUMNS}
    if "NAME" in lookup:
        column_names["NAME"] = lookup["NAME"]

    for value_column, error_column in OPTIONAL_COLUMN_PAIRS:
        supplied = [column for column in (value_column, error_column) if column in lookup]
        if len(supplied) == 1:
            missing_column = error_column if supplied[0] == value_column else value_column
            raise ValueError(
                f"Input CSV column {supplied[0]} requires companion column {missing_column}"
            )
        if len(supplied) == 2:
            column_names[value_column] = lookup[value_column]
            column_names[error_column] = lookup[error_column]

    return column_names


def association_results(targets, output, min_probability=0.01):
    """Return long-form absolute probabilities above ``min_probability``.

    Field hypotheses are excluded because they are not associations. Results are
    ordered by input row and then by decreasing association probability.
    """
    conflicting = [column for column in RESULT_COLUMNS if column in targets.columns]
    if conflicting:
        raise ValueError(
            "Input CSV uses reserved result columns: " + ", ".join(conflicting)
        )

    probabilities = output["ALL"]
    association_columns = [
        column for column in probabilities.columns if "FIELD" not in str(column).upper()
    ]

    rows = []
    for row_position in range(len(targets)):
        row_probabilities = probabilities.iloc[row_position][association_columns]
        selected = row_probabilities[row_probabilities > min_probability].sort_values(
            ascending=False
        )
        target_values = targets.iloc[row_position].to_dict()
        for association, probability in selected.items():
            rows.append(
                {
                    "INPUT_ROW": row_position + 1,
                    **target_values,
                    "ASSOCIATION": association,
                    "PROBABILITY": float(probability),
                }
            )

    columns = ["INPUT_ROW", *targets.columns, "ASSOCIATION", "PROBABILITY"]
    return pd.DataFrame(rows, columns=columns)


def classify_csv(input_path, output_path, min_probability=0.01):
    """Classify targets from ``input_path`` and write filtered results."""
    targets = pd.read_csv(input_path)
    if targets.empty:
        raise ValueError("Input CSV contains no targets")

    column_names = _column_names(targets)
    output = membership_probability(
        stars_data=targets,
        column_names=column_names,
        use_rv="RV" in column_names,
        use_plx="PLX" in column_names,
        use_dist="DIST" in column_names,
    )
    results = association_results(targets, output, min_probability=min_probability)
    results.to_csv(output_path, index=False)
    return results


def build_parser():
    parser = argparse.ArgumentParser(
        prog="banyan-sigma",
        description=(
            "Classify targets from a CSV and write only non-field associations "
            "above an absolute-probability threshold."
        ),
    )
    parser.add_argument("input_csv", type=Path, help="CSV containing target observables")
    parser.add_argument("output_csv", type=Path, help="destination results CSV")
    parser.add_argument(
        "--min-probability",
        type=_probability,
        default=0.01,
        metavar="P",
        help="strict absolute-probability threshold in [0, 1] (default: 0.01)",
    )
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        results = classify_csv(
            args.input_csv,
            args.output_csv,
            min_probability=args.min_probability,
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        parser.error(str(exc))

    print(f"Wrote {len(results)} association result(s) to {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
