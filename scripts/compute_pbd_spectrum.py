"""Compute the asymptotic PBD(m,{4,5}) block-number spectrum.

The main theorem in Zarankiewicz_p4_manuscript.tex is asymptotic in m.  This
script performs only its exact arithmetic: it lists the profiles predicted by
the theorem, counts them in two independent ways, and reports the endpoint
profiles excluded by the sharp local replication conditions.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


def first_at_least(lower: int, residue: int, modulus: int) -> int:
    return lower + (residue - lower) % modulus


def last_at_most(upper: int, residue: int, modulus: int) -> int:
    return upper - (upper - residue) % modulus


@dataclass(frozen=True)
class SpectrumSummary:
    order: int
    arithmetic_count: int
    realised_count: int
    excluded_count: int
    b_regular_lower: int
    b_regular_upper: int
    a_positive_floor: int
    steiner4_endpoint: bool
    steiner5_endpoint: bool
    second_order_error: Fraction


def spectrum_conditions(m: int, n4: int, n5: int) -> bool:
    if m % 4 not in (0, 1):
        return False
    alpha = (m - 1) % 3
    if alpha == 0:
        five_floor = n5 == 0 or n5 >= 12
    else:
        five_floor = 5 * n5 >= alpha * m
    four_floor = (
        n4 == 0 or n4 >= 13
        if m % 4 == 1
        else 4 * n4 >= m
    )
    return five_floor and four_floor


def arithmetic_profiles(m: int) -> list[tuple[int, int]]:
    """All nonnegative solutions of 6*n4+10*n5=binom(m,2)."""
    if m % 4 not in (0, 1):
        return []
    pairs = m * (m - 1) // 2
    profiles = []
    for n5 in range(pairs // 10 + 1):
        remainder = pairs - 10 * n5
        if remainder >= 0 and remainder % 6 == 0:
            profiles.append((remainder // 6, n5))
    return profiles


def realised_profiles(m: int) -> list[tuple[int, int]]:
    return [
        profile
        for profile in arithmetic_profiles(m)
        if spectrum_conditions(m, *profile)
    ]


def summary(m: int) -> SpectrumSummary:
    if m % 4 not in (0, 1):
        raise ValueError("the spectrum is empty unless m is 0 or 1 modulo 4")

    pairs = m * (m - 1) // 2
    alpha = (m - 1) % 3
    rho = pairs % 3
    sigma = pairs % 5

    if alpha == 0:
        b_lower = 12
    else:
        b_lower = first_at_least((alpha * m + 4) // 5, rho, 3)

    if m % 4 == 1:
        a_floor = first_at_least(13, sigma, 5)
    else:
        a_floor = first_at_least(m // 4, sigma, 5)

    b_upper_raw = (pairs - 6 * a_floor) // 10
    b_upper = last_at_most(b_upper_raw, rho, 3)
    interval_count = max(0, (b_upper - b_lower) // 3 + 1)

    steiner4 = alpha == 0
    steiner5 = m % 20 in (1, 5)
    realised_count = interval_count + int(steiner4) + int(steiner5)

    b_arithmetic_upper = last_at_most(pairs // 10, rho, 3)
    arithmetic_count = (b_arithmetic_upper - rho) // 3 + 1

    # On each fixed admissible residue class the exact count differs by O(1)
    # from this quadratic polynomial.
    linear_coefficient = (
        Fraction(1, 60)
        + Fraction(alpha, 15)
        + (Fraction(1, 20) if m % 4 == 0 else 0)
    )
    second_order_model = Fraction(m * m, 60) - linear_coefficient * m

    result = SpectrumSummary(
        order=m,
        arithmetic_count=arithmetic_count,
        realised_count=realised_count,
        excluded_count=arithmetic_count - realised_count,
        b_regular_lower=b_lower,
        b_regular_upper=b_upper,
        a_positive_floor=a_floor,
        steiner4_endpoint=steiner4,
        steiner5_endpoint=steiner5,
        second_order_error=Fraction(realised_count) - second_order_model,
    )

    return result


def write_csv(path: Path, rows: list[SpectrumSummary]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(SpectrumSummary.__dataclass_fields__)
        for row in rows:
            writer.writerow(
                [
                    row.order,
                    row.arithmetic_count,
                    row.realised_count,
                    row.excluded_count,
                    row.b_regular_lower,
                    row.b_regular_upper,
                    row.a_positive_floor,
                    int(row.steiner4_endpoint),
                    int(row.steiner5_endpoint),
                    str(row.second_order_error),
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("orders", nargs="*", type=int, default=[3601, 3604, 3605, 3608])
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--list", dest="list_order", type=int)
    parser.add_argument(
        "--residue-table",
        action="store_true",
        help="print the exact bounded term in (7) for each class modulo 60",
    )
    args = parser.parse_args()

    rows = [summary(m) for m in args.orders]
    for row in rows:
        print(
            f"m={row.order}: realised={row.realised_count}/"
            f"{row.arithmetic_count}, excluded={row.excluded_count}, "
            f"regular b-range=[{row.b_regular_lower},{row.b_regular_upper}], "
            f"positive a-floor={row.a_positive_floor}, "
            f"second-order error={row.second_order_error}"
        )

    if args.list_order is not None:
        print(f"\nTheorem-predicted profiles for m={args.list_order}:")
        for n4, n5 in realised_profiles(args.list_order):
            blocks = n4 + n5
            total_degree = 4 * n4 + 5 * n5
            print(f"n4={n4:8d} n5={n5:8d} n={blocks:8d} Z={total_degree:8d}")

    if args.csv:
        write_csv(args.csv, rows)

    if args.residue_table:
        print("\nresidue alpha eta linear_coefficient bounded_term")
        for residue in range(60):
            if residue % 4 not in (0, 1):
                continue
            result = summary(600 + residue)
            alpha = (residue - 1) % 3
            eta = int(residue % 4 == 0)
            coefficient = (
                Fraction(1, 60)
                + Fraction(alpha, 15)
                + Fraction(eta, 20)
            )
            print(
                f"{residue:2d} {alpha:5d} {eta:3d} "
                f"{str(coefficient):>18s} {str(result.second_order_error):>12s}"
            )


if __name__ == "__main__":
    main()
