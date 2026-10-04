"""Audit the arithmetic stitching in the complete block-number spectrum."""

from compute_pbd_spectrum import realised_profiles, summary


def spectrum_conditions(m, n4, n5):
    alpha = (m - 1) % 3
    quintuple_floor = (
        (n5 == 0 or n5 >= 12)
        if alpha == 0
        else (5 * n5 >= (alpha * m))
    )
    quadruple_floor = (
        (n4 == 0 or n4 >= 13) if m % 4 == 1 else (4 * n4 >= m)
    )
    return quintuple_floor and quadruple_floor


checked = 0
for m in range(70, 501):
    if m % 4 not in (0, 1):
        continue
    pair_count = m * (m - 1) // 2
    alpha = (m - 1) % 3
    beta = (1 - m) % 4
    assert beta in (0, 1)

    for n5 in range(pair_count // 10 + 1):
        remainder = pair_count - 10 * n5
        if remainder < 0 or remainder % 6:
            continue
        n4 = remainder // 6
        if not spectrum_conditions(m, n4, n5):
            continue

        checked += 1
        n = n4 + n5
        assert n4 % 5 == pair_count % 5
        assert n5 % 3 == pair_count % 3
        assert 10 * n == pair_count + 4 * n4
        assert 6 * n == pair_count - 4 * n5

        # Exact integer forms of the two CHM2 middle-interval comparisons.
        below_middle = 10 * n < pair_count + beta * m + 680
        above_middle = 30 * n > 5 * pair_count - 4 * alpha * m - 1560
        in_middle = not below_middle and not above_middle
        assert sum((below_middle, in_middle, above_middle)) == 1

        if below_middle:
            assert 4 * n4 < beta * m + 680
            if m % 4 == 1:
                if n4 == 0:
                    assert m % 20 in (1, 5)
                else:
                    nu = next(t for t in range(13, 18) if t % 5 == pair_count % 5)
                    assert n4 >= nu
            else:
                mu = next(
                    t for t in range(m // 4, m // 4 + 5)
                    if t % 5 == pair_count % 5
                )
                assert n4 >= mu

        if above_middle:
            assert 5 * n5 < alpha * m + 390
            if alpha == 0:
                assert n5 == 0 or (n5 >= 12 and n5 % 3 == 0)
            elif alpha == 1:
                assert 5 * n5 >= m and n5 % 3 == 1
            else:
                assert 5 * n5 >= 2 * m and n5 % 3 == 0


assert checked > 100_000

# Check the exact cardinality formula and its residue-wise O(1) remainder.
cardinality_checks = 0
errors_by_residue = {r: [] for r in range(60)}
for m in range(13, 501):
    if m % 4 not in (0, 1):
        continue
    result = summary(m)
    assert result.realised_count == len(realised_profiles(m))
    cardinality_checks += 1

for m in range(60, 50_001):
    if m % 4 not in (0, 1):
        continue
    result = summary(m)
    errors_by_residue[m % 60].append(result.second_order_error)

for residue, errors in errors_by_residue.items():
    if not errors:
        continue
    # Beyond the tiny overlapping-endpoint cases, the exact error is periodic
    # once an admissible residue modulo 60 is fixed.
    assert len(set(errors)) == 1, (residue, set(errors))

print(f"Complete-spectrum stitching verified for {checked} parameter pairs.")
print(
    "Exact spectrum-cardinality formula verified for "
    f"{cardinality_checks} orders; the second-order remainder is constant "
    "on each admissible class modulo 60."
)
