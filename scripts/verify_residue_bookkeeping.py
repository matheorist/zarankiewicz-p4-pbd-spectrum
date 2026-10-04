"""Verify the finite modulo-60 table in the residue-bookkeeping lemma."""


TABLE = {
    5: [0, 3, 36, 24, 12],
    8: [12, 0, 33, 21, 9],
    1: [24, 12, 0, 3, 36],
    4: [21, 9, 12, 0, 33],
    0: [5, 8, 1, 4, 17],
    9: [11, 4, 7, 20, 8],
}

SAME_PAIR = {
    5: ({1, 4}, {1, 4}),
    8: ({1, 4}, {1, 4}),
    1: ({5, 8}, {5, 8}),
    4: ({5, 8}, {5, 8}),
}
OPPOSITE_PAIRS = [({1, 4}, {5, 8}), ({5, 8}, {1, 4})]


checked = 0
for residue_mod_12, differences in TABLE.items():
    for j, difference in enumerate(differences):
        residue_mod_60 = residue_mod_12 + 12 * j
        assert (residue_mod_60 + difference) % 5 == 0
        q = (residue_mod_60 + difference) // 5
        w = q - difference

        if residue_mod_12 in SAME_PAIR:
            admissible_pairs = [SAME_PAIR[residue_mod_12]]
            assert difference <= 36
        else:
            admissible_pairs = OPPOSITE_PAIRS
            assert difference <= 20

        assert any(
            q % 12 in q_set and w % 12 in w_set
            for q_set, w_set in admissible_pairs
        )
        assert (4 * q + w - residue_mod_60) % 60 == 0
        checked += 1

assert checked == 30
print("Residue-bookkeeping table verified in all 30 admissible classes modulo 60.")
