"""Prefix-free codes for the census grammars.

Every choice in a program is written with a prefix-free code, so a program's length in bits is well defined and the set
of programs of each exact length can be enumerated.

  truncated binary  - a uniform choice among n options (n=1 costs 0 bits)
  hierarchical      - a parameter drawn from coarse-to-fine levels: level l is announced by (l-1) ones and a closing
                      zero (no closing zero on the last level), then the value is a truncated-binary index inside
                      the level. Coarse values are therefore cheaper than finely tuned ones.
  continuation bit  - '1' = another rule follows, '0' = end of program.
"""


def tb_len(i, n):
    """bits used by truncated-binary code of index i among n options."""
    if n <= 1: return 0
    k = n.bit_length() - 1
    u = (1 << (k + 1)) - n
    return k if i < u else k + 1


def tb_encode(i, n):
    assert 0 <= i < n, (i, n)
    if n <= 1: return ""
    k = n.bit_length() - 1
    u = (1 << (k + 1)) - n
    return format(i, f"0{k}b") if (i < u and k > 0) else ("" if i < u else format(i + u, f"0{k + 1}b"))


def tb_decode(bits, pos, n):
    """returns (index, new_pos)."""
    if n <= 1: return 0, pos
    k = n.bit_length() - 1
    u = (1 << (k + 1)) - n
    v = int(bits[pos:pos + k], 2) if k > 0 else 0
    if v < u: return v, pos + k
    v = int(bits[pos:pos + k + 1], 2)
    return v - u, pos + k + 1


def tb_all(n):
    """[(index, codeword)] for every option — used by the enumerator."""
    return [(i, tb_encode(i, n)) for i in range(n)]


def hier_codes(levels):
    """levels = [[values of level 1], [values of level 2], ...] -> [(value, codeword)] in coarse-to-fine order."""
    out, L = [], len(levels)
    for li, vals in enumerate(levels):
        prefix = "1" * li + ("0" if li < L - 1 else "")
        for i, v in enumerate(vals):
            out.append((v, prefix + tb_encode(i, len(vals))))
    return out


def hier_decode(bits, pos, levels):
    L, li = len(levels), 0
    while li < L - 1 and bits[pos] == "1": li += 1; pos += 1
    if li < L - 1: pos += 1                                   # closing zero
    i, pos = tb_decode(bits, pos, len(levels[li]))
    return levels[li][i], pos
