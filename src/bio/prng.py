"""Deterministic random numbers shared with the browser version (web/recipe.js), so both produce identical shuffles.

mulberry32 seeded by FNV-1a(32-bit) of a task key (e.g. "patient|field|T|neut"), so results don't depend on
processing order. choose() is a partial Fisher-Yates over the pool, consuming one draw per pick.
"""
import numpy as np

M = 0xFFFFFFFF


def fnv1a(s):
    h = 0x811C9DC5
    for b in s.encode("utf-8"):
        h = ((h ^ b) * 0x01000193) & M
    return h


class Mulberry32:
    def __init__(self, seed):
        self.a = seed & M

    def random(self):
        self.a = (self.a + 0x6D2B79F5) & M
        a = self.a
        t = ((a ^ (a >> 15)) * (1 | a)) & M
        t = ((t + (((t ^ (t >> 7)) * (61 | t)) & M)) & M) ^ t
        return ((t ^ (t >> 14)) & M) / 4294967296.0


def choose(rng, pool, k):
    """k distinct items from pool (1-D int array), partial Fisher-Yates."""
    p = np.array(pool, copy=True)
    n = len(p)
    for i in range(k):
        j = i + int(rng.random() * (n - i))
        p[i], p[j] = p[j], p[i]
    return p[:k]
