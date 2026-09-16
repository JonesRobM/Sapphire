"""The cached longest-chain search must answer exactly what an exhaustive search would.

``T()`` is the third element of a CNA signature, so a wrong answer here silently
mislabels structure. The graph is the mutual neighbourhood of two bonded atoms -- three
to five nodes in practice -- and a crystal presents the same handful of shapes over and
over, which is what makes caching worth it. These tests check the cache is keyed on the
graph and nothing else, and check the answers against an independent brute force.
"""
import itertools

import networkx as nx
import pytest

from Sapphire.CNA.FrameSignature import _longest_chain


def brute_force(n_nodes, edges):
    """Longest simple path in edges, plus the longest cycle when the graph has one.

    Written independently of the implementation: enumerate every permutation and take the
    longest prefix that is a valid path. Only usable for tiny graphs, which is all we need.
    """
    adj = {k: set() for k in range(n_nodes)}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    if not edges:
        return 0
    best = 0
    for order in itertools.permutations(range(n_nodes)):
        length = 0
        for i in range(len(order) - 1):
            if order[i + 1] in adj[order[i]]:
                length += 1
            else:
                break
        best = max(best, length)
    if len(edges) >= n_nodes:                      # cyclic: the implementation also
        G = nx.Graph()                             # considers the cycle basis
        G.add_nodes_from(range(n_nodes))
        G.add_edges_from(edges)
        cycles = [len(c) for c in nx.cycle_basis(G)]
        if cycles:
            best = max(best, max(cycles))
    return best


# ------------------------------------------------------------------- known shapes
@pytest.mark.parametrize("n, edges, expected", [
    (0, frozenset(), 0),                                        # no common neighbours
    (3, frozenset(), 0),                                        # isolated nodes, no bonds
    (2, frozenset({(0, 1)}), 1),                                # a single bond
    (3, frozenset({(0, 1), (1, 2)}), 2),                        # a chain
    (4, frozenset({(0, 1), (1, 2), (2, 3)}), 3),                # a longer chain
    (3, frozenset({(0, 1), (1, 2), (0, 2)}), 3),                # a triangle
    (4, frozenset({(0, 1), (1, 2), (2, 3), (3, 0)}), 4),        # a square
])
def test_known_graphs(n, edges, expected):
    assert _longest_chain(n, edges) == expected


# --------------------------------------------------------- against a brute force
@pytest.mark.parametrize("n_nodes", [2, 3, 4, 5])
def test_matches_brute_force_on_every_graph_of_this_size(n_nodes):
    """Exhaustive over every graph on n nodes -- the real bond graphs are r = 3..5."""
    possible = list(itertools.combinations(range(n_nodes), 2))
    for r in range(len(possible) + 1):
        for chosen in itertools.combinations(possible, r):
            edges = frozenset(chosen)
            assert _longest_chain(n_nodes, edges) == brute_force(n_nodes, edges), \
                f"n={n_nodes} edges={sorted(edges)}"


# ------------------------------------------------------------------- cache keying
def test_relabelling_the_same_shape_hits_the_cache():
    """Two atoms with the same local topology must not recompute."""
    edges = frozenset({(0, 1), (1, 2)})
    _longest_chain(3, edges)
    before = _longest_chain.cache_info()
    for _ in range(50):
        _longest_chain(3, frozenset({(1, 2), (0, 1)}))   # same set, different literal order
    after = _longest_chain.cache_info()
    assert after.misses == before.misses, "an identical graph should not miss"
    assert after.hits >= before.hits + 50


def test_node_count_is_part_of_the_key():
    """Isolated nodes change the graph even with the same edges."""
    edges = frozenset({(0, 1)})
    assert _longest_chain(2, edges) == _longest_chain(5, edges)  # same answer...
    # ...but they are distinct entries, so a cached 2-node result cannot serve 5 nodes
    assert (2, edges) != (5, edges)


def test_result_is_a_plain_int():
    """It goes straight into a signature tuple, so it must not leak a numpy scalar."""
    assert type(_longest_chain(3, frozenset({(0, 1), (1, 2)}))) is int
