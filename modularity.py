"""
:Module: modularity
:Author: d0hm4t06 3lv15 d0p91m4 <elvis[dot]dohmatob[at]inria[dot]fr>
:Synopsis: Routine functions for Markov zooming and modularity business

"""

import numpy as np
import scipy.sparse
import scipy.linalg
import scipy.stats
import networkx as nx


def wikipedia_example_graph():
    import itertools

    G = nx.Graph()
    for a in xrange(0, 9, 3):
        nodes = xrange(a, a + 3)
        G.add_edges_from(itertools.combinations(nodes, 2))
        G.add_edge(a, 9)

    return G


def rosval_mapequation_example_graph_1():
    G = nx.DiGraph()
    for a in xrange(0, 16, 4):
        nodes = xrange(a, a + 4)
        H = nx.DiGraph()
        H.add_nodes_from(nodes)
        H.add_edges_from([(x, ((x + 1) % 4 + a)) for x in nodes])
        G = nx.union(G, H)

    G.add_edges_from([(1, 5),
                      (2, 4),
                      (6, 10),
                      (7, 9),
                      (8, 14),
                      (11, 15),
                      (13, 3),
                      (12, 0)
                      ])

    return G


def compute_adjacency_matrix(G):
    """Computes the adjacency matrix of a graph or di-graph.

    Returns
    -------
    2D array

    """

    iG = nx.convert_node_labels_to_integers(G)
    adj_list = iG.adjacency_list()
    n_nodes = len(iG.nodes())

    adj_mat = np.zeros((n_nodes, n_nodes))
    for x in xrange(n_nodes):
        adj_mat[x, adj_list[x]] = 1

    return adj_mat


def compute_modularity_matrix_from_adjacency_matrix(A):
    in_degree = np.sum(A, axis=0)
    out_degree = np.sum(A, axis=1)
    m = np.sum(A)

    return A - 1. * np.outer(out_degree, in_degree) / m


# def compute_modularity_matrix(G):
#     return compute_modularity_matrix_from_adjacency_matrix(
#         compute_adjacency_matrix(G))


# def compute_steady_state_pi(adj_mat):
#     """Computes the steady-stead distribution of a random walk on an undirected
#     graph, Given its adjacency matrix.

#     Parameters
#     ----------
#     adj_mat: 2D array
#         The adjacency matrix of the underlying graph

#     Returns
#     -------
#     1D array

#     """

#     return 1. * np.sum(adj_mat, axis=0) / np.sum(adj_mat)  # d_j / 2|E|


def compute_transition_rates(adj_mat):
    """Computes the transition matrix (of flow rates) for the embedded Markov
    chain associeted with a random walk on a graph. The computed matrix embeds
    transition rates, not necessary
    probs!!!

    Parameters
    ----------
    adj_mat: 2D array
        The adjacency matrix of the underlying graph

    Returns
    -------
    2D array

    """

    # compute the node degrees
    degrees = np.sum(adj_mat, axis=0)
    degrees[degrees == 0] = adj_mat.shape[0]

    # compute D = diag(1 / degree_i), for each node i
    D = scipy.sparse.csr_matrix(np.diag(1. / degrees))

    # compute the transition matrix M = A\D
    M = D.dot(adj_mat)

    return M


def compute_logtime_flow(adj_mat, n_timepoints, time_resolution=0.005,
                         logstep=2,
                         threshold=1e-6,
                         ):
    """Computes a series of 2D arrays for the diffusion aasiociated with a
    random walk on a graph, at different (log-spaced) multiples of a given
    base time resolution. This corresponds exactly a log-spaced time sampled
    version of the CTMC of Schaub, the sampling frequency being 1 / resolution

    The sample CTMC is a DTMC with (1-step) transition matrix P given by:

        P_ij = 1 - q_ij x time_resoultion, if i != j
        P_ii = q_ii x time_resolution,

    where Q =(q_ij) is the (infinitesimal) generator matrix of the CTMC.

    Parameters
    ----------
    adj_mat: 2D array
        the adjacency matrix of the underlying graph

    n_timepoints: int (optional, default 2)
        number of timepoints to consider. For example, if this equals 5, then
        we'll only compute flows corresponding to the 1-, 2 ^ (logstep)-, ..,
        2 ^ (4 x logstep)-step transition matrices

    time_resolution: float (optional, default 0.01)
        1 / time_resolution corresponds to the sampling frequency of the CTMC

    logstep:
        linear spacing between timepoints in log-space. For example if this
        equals 1 (resp. 2), then the time steps will be 0, 2, 4, 8, 16, ..
        (0, 4, 16, 32, 64, .., resp.)

    threshold: float (optional, default 1e-6)
        threshold for elements in computed transition matrices; at each stage,
        all elements less than this threshold will be set to zero

    Returns
    -------
    3D array, corrsponding to the flow (one per time point) orderly packed
    together into an film

    """

    n = adj_mat.shape[0]

    # compute D = diag(degree_i), for each node i
    degrees = np.sum(adj_mat, axis=0)
    D = scipy.sparse.csr_matrix(np.diag(degrees))

    flow = np.zeros((n, n, n_timepoints))

    # compute transition rates
    trans_mat = compute_transition_rates(adj_mat)

    # compute generator matrix
    U = trans_mat * time_resolution
    P = np.eye(n) - np.diag(U.sum(axis=1)) + U

    # iteratively compute flow
    step = 2 ** logstep
    Q = scipy.sparse.csr_matrix(P)
    for k in xrange(n_timepoints):
        print Q
        print

        # compute flow ate time[k]
        flow[:, :, k] = D.dot(Q).toarray()

        # # threshold the flow
        # for i in xrange(n):
        #     threshold = scipy.stats.scoreatpercentile(flow[i, :, k], 7)
        #     flow[i, flow[i, :, k] < threshold, k] = 0.

        # update: Q = Q^(2^logstep)
        Q = scipy.sparse.csr_matrix(np.linalg.matrix_power(
                Q.toarray(),
                step))

    return flow


# def compute_flow_graph(flow_mat):
#     """Computes the (directed) flow graph associated with a random walk "
#     "on a graph under equilibrium.

#     Parameters
#     ----------
#     adj_mat: 2D array
#         the adjacency matrix of the underlying graph
#     t: float
#         upper limit of the time window [0, t) we're interested in.

#     Returns
#     -------
#     networkx.DiGraph object

#     """

#     # convert flow matrix to di-graph
#     flow_graph = nx.DiGraph(flow_mat)

#     return flow_graph


# def compute_dop_B(adj_mat, t):
#     pi = compute_steady_state_pi(adj_mat)

#     return  compute_flow_matrix(adj_mat, t) - np.dot(pi.T, pi)


def compute_modularity_from_adjacency_matrix(adj_mat, partition,
                                             directed=False):
    m = len(partition)
    n_nodes = adj_mat.shape[0]
    n_edges = np.sum(adj_mat)
    if directed:
        n_edges *= 2

    # S[c, n] is 1 if node n is in community c; 0 otherwise
    S = np.zeros((m, n_nodes))
    for j, c in zip(xrange(m), partition):
        S[j, partition[c]] = 1

    B = compute_modularity_matrix_from_adjacency_matrix(adj_mat)

    return np.trace(np.dot(S, np.dot(B, S.T))) / n_edges


def compute_modularity(G, partition):
    # total edge count (undirected edges double countedd)
    n_edges = len(G.edges())
    if not G.is_directed():
        n_edges *= 2

    adj_mat = compute_adjacency_matrix(G)

    return compute_modularity_from_adjacency_matrix(
        adj_mat,
        partition
        )


# def improved_2_split(B_g, s_g):
#     n_g = len(s_g)
#     s_g = np.array(s_g)
#     unmoved = list(xrange(n_g))

#     for i in xrange(n_g):
#         score = np.zeros(len(unmoved))
#         Q_0 = 0.5 * np.dot(s_g.T, np.dot(B_g, s_g))
#         for k in xrange(len(unmoved)):
#             s_g[unmoved[k]] *= -1
#             score[k] = 0.5 * np.dot(s_g.T, np.dot(B_g, s_g)) - Q_0
#             s_g[unmoved[k]] *= -1

#         j_ = np.argmax(score)
#         if score[j_] > 0:
#             s_g[j_] *= -1
#             unmoved.remove(unmoved[j_])
#         else:
#             break

#     return s_g


# def compute_B_g(B, g):
#     B_g = B[g, :][:, g]
#     B_g -= np.diag(B_g.sum(axis=0))

#     return B_g


# def split_module(B, g):
#     """Tries to split split a module g.

#     Parameters
#     ----------
#     B: numpy 2D array
#         modularity matrix of underlying graph/network
#     g: array-like
#         module to be splitted


#     Returns
#     -------
#     A piece of g (or g, if can't split)

#     """

#     # slice B with g x g
#     n_g = len(g)
#     B_g = compute_B_g(B, g)

#     # compute the spectrum (sorted by ascending eigen values) of B_g
#     eig_vals, eig_vecs = np.linalg.eig(B_g)
#     sort_perm = np.argsort(eig_vals)
#     eig_vals = eig_vals[sort_perm]
#     eig_vecs = eig_vecs[sort_perm]

#     # try to split g
#     if eig_vals[-1] > 0:
#         s_g = np.ones(n_g)
#         s_g[eig_vecs[-1] < 0] = -1
#         s_g = improved_2_split(B_g, s_g)
#         dQ = 0.5 * np.dot(s_g.T, np.dot(B_g, s_g))
#         if dQ > 0:
#             g_1 = np.array(g)[np.nonzero(s_g > 0)[0]]

#             # return chopped-off component of g
#             return g_1

#     # g is indivisible
#     return g


# def newman_split(G):
#     B = compute_modularity_matrix(G)
#     P = {"divisible": [G.nodes()], "indivisible": []}
#     k = 0
#     while P["divisible"]:
#         k += 1
#         g = P["divisible"][0]
#         g_1 = split_module(B, g)
#         if len(g_1) in [0, len(g)]:
#             P["divisible"].remove(g)
#             P["indivisible"].append(g)
#         else:
#             g_2 = [x for x in g if not x in g_1]
#             P["divisible"].remove(g)
#             for g in [g_1, g_2]:
#                 if len(g_1) == 1:
#                     P["indivisible"].append(g)
#                 else:
#                     P["divisible"].append(g)

#     return P

# # A = np.zeros((5, 5))
# # A[0,1] = A[2,3] = A[2,4] = A[3,4] = A[4,3] = A[4,2] = A[3,2] = A[1,0] = 1
# # B = compute_modularity_matrix_from_adjacency_matrix(A)
# # split_module(B, [0, 3, 4])
# H = rosval_mapequation_example_graph_1()
# B = compute_modularity_matrix(H)
# print newman_split(H)
# # g = [0, 1, 2, 9]
# # print improved_2_split(B, [-1, 1, 1, -1, -1, -1, -1, -1, -1, -1])
# # print compute_modularity(H, [[0, 1, 2], xrange(3, 10)])
# # print compute_modularity(H, [[3, 4, 5, 6, 7, 8], [0, 1, 2, 9]])
