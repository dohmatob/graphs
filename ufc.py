"""
:Author: dohmatob elvis dopgima elvis[dot]dohmatob[at]inria[dot]fr
:Synopsis: A new community detection algorithm using "Markov diffusion"
:Module: ufc (Ultimate Fighter Champion)

"""

import sys
import scipy.linalg
import pylab as pl
import numpy as np
import modularity
import networkx as nx


def compute_modularity(H, pi, expmtQ):
    pi = pi.reshape((1, len(pi)))

    R = np.dot(H,
               np.dot(np.dot(np.diag(pi[0, :]), expmtQ) - np.dot(pi.T, pi),
                      H.T)
               )

    return np.trace(R)


def do_dohmatob(G, times=None):
    # prepare timepoints
    if times is None:
        times = np.arange(0, 1000, .01)

    # sanitize graph
    # G = nx.Graph(G)
    G = nx.convert_node_labels_to_integers(G)
    n_nodes = len(G.nodes())

    # compute adjacency matrix
    adj_mat = modularity.compute_adjacency_matrix(G)

    # compute infinitesimal generator of associated CTMC
    generator_mat = modularity.compute_infinitesimal_generator(adj_mat)
    generator_mat = -1. * modularity.compute_laplacian(adj_mat,
                                                       normalized=True)

    # generator_mat = np.abs(generator_mat)

    # # # compute steady-state distribution
    steady_state_pi = modularity.compute_steady_state_pi(adj_mat)

    # # XXX rm the following code stub
    # H = np.zeros((n_nodes, n_nodes))
    # H[:, 0] = 1
    # expmtQ = scipy.linalg.expm(1000000 * generator_mat.T)
    # print compute_modularity(H, steady_state_pi, expmtQ)

    # return
    # run [my] perfectly asynchronous LPA
    Q_vec = np.zeros(len(times))
    clustering_mat = np.zeros((n_nodes, len(times)))
    paint_0  = np.eye(n_nodes)

    # paint_0 = np.zeros((n_nodes, n_nodes))
    # paint_0[0, [0, 1, 2]] = 1
    # paint_0[1, [3, 4, 5]] = 1
    # paint_0[2, [6, 7, 8, 9]] = 1

    # expmtQ = scipy.linalg.expm(1000 * generator_mat)
    # paint = np.dot(paint_0, expmtQ)
    # print compute_modularity(paint,
    #                          steady_state_pi,
    #                          expmtQ)
    # return

    # def f(j):
    #     t = times[int(j)]

    #     expmtQ = scipy.linalg.expm(t * generator_mat)
    #     paint = np.dot(paint_0, expmtQ)

    #     print "time: %.3f" % t

    #     Q = compute_modularity(paint, steady_state_pi,
    #                            expmtQ)

    #     # C = np.argmax(paint, axis=0)
    #     # partition = dict((c, [n for n in xrange(n_nodes) if C[n] == c])
    #     #                  for c in set(C))

    #     # # Q = compute_modularity(
    #     # #     paint,
    #     # #     steady_state_pi,
    #     # #     expmtQ)

    #     # B = modularity.compute_modularity_matrix_from_adjacency_matrix(
    #     #     adj_mat)

    #     # Q = np.trace(np.dot(paint, np.dot(B, paint.T)))
    #     # # Q = modularity.compute_modularity_from_adjacency_matrix(
    #     # #     adj_mat,
    #     # #     partition,
    #     # #     )

    #     return -Q

    # def fprime(j):
    #     t = times[int(j)]

    #     expmtQ = scipy.linalg.expm(t * generator_mat)
    #     paint = np.dot(paint_0, expmtQ)

    #     print "time: %.3f" % t

    #     B = modularity.compute_modularity_matrix_from_adjacency_matrix(
    #         adj_mat)
    #     B = np.dot(B, generator_mat.T) + np.dot(generator_mat,
    #                                             B)
    #     gradQ = np.trace(np.dot(paint, np.dot(B, paint.T)))

    #     return -gradQ

    # import scipy.optimize
    # guess = [0]
    # bounds = [(0, 100)]
    # print scipy.optimize.fmin_cg(f, guess, fprime=None,
    #                              gtol=1e-10,
    #                              full_output=True)

    # return

    for j in xrange(len(times)):
        t = times[j]

        expmtQ = scipy.linalg.expm(
            t * generator_mat)
        paint = np.dot(paint_0, expmtQ)
        paint_rs = paint.sum(axis=1)
        paint_ = paint / paint_rs[:, np.newaxis]
        paint = paint_

        print "time: %.3f" % t
        C = np.argmax(paint, axis=0)

        partition = dict((c, [n for n in xrange(n_nodes) if C[n] == c])
                         for c in set(C))

        # Q = compute_modularity(
        #     paint,
        #     steady_state_pi,
        #     expmtQ)

        # Q = compute_modularity(paint,
        #                        steady_state_pi,
        #                        expmtQ,
        #                        )

        Q = modularity.compute_modularity_from_adjacency_matrix(
            adj_mat,
            partition)

        Q_vec[j] = Q
        clustering_mat[:, j] = C

    # QA
    perm = np.argsort(Q_vec)
    best_t = times[perm][-1]
    best_Q = Q_vec[perm][-1]

    # for j in xrange(len(times)):
    #     C = clustering_mat[:, perm][:, j - len(times)]
    #     if len(set(C)) == 5:
    #         best_C = C
    #         print "BEST:", set(best_C)
    #         break

    best_C = clustering_mat[:, perm][:, -1]

    print times[perm]
    print Q_vec[perm]
    print clustering_mat[:, perm]
    print "best time: %.3f, best Q: %.3f, best C: %s" % (
        best_t, best_Q, best_C)
    nx.draw_graphviz(G, node_color=best_C,
                     # node_size=30,
                     with_labels=False,
                     )

    pl.figure()
    pl.plot(times, Q_vec)
    pl.title("Modularity (Girvan-Newman) over 'Markov time'")
    pl.xlabel("time")
    pl.ylabel("modularity")
    pl.show()

if __name__ == '__main__':
    # sanitize cmd-line
    if len(sys.argv) < 2:
        print ("\r\n\tUsage: python %s <path_to_network_pajek_or_graphml_file>"
               ) % sys.argv[0]
        print ("\tExample: python %s "
               "example_data/flow.net\r\n" % sys.argv[0])
        sys.exit(1)

    if sys.argv[1].endswith(".net") or sys.argv[1].endswith(".paj")\
            or sys.argv[1].endswith(".NET"):
        G = nx.read_pajek(sys.argv[1])
    elif sys.argv[1].endswith(".graphml"):
        G = nx.read_graphml(sys.argv[1])
    elif sys.argv[1].endswith(".gml"):
        G = nx.read_gml(sys.argv[1])

    if not ''.join(sys.argv[1:]).endswith('--directed'):
        G = nx.Graph(G)

    # G = nx.DiGraph(modularity.rosval_mapequation_example_graph_1())
    do_dohmatob(G,
                times=np.arange(0, 1000, .1))
