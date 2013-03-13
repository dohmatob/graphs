"""
:Author: dohmatob elvis dopgima

"""

import networkx as nx
import matplotlib.pyplot as plt
import re
import sys
import modularity
import commands
import numpy as np
import random

"""regexp pattern for a line in an infomap .tree file"""
INFOMAPTREE_LINE_PATTERN = ('\\n(?P<node_comm_hie>\d+?(?::\d+)*?) '
                            '(?P<pagerank>.+?) "(?P<node_id>.*?)"')


def parse_infomaptree(filename):
    """Functiopn parses a .tree file (as produced by say Rosval et
    al.'s infomap).

    Returns
    -------
    A list of triplets: last coord is node label, middle coord is pagerank,
    0th coord is node's hierachical community appartenance.

    """

    with open(filename, 'r') as fd:
        # read file
        dump = fd.read()  # XXX can't read very large files this way!!!

        # parse and return
        return re.findall(INFOMAPTREE_LINE_PATTERN, dump, re.DOTALL)


def do_markov_zoom(adj_mat, filename, time=None, n_attempts=10,
                compute_modularities=True, do_plots=True):
    """Runs the Michael Schaub et al. Markov Zooming algorithm.

    Parameters
    ----------
    adj_mat: np.ndarray (2 x 2)
        adjacency matrix of underlying network/graph

    filename: string
        prefix for .net, .clu, .tree, etc., output files

    time: iterable of object of floats (list, array, etc., optional;
    default None)
        series of Markov time gaps (resolutions) we're interested in

    n_attempts: int (optional, default 10)
        number of attempts made to partition the network (passed
        used in the Rosval et al.'s mapequation back-end

    compute_modularities: boolean (optional, default True)
        if true, the modularity of each parition (one per time point)  will
        be computed

    do_plots: boolean (optional, default True)
        if true, plots of the clustering will be generated

    Returns
    -------
    clustering_matrix: 2 x 2 np.ndarray of shape of ints (#nodes, #time points)
        clustering matrix for the network; each column contains the cluster
        membership of nodes to modules, in at each markov time gap (resolution)

    time: np.array
        Markov time gaps/resolutions at which zooming has been done

    modularity: np.array
        array of modularities of the partitionings described in clustering_mat

    """

    filename = filename.split(".")[0]

    G = nx.Graph(adj_mat)

    G_nodes = G.nodes()
    G_nodes_map = dict((G_nodes[j], j) for j in xrange(len(G_nodes)))
    if not isinstance(adj_mat, np.ndarray):
        adj_mat = modularity.compute_adjacency_matrix(G)

    nx.write_pajek(G, filename + '.net')

    # compute connected components
    connected_components = nx.connected_component_subgraphs(G)

    if time is None:
        time = np.logspace(-2, 2, num=9, base=10.)

    modularities = None
    if compute_modularities:
        modularities = []

    clustering_matrix = []
    for j in xrange(len(time)):
        t = time[j]
        print
        print "+" * 80
        print ('Zooming graph/network at Marvov time gap/resolution '
               'of %f ..') % t
        print "+" * 80

        # compute flow matrix
        colors = {}
        color_count = 0
        fufu = []
        for k in xrange(len(connected_components)):
            H = connected_components[k]
            H_adj_mat = modularity.compute_adjacency_matrix(H)
            input_markov_flow_filename = filename + ('%i_component_%i_markov_'
                                                     'flow') % (j, k)
            flow_H = modularity.compute_flow_graph(H_adj_mat, t)

            # derive a directed grapth representing the computed flow dynamics
            H_nodes = H.nodes()
            H_node_labels_map = dict((j, H_nodes[j])
                                     for j in xrange(len(H_nodes)))
            flow_H = nx.relabel_nodes(flow_H, H_node_labels_map)
            H_nodes = flow_H.nodes()
            nx.write_pajek(flow_H, input_markov_flow_filename + '.net')

            # run the Map Equation algorithm on the flow graph
            markov_flow_infomaptree_filename = input_markov_flow_filename
            markov_flow_infomaptree_filename += '.tree'
            seed = random.randint(0, 10000)
            commands.getoutput("infohiermap %s %s %s" % (
                    seed,
                    input_markov_flow_filename + '.net',
                    n_attempts))

            # parse output community structure
            infomaptree = parse_infomaptree(markov_flow_infomaptree_filename)
            assert len(infomaptree) == len(H_nodes)

            _colors = dict((G_nodes_map[x[2]], int(x[0].split(
                            ':')[0]) + color_count)
                          for x in infomaptree)
            assert not [c for c in _colors.values() if c in fufu]

            color_count += len(set(_colors.values()))
            fufu += list(set(_colors.values()))

            colors.update(_colors)

        clustering_matrix.append(colors.values())
        C = dict((c, [node
                      for node, com in colors.iteritems() if com == c])
                 for c in colors.values())

        # compute modularity
        if compute_modularities:
            Q = modularity.compute_modularity_from_adjacency_matrix(
                adj_mat,
                C,
                directed=G.is_directed())

            modularities.append(Q)

    # generate plots
    if do_plots:
        print "Generating plots .."
        display_markov_zoom(
            G, np.array(clustering_matrix).T, time, modularities)

    return clustering_matrix, time, modularities


def display_markov_zoom(G, clustering_matrix, time, modularities):
    n, m = clustering_matrix.shape

    assert len(time) == m
    assert len(modularities) == m

    G_nodes = G.nodes()

    n_cols = 4
    pos = nx.graphviz_layout(G)
    n_rows = (len(time) - 1) / n_cols
    n_rows += 1  # because subplots are counter from 1, not 0
    n_rows += 1  # becase we'll add a modularity plot at the end
    fig = plt.figure()
    fig.suptitle((".oO Markov Zooming: Optimal community structures for "
                  "different timescales (log-spaced) Oo."), fontsize=14)
    for j in xrange(m):
        t = time[j]
        Q = modularities[j]

        # generate color for nodes in each module
        node_color = np.zeros(len(G_nodes))
        modules = set(clustering_matrix[:, j])
        n_modules = len(modules)
        for color, c in zip(xrange(n_modules), modules):
            node_color[clustering_matrix[:, j] == c] = color

        # display colored graph proper
        ax = plt.subplot2grid((n_rows, n_cols),
                              np.unravel_index(j, (n_rows, n_cols)))
        nx.draw(G, pos, node_color=node_color, with_labels=False,
                node_size=30)

        ax.axis('off')
        subplot_title = 'Markov time gap: %.3f\nCommunities: %i' % (
            t, n_modules)
        if not modularities is None:
            subplot_title += '\nModularity: %.3f' % Q
        ax.set_title(subplot_title, fontsize=10)

    # plot modularity curve over timescales
    ax = plt.subplot2grid((n_rows, n_cols), (n_rows - 1, 0), colspan=n_cols,
                         rowspan=2)

    if not modularities is None:
        ax.plot(time, modularities)
        ax.set_title("Modularity curve", fontsize=12)
        ax.set_ylabel("modularity", fontsize=10)
        ax.set_xlabel("time span", fontsize=10)

    # show plots
    plt.tight_layout(pad=0.5, w_pad=0, h_pad=0)
    plt.show()

if __name__ == '__main__':
    G = nx.Graph(nx.Graph(nx.read_pajek(sys.argv[2])))
    time = np.logspace(-2, 2, num=12, base=10.)

    do_markov_zoom(G, sys.argv[2], time)
