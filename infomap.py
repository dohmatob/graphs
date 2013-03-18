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
                   is_connected=False,
                   n_timepoints=20,
                   time_resolution=0.005,
                   logstep=2,
                   do_plots=True):
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

    nx.write_graphml(G, filename + '.graphml')

    # compute connected components
    if not is_connected:
        connected_components = nx.connected_component_subgraphs(G)
    else:
        connected_components = [G]

    logtime = xrange(0, n_timepoints * logstep, logstep)

    clustering_matrix = np.zeros((len(G_nodes), n_timepoints), dtype='int')
    for k in xrange(len(connected_components)):
        print
        print "Zooming connected component number %i .." % k
        print
        H = connected_components[k]
        H_adj_mat = modularity.compute_adjacency_matrix(H)

        logtime_flow = modularity.compute_logtime_flow(H_adj_mat, n_timepoints,
                                                       logstep=logstep)
        for j in xrange(n_timepoints):
            t = logtime[j]
            print
            print '\tZooming log time step %i ' % t

            # compute flow matrix
            input_markov_flow_filename = filename + ('%i_component_%i_markov_'
                                                     'flow') % (j, k)
            flow_H = nx.DiGraph(logtime_flow[:, :, j])

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
            print commands.getoutput("infohiermap %s %s %s" % (
                    seed,
                    input_markov_flow_filename + '.net',
                    n_attempts))

            # parse output community structure
            infomaptree = parse_infomaptree(markov_flow_infomaptree_filename)
            assert len(infomaptree) == len(H_nodes)

            try:
                _colors = dict((G_nodes_map[x[2]], int(x[0].split(
                                ':')[0]) + clustering_matrix[:, j].max())
                               for x in infomaptree)
            except KeyError:
                _colors = dict((G_nodes_map[int(x[2])], int(x[0].split(
                                ':')[0]) + clustering_matrix[:, j].max())
                               for x in infomaptree)

            clustering_matrix[_colors.keys(), j] = _colors.values()

    # generate plots
    title = ("Markov Zooming: Optimal community structures for different "
             "timescales (log-spaced), CTMC time-sampled once every %.3fs "
             ) % time_resolution

    if do_plots:
        print "Generating plots .."
        display_markov_zoom(
            G, clustering_matrix, logtime, title=title)

    return clustering_matrix


def display_markov_zoom(G, clustering_matrix, time, title=None):
    """Function to display results of do_markov_zoom(..)

    Parameters
    ----------
    G: nx.Graph object
        Graph under consideration

    clustering_matrix: 2 x 2 np.ndarray of shape of ints (#nodes, #time points)
        clustering matrix for the network; each column contains the cluster
        membership of nodes to modules, in at each markov time gap (resolution)

    time: np.array
        Markov time gaps/resolutions at which zooming has been done

    title: string (optional, default None)
        supertitle for all generated plots/subplots

    """

    n, m = clustering_matrix.shape
    G_nodes = G.nodes()

    Q_list = []
    n_cols = 4
    pos = nx.graphviz_layout(G)
    n_rows = (len(time) - 1) / n_cols
    n_rows += 1  # because subplots are counter from 1, not 0
    n_rows += 1  # becase we'll add a modularity plot at the end
    fig = plt.figure()
    title = title if not title is None else (
        "Markov Zooming: Optimal community structures for "
        "different timescales (log-spaced) ")

    fig.suptitle(title, fontsize=14)
    for j in xrange(m):
        t = time[j]

        modules = list(set(clustering_matrix[:, j]))
        modules = dict((c, np.nonzero(clustering_matrix[:, j] == c)[0])
                       for c in modules)

        n_modules = len(modules)

        # compute modularity
        Q = modularity.compute_modularity(G, modules)
        Q_list.append(Q)

        # generate color for nodes in each module
        node_color = np.zeros(len(G_nodes))
        for color, c in zip(xrange(n_modules), modules):
            node_color[clustering_matrix[:, j] == c] = color

        # display colored graph proper
        ax = plt.subplot2grid((n_rows, n_cols),
                              np.unravel_index(j, (n_rows, n_cols)))
        nx.draw(G, pos, node_color=node_color,
                with_labels=False,
                node_size=30,
                )

        ax.axis('off')
        subplot_title = 'Log time step: %.3f\nCommunities: %i' % (
            t, n_modules)
        subplot_title += '\nModularity: %.3f' % Q
        ax.set_title(subplot_title, fontsize=10)

    # plot modularity curve over timescales
    ax = plt.subplot2grid((n_rows, n_cols), (n_rows - 1, 0), colspan=n_cols,
                         rowspan=2)

    ax.plot(time, Q_list)
    ax.set_title("Modularity curve", fontsize=12)
    ax.set_ylabel("modularity", fontsize=10)
    ax.set_xlabel("time span (logspace)", fontsize=10)

    # show plots
    plt.tight_layout(pad=0.5, w_pad=0, h_pad=0)
    plt.show()

if __name__ == '__main__':
    # sanitize cmd-line
    if len(sys.argv) < 2:
        print ("\r\n\tUsage: python %s <path_to_network_pajek_or_graphml_file>"
               ) % sys.argv[0]
        print ("\tExample: python infomap.py "
               "~/Downloads/univ_dataset_TSPE.graphml\r\n")
        sys.exit(1)

    if sys.argv[1].endswith(".net") or sys.argv[1].endswith(".paj")\
            or sys.argv[1].endswith(".NET"):
        G = nx.Graph(nx.Graph(nx.read_pajek(sys.argv[1])))
    elif sys.argv[1].endswith(".graphml"):
        G = nx.Graph(nx.Graph(nx.read_graphml(sys.argv[1])))
    elif sys.argv[1].endswith(".gml"):
        G = nx.Graph(nx.Graph(nx.read_gml(sys.argv[1])))

    # real business here
    do_markov_zoom(G, sys.argv[1], n_timepoints=16,
                   logstep=1,
                   )
