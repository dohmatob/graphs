"""
:Module: parsers
:Synopsis: Tools for parsing community files from MapEquation sofware
(infomap_dir, infohiermap_dir, etc.)
:Author: dohmatob elvis dopgima elvis[dot]dohmatob[at]inria[dot]fr

"""

import re

"""regexp pattern for a line in an infomap .tree file"""
INFOMAPTREE_LINE_PATTERN = ('\\n(?P<node_comm_hie>\d+?(?::\d+)*?) '
                            '(?P<pagerank>.+?) "(?P<node_id>.*?)"')


def parse_infomapclu(filename):
    with open(filename, 'r') as fd:
        clu = fd.read()

        comm = re.findall('\r\n(\d+)', clu)

        partitioning = [[j for j in xrange(len(comm)) if comm[j] == c]
                        for c in comm]

        return partitioning


def parse_infohiermaptree(filename):
    """Function parses a .tree file, produced by Rosval et al.'s
    infohiermap_dir).

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


def get_structure_onelevel(stuff):
    """Extracts one-level community structure

    Parameters
    ----------
    stuff: list
        tree structure, list of tuples; one per node

    Returns
    -------
    list of tuples representing collapsed tree structure, one
    tuple per node; last coordinate of tuple is node label

    """

    structure = []
    while stuff:
        x = stuff.pop()
        x_node = x[2]
        x_comm = ':'.join(x[0].split(':')[:-1])
        structure.append((x_comm, 0, [x_node]))

        _stuff = list(stuff)
        for y in _stuff:
            y_comm = ':'.join(y[0].split(':')[:-1])
            if y_comm == x_comm:
                y_node = y[2]
                structure[-1][2].append(y_node)
                stuff.remove(y)

    return structure


def is_terminal_node(x):
    """Determines whether a given node is terminal

    Parameters
    ----------
    x: tuple or list, or string
        node under inspection

    """

    if not isinstance(x, basestring):
        if isinstance(x[0], basestring):
            return True

    return False


def get_structure(stuff, simplified=True, only_terminal_nodes=False):
    """Extracts community structure

    Parameters
    ----------
    stuff: list
        tree structure, list of tuples; one per node

    simplified: bool (default True):
        if true, then output will be a list of list, each list being the
        node members of a cluster/community

    only_terminal_nodes: bool (default False)
        if True, community structure for only terminal nodes will be sought-for

    Returns
    -------
    list of tuples representing collapsed tree structure, one
    tuple per node; last coordinate of tuple is node label

    """

    # loop for community structure
    structure = []
    onelevel_structure = stuff
    while True:
        structure += onelevel_structure

        if len(onelevel_structure) < 2:
            break

        onelevel_structure = get_structure_onelevel(
            onelevel_structure)

    # simplify theze thingz
    if simplified:
        only_terminal_nodes = True

    # only terminal nodes, plz
    if only_terminal_nodes:
        structure = [x for x in structure if is_terminal_node(x[2])]

    if simplified:
        structure = [[int(node) for node in x[2]] for x in structure]

    return structure
