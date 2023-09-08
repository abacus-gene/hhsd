'''
THIS MODULE CONTAINS CODE TAKEN FROM THE ORIGINAL ETE3 REPOSITORY:

https://github.com/etetoolkit/ete/blob/ete3/ete3/parser/newick.py
https://github.com/etetoolkit/ete/blob/ete3/ete3/coretype/tree.py

THE CODE WAS MODIFIED TO REMOVE DEPENDENCIES ON PYQT5, SIX, AND ETE3 ITSELF.
ONLY FUNCTIONS AND METHODS RELEVANT TO HHMS WERE KEPT, OTHERS WERE REMOVED.
'''

# #START_LICENSE###########################################################
#
#
# This file is part of the Environment for Tree Exploration program
# (ETE).  http://etetoolkit.org
#
# ETE is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ETE is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ETE.  If not, see <http://www.gnu.org/licenses/>.
#
#
#                     ABOUT THE ETE PACKAGE
#                     =====================
#
# ETE is distributed under the GPL copyleft license (2008-2015).
#
# If you make use of ETE in published work, please cite:
#
# Jaime Huerta-Cepas, Joaquin Dopazo and Toni Gabaldon.
# ETE: a python Environment for Tree Exploration. Jaime BMC
# Bioinformatics 2010,:24doi:10.1186/1471-2105-11-24
#
# Note that extra references to the specific methods implemented in
# the toolkit may be available in the documentation.
#
# More info at http://etetoolkit.org. Contact: huerta@embl.de
#
#
# #END_LICENSE#############################################################


import re
import os
import copy
from collections import deque
from functools import cmp_to_key

__all__ = ["read_newick", "write_newick", "print_supported_formats", "Tree", "TreeNode"]


'''
NEWICK RELATED FUNCTIONS FROM ETE3
'''


ITERABLE_TYPES = set([list, set, tuple, frozenset])

# Regular expressions used for reading newick format
_ILEGAL_NEWICK_CHARS = ":;(),\[\]\t\n\r="
_NON_PRINTABLE_CHARS_RE = "[\x00-\x1f]+"

_NHX_RE = "\[&&NHX:[^\]]*\]"
_FLOAT_RE = "\s*[+-]?\d+\.?\d*(?:[eE][-+]?\d+)?\s*"
#_FLOAT_RE = "[+-]?\d+\.?\d*"
#_NAME_RE = "[^():,;\[\]]+"
_NAME_RE = "[^():,;]+?"

# thanks to: http://stackoverflow.com/a/29452781/1006828
_QUOTED_TEXT_RE = r"""((?=["'])(?:"[^"\\]*(?:\\[\s\S][^"\\]*)*"|'[^'\\]*(?:\\[\s\S][^'\\]*)*'))"""
#_QUOTED_TEXT_RE = r"""["'](?:(?<=")[^"\\]*(?s:\\.[^"\\]*)*"|(?<=')[^'\\]*(?s:\\.[^'\\]*)*')""]"]"""
#_QUOTED_TEXT_RE = r"""(?=["'])(?:"[^"\\]*(?:\\[\s\S][^"\\]*)*"|'[^'\\]*(?:\\[\s\S][^'\\]*)*')]"]")"]"""

_QUOTED_TEXT_PREFIX='ete3_quotref_'

DEFAULT_DIST = 1.0
DEFAULT_NAME = ''
DEFAULT_SUPPORT = 1.0
FLOAT_FORMATTER = "%0.6g"
#DIST_FORMATTER = ":"+FLOAT_FORMATTER
NAME_FORMATTER = "%s"

def set_float_format(formatter):
    ''' Set the conversion format used to represent float distances and support
    values in the newick representation of trees.

    For example, use set_float_format('%0.32f') to specify 32 decimal numbers
    when exporting node distances and bootstrap values.

    Scientific notation (%e) or any other custom format is allowed. The
    formatter string should not contain any character that may break newick
    structure (i.e.: ":;,()")

    '''
    global FLOAT_FORMATTER
    FLOAT_FORMATTER = formatter
    #DIST_FORMATTER = ":"+FLOAT_FORMATTER

# Allowed formats. This table is used to read and write newick using
# different convenctions. You can also add your own formats in an easy way.
#
#
# FORMAT: [[LeafAttr1, LeafAttr1Type, Strict?], [LeafAttr2, LeafAttr2Type, Strict?],\
#    [InternalAttr1, InternalAttr1Type, Strict?], [InternalAttr2, InternalAttr2Type, Strict?]]
#
# Attributes are placed in the newick as follows:
#
# .... ,LeafAttr1:LeafAttr2)InternalAttr1:InternalAttr2 ...
#
#
#           /-A
# -NoName--|
#          |          /-B
#           \C-------|
#                    |          /-D
#                     \E-------|
#                               \-G
#
# Format 0 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)1.000000:0.642905)1.000000:0.567737);
# Format 1 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)E:0.642905)C:0.567737);
# Format 2 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)1.000000:0.642905)1.000000:0.567737);
# Format 3 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)E:0.642905)C:0.567737);
# Format 4 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)));
# Format 5 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729):0.642905):0.567737);
# Format 6 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)E)C);
# Format 7 = (A,(B,(D,G)E)C);
# Format 8 = (A,(B,(D,G)));
# Format 9 = (,(,(,)));

NW_FORMAT = {
  0: [['name', str, True],  ["dist", float, True],    ['support', float, True],   ["dist", float, True]], # Flexible with support
  1: [['name', str, True],  ["dist", float, True],    ['name', str, True],      ["dist", float, True]], # Flexible with internal node names
  2: [['name', str, False], ["dist", float, False],   ['support', float, False],  ["dist", float, False]],# Strict with support values
  3: [['name', str, False], ["dist", float, False],   ['name', str, False],     ["dist", float, False]], # Strict with internal node names
  4: [['name', str, False], ["dist", float, False],   [None, None, False],        [None, None, False]],
  5: [['name', str, False], ["dist", float, False],   [None, None, False],        ["dist", float, False]],
  6: [['name', str, False], [None, None, False],      [None, None, False],        ["dist", float, False]],
  7: [['name', str, False], ["dist", float, False],   ["name", str, False],       [None, None, False]],
  8: [['name', str, False], [None, None, False],      ["name", str, False],       [None, None, False]],
  9: [['name', str, False], [None, None, False],      [None, None, False],        [None, None, False]], # Only topology with node names
  100: [[None, None, False],  [None, None, False],      [None, None, False],        [None, None, False]] # Only Topology
}


def format_node(node, node_type, format, dist_formatter=None,
                support_formatter=None, name_formatter=None,
                quoted_names=False):
    if dist_formatter is None: dist_formatter = FLOAT_FORMATTER
    if support_formatter is None: support_formatter = FLOAT_FORMATTER
    if name_formatter is None: name_formatter = NAME_FORMATTER

    if node_type == "leaf":
        container1 = NW_FORMAT[format][0][0] # name
        container2 = NW_FORMAT[format][1][0] # dists
        converterFn1 = NW_FORMAT[format][0][1]
        converterFn2 = NW_FORMAT[format][1][1]
        flexible1 = NW_FORMAT[format][0][2]
    else:
        container1 = NW_FORMAT[format][2][0] #support/name
        container2 = NW_FORMAT[format][3][0] #dist
        converterFn1 = NW_FORMAT[format][2][1]
        converterFn2 = NW_FORMAT[format][3][1]
        flexible1 = NW_FORMAT[format][2][2]

    if converterFn1 == str:
        try:
            if not quoted_names:
                FIRST_PART = re.sub("["+_ILEGAL_NEWICK_CHARS+"]", "_", \
                                    str(getattr(node, container1)))
            else:
                FIRST_PART = str(getattr(node, container1))
            if not FIRST_PART and container1 == 'name' and not flexible1:
                FIRST_PART = "NoName"

        except (AttributeError, TypeError):
            FIRST_PART = "?"

        FIRST_PART = name_formatter %FIRST_PART
        if quoted_names:
            #FIRST_PART = '"%s"' %FIRST_PART.decode('string_escape').replace('"', '\\"')
            FIRST_PART = '"%s"' %FIRST_PART

    elif converterFn1 is None:
        FIRST_PART = ""
    else:
        try:
            FIRST_PART = support_formatter %(converterFn2(getattr(node, container1)))
        except (ValueError, TypeError):
            FIRST_PART = "?"

    if converterFn2 == str:
        try:
            SECOND_PART = ":"+re.sub("["+_ILEGAL_NEWICK_CHARS+"]", "_", \
                                  str(getattr(node, container2)))
        except (ValueError, TypeError):
            SECOND_PART = ":?"
    elif converterFn2 is None:
        SECOND_PART = ""
    else:
        try:
            #SECOND_PART = ":%0.6f" %(converterFn2(getattr(node, container2)))
            SECOND_PART = ":%s" %(dist_formatter %(converterFn2(getattr(node, container2))))
        except (ValueError, TypeError):
            SECOND_PART = ":?"

    return "%s%s" %(FIRST_PART, SECOND_PART)


class NewickError(Exception):
    """Exception class designed for NewickIO errors."""
    def __init__(self, value):
        if value is None:
            value = ''
        value += "\nYou may want to check other newick loading flags like 'format' or 'quoted_node_names'."
        Exception.__init__(self, value)

def read_newick(newick, root_node=None, format=0, quoted_names=False):
    """ Reads a newick tree from either a string or a file, and returns
    an ETE tree structure.

    A previously existent node object can be passed as the root of the
    tree, which means that all its new children will belong to the same
    class as the root(This allows to work with custom TreeNode
    objects).

    You can also take advantage from this behaviour to concatenate
    several tree structures.
    """

    if root_node is None:
        #from ..coretype.tree import TreeNode
        root_node = TreeNode()

    if isinstance(newick, str):

        nw = newick
        matcher = compile_matchers(formatcode=format)
        nw = nw.strip()

        if not nw.startswith('(') and nw.endswith(';'):
            return _read_newick_from_string(nw, root_node, matcher, format, quoted_names)
        elif not nw.startswith('(') or not nw.endswith(';'):
            raise NewickError('Malformed newick tree structure.')
        else:
            return _read_newick_from_string(nw, root_node, matcher, format, quoted_names)

    else:
        raise NewickError("'newick' argument must be either a filename or a newick string.")

def _read_newick_from_string(nw, root_node, matcher, formatcode, quoted_names):
    """ Reads a newick string in the New Hampshire format. """

    if quoted_names:
        # Quoted text is mapped to references
        quoted_map = {}
        unquoted_nw = ''
        counter = 0
        for token in re.split(_QUOTED_TEXT_RE, nw):
            counter += 1
            if counter % 2 == 1 : # normal newick tree structure data
                unquoted_nw += token
            else: # quoted text, add to dictionary and replace with reference
                quoted_ref_id= _QUOTED_TEXT_PREFIX + str(int(counter/2))
                unquoted_nw += quoted_ref_id
                quoted_map[quoted_ref_id]=token[1:-1]  # without the quotes
        nw = unquoted_nw

    if not nw.startswith('(') and nw.endswith(';'):
        _read_node_data(nw[:-1], root_node, "single", matcher, format)
        if quoted_names:
            if root_node.name.startswith(_QUOTED_TEXT_PREFIX):
                root_node.name = quoted_map[root_node.name]
        return root_node

    if nw.count('(') != nw.count(')'):
        raise NewickError('Parentheses do not match. Broken tree structure?')

    # white spaces and separators are removed
    nw = re.sub("[\n\r\t]+", "", nw)

    current_parent = None
    # Each chunk represents the content of a parent node, and it could contain
    # leaves and closing parentheses.
    # We may find:
    # leaf, ..., leaf,
    # leaf, ..., leaf))),
    # leaf)), leaf, leaf))
    # leaf))
    # ) only if formatcode == 100

    for chunk in nw.split("(")[1:]:
        # If no node has been created so far, this is the root, so use the node.
        current_parent = root_node if current_parent is None else current_parent.add_child()

        subchunks = [ch.strip() for ch in chunk.split(",")]
        # We should expect that the chunk finished with a comma (if next chunk
        # is an internal sister node) or a subchunk containing closing parenthesis until the end of the tree.
        #[leaf, leaf, '']
        #[leaf, leaf, ')))', leaf, leaf, '']
        #[leaf, leaf, ')))', leaf, leaf, '']
        #[leaf, leaf, ')))', leaf), leaf, 'leaf);']
        if subchunks[-1] != '' and not subchunks[-1].endswith(';'):
            raise NewickError('Broken newick structure at: %s' %chunk)

        # lets process the subchunks. Every closing parenthesis will close a
        # node and go up one level.
        for i, leaf in enumerate(subchunks):
            if leaf.strip() == '' and i == len(subchunks) - 1:
                continue # "blah blah ,( blah blah"
            closing_nodes = leaf.split(")")

            # first part after splitting by ) always contain leaf info
            _read_node_data(closing_nodes[0], current_parent, "leaf", matcher, formatcode)

            # next contain closing nodes and data about the internal nodes.
            if len(closing_nodes)>1:
                for closing_internal in closing_nodes[1:]:
                    closing_internal =  closing_internal.rstrip(";")
                    # read internal node data and go up one level
                    _read_node_data(closing_internal, current_parent, "internal", matcher, formatcode)
                    current_parent = current_parent.up

    # references in node names are replaced with quoted text before returning
    if quoted_names:
        for node in root_node.traverse():
            if node.name.startswith(_QUOTED_TEXT_PREFIX):
                node.name = quoted_map[node.name]

    return root_node

def _parse_extra_features(node, NHX_string):
    """ Reads node's extra data form its NHX string. NHX uses this
    format:  [&&NHX:prop1=value1:prop2=value2] """
    NHX_string = NHX_string.replace("[&&NHX:", "")
    NHX_string = NHX_string.replace("]", "")
    for field in NHX_string.split(":"):
        try:
            pname, pvalue = field.split("=")
        except ValueError as e:
            raise NewickError('Invalid NHX format %s' %field)
        node.add_feature(pname, pvalue)

def compile_matchers(formatcode):
    matchers = {}
    for node_type in ["leaf", "single", "internal"]:
        if node_type == "leaf" or node_type == "single":
            container1 = NW_FORMAT[formatcode][0][0]
            container2 = NW_FORMAT[formatcode][1][0]
            converterFn1 = NW_FORMAT[formatcode][0][1]
            converterFn2 = NW_FORMAT[formatcode][1][1]
            flexible1 = NW_FORMAT[formatcode][0][2]
            flexible2 = NW_FORMAT[formatcode][1][2]
        else:
            container1 = NW_FORMAT[formatcode][2][0]
            container2 = NW_FORMAT[formatcode][3][0]
            converterFn1 = NW_FORMAT[formatcode][2][1]
            converterFn2 = NW_FORMAT[formatcode][3][1]
            flexible1 = NW_FORMAT[formatcode][2][2]
            flexible2 = NW_FORMAT[formatcode][3][2]

        if converterFn1 == str:
            FIRST_MATCH = "("+_NAME_RE+")"
        elif converterFn1 == float:
            FIRST_MATCH = "("+_FLOAT_RE+")"
        elif converterFn1 is None:
            FIRST_MATCH = '()'

        if converterFn2 == str:
            SECOND_MATCH = "(:"+_NAME_RE+")"
        elif converterFn2 == float:
            SECOND_MATCH = "(:"+_FLOAT_RE+")"
        elif converterFn2 is None:
            SECOND_MATCH = '()'

        if flexible1 and node_type != 'leaf':
            FIRST_MATCH += "?"
        if flexible2:
            SECOND_MATCH += "?"


        matcher_str= '^\s*%s\s*%s\s*(%s)?\s*$' % (FIRST_MATCH, SECOND_MATCH, _NHX_RE)
        compiled_matcher = re.compile(matcher_str)
        matchers[node_type] = [container1, container2, converterFn1, converterFn2, compiled_matcher]

    return matchers

def _read_node_data(subnw, current_node, node_type, matcher, formatcode):
    """ Reads a leaf node from a subpart of the original newick
    tree """

    if node_type == "leaf" or node_type == "single":
        if node_type == "leaf":
            node = current_node.add_child()
        else:
            node = current_node
    else:
        node = current_node

    subnw = subnw.strip()

    if not subnw and node_type == 'leaf' and formatcode != 100:
        raise NewickError('Empty leaf node found')
    elif not subnw:
        return

    container1, container2, converterFn1, converterFn2, compiled_matcher = matcher[node_type]
    data = re.match(compiled_matcher, subnw)
    if data:
        data = data.groups()
        # This prevents ignoring errors even in flexible nodes:
        if subnw and data[0] is None and data[1] is None and data[2] is None:
            raise NewickError("Unexpected newick format '%s'" %subnw)

        if data[0] is not None and data[0] != '':
            node.add_feature(container1, converterFn1(data[0].strip()))

        if data[1] is not None and data[1] != '':
            node.add_feature(container2, converterFn2(data[1][1:].strip()))

        if data[2] is not None \
                and data[2].startswith("[&&NHX"):
            _parse_extra_features(node, data[2])
    else:
        raise NewickError("Unexpected newick format '%s' " %subnw[0:50])
    return

def write_newick(rootnode, features=None, format=1, format_root_node=True,
                 is_leaf_fn=None, dist_formatter=None, support_formatter=None,
                 name_formatter=None, quoted_names=False):
    """ Iteratively export a tree structure and returns its NHX
    representation. """
    newick = []
    leaf = is_leaf_fn if is_leaf_fn else lambda n: not bool(n.children)
    for postorder, node in rootnode.iter_prepostorder(is_leaf_fn=is_leaf_fn):
        if postorder:
            newick.append(")")
            if node.up is not None or format_root_node:
                newick.append(format_node(node, "internal", format,
                                          dist_formatter=dist_formatter,
                                          support_formatter=support_formatter,
                                          name_formatter=name_formatter,
                                          quoted_names=quoted_names))
                newick.append(_get_features_string(node, features))
        else:
            if node is not rootnode and node != node.up.children[0]:
                newick.append(",")

            if leaf(node):
                newick.append(format_node(node, "leaf", format,
                                          dist_formatter=dist_formatter,
                                          support_formatter=support_formatter,
                                          name_formatter=name_formatter,
                                          quoted_names=quoted_names))
                newick.append(_get_features_string(node, features))
            else:
                newick.append("(")

    newick.append(";")
    return ''.join(newick)

def _get_features_string(self, features=None):
    """ Generates the extended newick string NHX with extra data about
    a node. """
    string = ""
    if features is None:
        features = []
    elif features == []:
        features = sorted(self.features)

    for pr in features:
        if hasattr(self, pr):
            raw = getattr(self, pr)
            if type(raw) in ITERABLE_TYPES:
                raw = '|'.join(map(str, raw))
            elif type(raw) == dict:
                raw = '|'.join(map(lambda x,y: "%s-%s" %(x, y), raw.items()))
            elif type(raw) == str:
                pass
            else:
                raw = str(raw)

            value = re.sub("["+_ILEGAL_NEWICK_CHARS+"]", "_", \
                             raw)
            if string != "":
                string +=":"
            string +="%s=%s"  %(pr, str(value))
    if string != "":
        string = "[&&NHX:"+string+"]"

    return string

'''
TREE CLASS FROM ETE3
'''

DEFAULT_COMPACT = False
DEFAULT_SHOWINTERNAL = False
DEFAULT_DIST = 1.0
DEFAULT_SUPPORT = 1.0
DEFAULT_NAME = ""

class TreeError(Exception):
    """
    A problem occurred during a TreeNode operation
    """
    def __init__(self, value=''):
        self.value = value
    def __str__(self):
        return repr(self.value)

class TreeNode(object):
    """
    TreeNode (Tree) class is used to store a tree structure. A tree
    consists of a collection of TreeNode instances connected in a
    hierarchical way. Trees can be loaded from the New Hampshire Newick
    format (newick).

    :argument newick: Path to the file containing the tree or, alternatively,
       the text string containing the same information.

    :argument 0 format: subnewick format

      .. table::

          ======  ==============================================
          FORMAT  DESCRIPTION
          ======  ==============================================
          0        flexible with support values
          1        flexible with internal node names
          2        all branches + leaf names + internal supports
          3        all branches + all names
          4        leaf branches + leaf names
          5        internal and leaf branches + leaf names
          6        internal branches + leaf names
          7        leaf branches + all names
          8        all names
          9        leaf names
          100      topology only
          ======  ==============================================

    :returns: a tree node object which represents the base of the tree.

    **Examples:**

    ::

        t1 = Tree() # creates an empty tree
        t2 = Tree('(A:1,(B:1,(C:1,D:1):0.5):0.5);')
        t3 = Tree('/home/user/myNewickFile.txt')
    """

    def _get_dist(self):
        return self._dist
    def _set_dist(self, value):
        try:
            self._dist = float(value)
        except ValueError:
            raise TreeError('node dist must be a float number')

    def _get_support(self):
        return self._support
    def _set_support(self, value):
        try:
            self._support = float(value)
        except ValueError:
            raise TreeError('node support must be a float number')

    def _get_up(self):
        return self._up
    def _set_up(self, value):
        if type(value) == type(self) or value is None:
            self._up = value
        else:
            raise TreeError("bad node_up type")

    def _get_children(self):
        return self._children
    def _set_children(self, value):
        if type(value) == list:
            for n in value:
                if type(n) != type(self):
                    raise TreeError("Incorrect child node type")
            self._children = value
        else:
            raise TreeError("Incorrect children type")


    #: Branch length distance to parent node. Default = 0.0
    dist = property(fget=_get_dist, fset=_set_dist)
    #: Branch support for current node
    support = property(fget=_get_support, fset=_set_support)
    #: Pointer to parent node
    up = property(fget=_get_up, fset=_set_up)
    #: A list of children nodes
    children = property(fget=_get_children, fset=_set_children)

    def __init__(self, newick=None, format=0, dist=None, support=None,
                 name=None, quoted_node_names=False):
        self._children = []
        self._up = None
        self._dist = DEFAULT_DIST
        self._support = DEFAULT_SUPPORT
        self.features = set([])
        # Add basic features
        self.features.update(["dist", "support", "name"])
        if dist is not None:
            self.dist = dist
        if support is not None:
            self.support = support

        self.name = name if name is not None else DEFAULT_NAME

        # Initialize tree
        if newick is not None:
            self._dist = 0.0
            read_newick(newick, root_node = self, format=format,
                        quoted_names=quoted_node_names)


    def __nonzero__(self):
        return True

    def __bool__(self):
        """
        Python3's equivalent of __nonzero__
        If this is not defined bool(class_instance) will call
        __len__ in python3
        """
        return True

    def __repr__(self):
        return "Tree node '%s' (%s)" %(self.name, hex(self.__hash__()))

    def __and__(self, value):
        """ This allows to execute tree&'A' to obtain the descendant node
        whose name is A"""
        value=str(value)
        try:
            first_match = next(self.iter_search_nodes(name=value))
            return first_match
        except StopIteration:
            raise TreeError("Node not found")

    def __add__(self, value):
        """ This allows to sum two trees."""
        # Should a make the sum with two copies of the original trees?
        if type(value) == self.__class__:
            new_root = self.__class__()
            new_root.add_child(self)
            new_root.add_child(value)
            return new_root
        else:
            raise TreeError("Invalid node type")

    def __str__(self):
        """ Print tree in newick format. """
        return self.get_ascii(compact=DEFAULT_COMPACT, \
                                show_internal=DEFAULT_SHOWINTERNAL)

    def __contains__(self, item):
        """ Check if item belongs to this node. The 'item' argument must
        be a node instance or its associated name."""
        if isinstance(item, self.__class__):
            return item in set(self.get_descendants())
        elif type(item)==str:
            return item in set([n.name for n in self.traverse()])

    def __len__(self):
        """Node len returns number of children."""
        return len(self.get_leaves())

    def __iter__(self):
        """ Iterator over leaf nodes"""
        return self.iter_leaves()

    def add_feature(self, pr_name, pr_value):
        """
        Add or update a node's feature.
        """
        setattr(self, pr_name, pr_value)
        self.features.add(pr_name)

    def add_features(self, **features):
        """
        Add or update several features. """
        for fname, fvalue in features.items():
            setattr(self, fname, fvalue)
            self.features.add(fname)

    def del_feature(self, pr_name):
        """
        Permanently deletes a node's feature.
        """
        if hasattr(self, pr_name):
            delattr(self, pr_name)
            self.features.remove(pr_name)

    # Topology management
    def add_child(self, child=None, name=None, dist=None, support=None):
        """
        Adds a new child to this node. If child node is not suplied
        as an argument, a new node instance will be created.

        :argument None child: the node instance to be added as a child.
        :argument None name: the name that will be given to the child.
        :argument None dist: the distance from the node to the child.
        :argument None support: the support value of child partition.

        :returns: The child node instance

        """
        if child is None:
            child = self.__class__()

        if name is not None:
            child.name = name
        if dist is not None:
            child.dist = dist
        if support is not None:
            child.support = support

        self.children.append(child)
        child.up = self
        return child

    def remove_child(self, child):
        """
        Removes a child from this node (parent and child
        nodes still exit but are no longer connected).
        """
        try:
            self.children.remove(child)
        except ValueError as e:
            raise TreeError("child not found")
        else:
            child.up = None
            return child

    def add_sister(self, sister=None, name=None, dist=None):
        """
        Adds a sister to this node. If sister node is not supplied
        as an argument, a new TreeNode instance will be created and
        returned.
        """
        if self.up is None:
            raise TreeError("A parent node is required to add a sister")
        else:
            return self.up.add_child(child=sister, name=name, dist=dist)

    def remove_sister(self, sister=None):
        """
        Removes a sister node. It has the same effect as
        **`TreeNode.up.remove_child(sister)`**

        If a sister node is not supplied, the first sister will be deleted
        and returned.

        :argument sister: A node instance

        :return: The node removed
        """
        sisters = self.get_sisters()
        if len(sisters) > 0:
            if sister is None:
                sister = sisters.pop(0)
            return self.up.remove_child(sister)

    def delete(self, prevent_nondicotomic=True, preserve_branch_length=False):
        """
        Deletes node from the tree structure. Notice that this method
        makes 'disappear' the node from the tree structure. This means
        that children from the deleted node are transferred to the
        next available parent.

        :param True prevent_nondicotomic: When True (default), delete
            function will be execute recursively to prevent
            single-child nodes.

        :param False preserve_branch_length: If True, branch lengths
            of the deleted nodes are transferred (summed up) to its
            parent's branch, thus keeping original distances among
            nodes.

        **Example:**

        ::

                / C
          root-|
               |        / B
                \--- H |
                        \ A

          > H.delete() will produce this structure:

                / C
               |
          root-|--B
               |
                \ A

        """
        parent = self.up
        if parent:
            if preserve_branch_length:
                if len(self.children) == 1:
                    self.children[0].dist += self.dist
                elif len(self.children) > 1:
                    parent.dist += self.dist

            for ch in self.children:
                parent.add_child(ch)

            parent.remove_child(self)

        # Avoids parents with only one child
        if prevent_nondicotomic and parent and\
              len(parent.children) < 2:
            parent.delete(prevent_nondicotomic=False,
                          preserve_branch_length=preserve_branch_length)


    def detach(self):
        """
        Detachs this node (and all its descendants) from its parent
        and returns the referent to itself.

        Detached node conserves all its structure of descendants, and can
        be attached to another node through the 'add_child' function. This
        mechanism can be seen as a cut and paste.
        """

        if self.up:
            self.up.children.remove(self)
            self.up = None
        return self


    def prune(self, nodes, preserve_branch_length=False):
        """Prunes the topology of a node to conserve only the selected list of leaf
        internal nodes. The minimum number of nodes that conserve the
        topological relationships among the requested nodes will be
        retained. Root node is always conserved.

        :var nodes: a list of node names or node objects that should be retained

        :param False preserve_branch_length: If True, branch lengths
          of the deleted nodes are transferred (summed up) to its
          parent's branch, thus keeping original distances among
          nodes.

        **Examples:**

        ::

          t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', format=1)
          t1.prune(['A', 'B'])


          #                /-A
          #          /D /C|
          #       /F|      \-B
          #      |  |
          #    /H|   \-E
          #   |  |                        /-A
          #-root  \-G                 -root
          #   |                           \-B
          #   |   /-I
          #    \K|
          #       \-J



          t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', format=1)
          t1.prune(['A', 'B', 'C'])

          #                /-A
          #          /D /C|
          #       /F|      \-B
          #      |  |
          #    /H|   \-E
          #   |  |                              /-A
          #-root  \-G                  -root- C|
          #   |                                 \-B
          #   |   /-I
          #    \K|
          #       \-J



          t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', format=1)
          t1.prune(['A', 'B', 'I'])


          #                /-A
          #          /D /C|
          #       /F|      \-B
          #      |  |
          #    /H|   \-E                    /-I
          #   |  |                      -root
          #-root  \-G                      |   /-A
          #   |                             \C|
          #   |   /-I                          \-B
          #    \K|
          #       \-J

          t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', format=1)
          t1.prune(['A', 'B', 'F', 'H'])

          #                /-A
          #          /D /C|
          #       /F|      \-B
          #      |  |
          #    /H|   \-E
          #   |  |                              /-A
          #-root  \-G                -root-H /F|
          #   |                                 \-B
          #   |   /-I
          #    \K|
          #       \-J

        """
        def cmp_nodes(x, y):
            # if several nodes are in the same path of two kept nodes,
            # only one should be maintained. This prioritize internal
            # nodes that are already in the to_keep list and then
            # deeper nodes (closer to the leaves).
            if n2depth[x] > n2depth[y]:
                return -1
            elif n2depth[x] < n2depth[y]:
                return 1
            else:
                return 0

        to_keep = set(_translate_nodes(self, *nodes))
        start, node2path = self.get_common_ancestor(to_keep, get_path=True)
        to_keep.add(self)

        # Calculate which kept nodes are visiting the same nodes in
        # their path to the common ancestor.
        n2count = {}
        n2depth = {}
        for seed, path in node2path.items():
            for visited_node in path:
                if visited_node not in n2depth:
                    depth = visited_node.get_distance(start, topology_only=True)
                    n2depth[visited_node] = depth
                if visited_node is not seed:
                    n2count.setdefault(visited_node, set()).add(seed)

        # if several internal nodes are in the path of exactly the same kept
        # nodes, only one (the deepest) should be maintain.
        visitors2nodes = {}
        for node, visitors in n2count.items():
            # keep nodes connection at least two other nodes
            if len(visitors)>1:
                visitor_key = frozenset(visitors)
                visitors2nodes.setdefault(visitor_key, set()).add(node)

        for visitors, nodes in visitors2nodes.items():
            if not (to_keep & nodes):
                sorted_nodes = sorted(nodes, key=cmp_to_key(cmp_nodes))
                to_keep.add(sorted_nodes[0])

        for n in self.get_descendants('postorder'):
            if n not in to_keep:
                if preserve_branch_length:
                    if len(n.children) == 1:
                        n.children[0].dist += n.dist
                    elif len(n.children) > 1 and n.up:
                        n.up.dist += n.dist

                n.delete(prevent_nondicotomic=False)


    def swap_children(self):
        """
        Swaps current children order.
        """
        if len(self.children)>1:
            self.children.reverse()


    # #####################
    # Tree traversing
    # #####################


    def get_children(self):
        """
        Returns an independent list of node's children.
        """
        return [ch for ch in self.children]

    def get_sisters(self):
        """
        Returns an independent list of sister nodes.
        """
        if self.up is not None:
            return [ch for ch in self.up.children if ch!=self]
        else:
            return []

    def iter_leaves(self, is_leaf_fn=None):
        """
        Returns an iterator over the leaves under this node.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        for n in self.traverse(strategy="preorder", is_leaf_fn=is_leaf_fn):
            if not is_leaf_fn:
                if n.is_leaf():
                    yield n
            else:
                if is_leaf_fn(n):
                    yield n

    def get_leaves(self, is_leaf_fn=None):
        """
        Returns the list of terminal nodes (leaves) under this node.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        return [n for n in self.iter_leaves(is_leaf_fn=is_leaf_fn)]

    def iter_leaf_names(self, is_leaf_fn=None):
        """
        Returns an iterator over the leaf names under this node.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        for n in self.iter_leaves(is_leaf_fn=is_leaf_fn):
            yield n.name

    def get_leaf_names(self, is_leaf_fn=None):
        """
        Returns the list of terminal node names under the current
        node.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        return [name for name in self.iter_leaf_names(is_leaf_fn=is_leaf_fn)]

    def iter_descendants(self, strategy="levelorder", is_leaf_fn=None):
        """
        Returns an iterator over all descendant nodes.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        for n in self.traverse(strategy=strategy, is_leaf_fn=is_leaf_fn):
            if n is not self:
                yield n

    def get_descendants(self, strategy="levelorder", is_leaf_fn=None):
        """
        Returns a list of all (leaves and internal) descendant nodes.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        return [n for n in self.iter_descendants(strategy=strategy, \
                                                 is_leaf_fn=is_leaf_fn)]

    def traverse(self, strategy="levelorder", is_leaf_fn=None):
        """
        Returns an iterator to traverse the tree structure under this
        node.

        :argument "levelorder" strategy: set the way in which tree
           will be traversed. Possible values are: "preorder" (first
           parent and then children) 'postorder' (first children and
           the parent) and "levelorder" (nodes are visited in order
           from root to leaves)

        :argument None is_leaf_fn: If supplied, ``is_leaf_fn``
           function will be used to interrogate nodes about if they
           are terminal or internal. ``is_leaf_fn`` function should
           receive a node instance as first argument and return True
           or False. Use this argument to traverse a tree by
           dynamically collapsing internal nodes matching
           ``is_leaf_fn``.
        """
        if strategy=="preorder":
            return self._iter_descendants_preorder(is_leaf_fn=is_leaf_fn)
        elif strategy=="levelorder":
            return self._iter_descendants_levelorder(is_leaf_fn=is_leaf_fn)
        elif strategy=="postorder":
            return self._iter_descendants_postorder(is_leaf_fn=is_leaf_fn)

    def iter_prepostorder(self, is_leaf_fn=None):
        """
        Iterate over all nodes in a tree yielding every node in both
        pre and post order. Each iteration returns a postorder flag
        (True if node is being visited in postorder) and a node
        instance.
        """
        to_visit = [self]
        if is_leaf_fn is not None:
            _leaf = is_leaf_fn
        else:
            _leaf = self.__class__.is_leaf

        while to_visit:
            node = to_visit.pop(-1)
            try:
                node = node[1]
            except TypeError:
                # PREORDER ACTIONS
                yield (False, node)
                if not _leaf(node):
                    # ADD CHILDREN
                    to_visit.extend(reversed(node.children + [[1, node]]))
            else:
                #POSTORDER ACTIONS
                yield (True, node)

    def _iter_descendants_postorder(self, is_leaf_fn=None):
        to_visit = [self]
        if is_leaf_fn is not None:
            _leaf = is_leaf_fn
        else:
            _leaf = self.__class__.is_leaf

        while to_visit:
            node = to_visit.pop(-1)
            try:
                node = node[1]
            except TypeError:
                # PREORDER ACTIONS
                if not _leaf(node):
                    # ADD CHILDREN
                    to_visit.extend(reversed(node.children + [[1, node]]))
                else:
                    yield node
            else:
                #POSTORDER ACTIONS
                yield node

    def _iter_descendants_levelorder(self, is_leaf_fn=None):
        """
        Iterate over all desdecendant nodes.
        """
        tovisit = deque([self])
        while len(tovisit)>0:
            node = tovisit.popleft()
            yield node
            if not is_leaf_fn or not is_leaf_fn(node):
                tovisit.extend(node.children)

    def _iter_descendants_preorder(self, is_leaf_fn=None):
        """
        Iterator over all descendant nodes.
        """
        to_visit = deque()
        node = self
        while node is not None:
            yield node
            if not is_leaf_fn or not is_leaf_fn(node):
                to_visit.extendleft(reversed(node.children))
            try:
                node = to_visit.popleft()
            except:
                node = None

    def iter_ancestors(self):
        '''versionadded: 2.2

        Iterates over the list of all ancestor nodes from current node
        to the current tree root.

        '''
        node = self
        while node.up is not None:
            yield node.up
            node = node.up

    def get_ancestors(self):
        '''versionadded: 2.2

        Returns the list of all ancestor nodes from current node to
        the current tree root.

        '''
        return [n for n in self.iter_ancestors()]

    def write(self, features=None, outfile=None, format=0, is_leaf_fn=None,
              format_root_node=False, dist_formatter=None, support_formatter=None,
              name_formatter=None, quoted_node_names=False):
        """
        Returns the newick representation of current node. Several
        arguments control the way in which extra data is shown for
        every node:

        :argument features: a list of feature names to be exported
          using the Extended Newick Format (i.e. features=["name",
          "dist"]). Use an empty list to export all available features
          in each node (features=[])

        :argument outfile: writes the output to a given file

        :argument format: defines the newick standard used to encode the
          tree. See tutorial for details.

        :argument False format_root_node: If True, it allows features
          and branch information from root node to be exported as a
          part of the newick text string. For newick compatibility
          reasons, this is False by default.

        :argument is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.

        **Example:**

        ::

             t.write(features=["species","name"], format=1)

        """

        nw = write_newick(self, features=features, format=format,
                          is_leaf_fn=is_leaf_fn,
                          format_root_node=format_root_node,
                          dist_formatter=dist_formatter,
                          support_formatter=support_formatter,
                          name_formatter=name_formatter,
                          quoted_names=quoted_node_names)

        if outfile is not None:
            with open(outfile, "w") as OUT:
                OUT.write(nw)
        else:
            return nw

    def get_tree_root(self):
        """
        Returns the absolute root node of current tree structure.
        """
        root = self
        while root.up is not None:
            root = root.up
        return root

    def get_common_ancestor(self, *target_nodes, **kargs):
        """
        Returns the first common ancestor between this node and a given
        list of 'target_nodes'.

        **Examples:**

        ::

          t = tree.Tree("(((A:0.1, B:0.01):0.001, C:0.0001):1.0[&&NHX:name=common], (D:0.00001):0.000001):2.0[&&NHX:name=root];")
          A = t.get_descendants_by_name("A")[0]
          C = t.get_descendants_by_name("C")[0]
          common =  A.get_common_ancestor(C)
          print common.name

        """

        get_path = kargs.get("get_path", False)

        if len(target_nodes) == 1 and type(target_nodes[0]) \
                in set([set, tuple, list, frozenset]):
            target_nodes = target_nodes[0]

        # Convert node names into node instances
        target_nodes = _translate_nodes(self, *target_nodes)


        if type(target_nodes) != list:
            # If only one node is provided and is the same as the seed node,
            # return itself
            if target_nodes is self:
                if get_path:
                    return self, {}
                else:
                    return self
            else:
                #Otherwise find the common ancestor of current seed node and
                #the target_node provided
                target_nodes = [target_nodes, self]

        n2path = {}
        reference = []
        ref_node = None
        for n in target_nodes:
            current = n
            while current:
                n2path.setdefault(n, set()).add(current)
                if not ref_node:
                    reference.append(current)
                current = current.up
            if not ref_node:
                ref_node = n

        common = None
        for n in reference:
            broken = False
            for node, path in n2path.items():
                if node is not ref_node and n not in path:
                    broken = True
                    break

            if not broken:
                common = n
                break
        if not common:
            raise TreeError("Nodes are not connected!")

        if get_path:
            return common, n2path
        else:
            return common

    def iter_search_nodes(self, **conditions):
        """
        Search nodes in an iterative way. Matches are yielded as they
        are being found. This avoids needing to scan the full tree
        topology before returning the first matches. Useful when
        dealing with huge trees.
        """

        for n in self.traverse():
            conditions_passed = 0
            for key, value in conditions.items():
                if hasattr(n, key) and getattr(n, key) == value:
                    conditions_passed +=1
            if conditions_passed == len(conditions):
                yield n

    def search_nodes(self, **conditions):
        """
        Returns the list of nodes matching a given set of conditions.

        **Example:**

        ::

          tree.search_nodes(dist=0.0, name="human")

        """
        matching_nodes = []
        for n in self.iter_search_nodes(**conditions):
            matching_nodes.append(n)
        return matching_nodes

    def get_leaves_by_name(self, name):
        """
        Returns a list of leaf nodes matching a given name.
        """
        return self.search_nodes(name=name, children=[])

    def is_leaf(self):
        """
        Return True if current node is a leaf.
        """
        return len(self.children) == 0

    def is_root(self):
        """
        Returns True if current node has no parent
        """
        if self.up is None:
            return True
        else:
            return False

    # ###########################
    # Distance related functions
    # ###########################
    def get_distance(self, target, target2=None, topology_only=False):
        """
        Returns the distance between two nodes. If only one target is
        specified, it returns the distance between the target and the
        current node.

        :argument target: a node within the same tree structure.

        :argument target2: a node within the same tree structure. If
          not specified, current node is used as target2.

        :argument False topology_only: If set to True, distance will
          refer to the number of nodes between target and target2.

        :returns: branch length distance between target and
          target2. If topology_only flag is True, returns the number
          of nodes between target and target2.

        """

        if target2 is None:
            target2 = self
            root = self.get_tree_root()
        else:
            # is target node under current node?
            root = self

        target, target2 = _translate_nodes(root, target, target2)
        ancestor = root.get_common_ancestor(target, target2)

        dist = 0.0
        for n in [target2, target]:
            current = n
            while current != ancestor:
                if topology_only:
                    if  current!=target:
                        dist += 1
                else:
                    dist += current.dist
                current = current.up
        return dist

    def get_farthest_node(self, topology_only=False):
        """
        Returns the node's farthest descendant or ancestor node, and the
        distance to it.

        :argument False topology_only: If set to True, distance
          between nodes will be referred to the number of nodes
          between them. In other words, topological distance will be
          used instead of branch length distances.

        :return: A tuple containing the farthest node referred to the
          current node and the distance to it.

        """
        # Init farthest node to current farthest leaf
        farthest_node, farthest_dist = self.get_farthest_leaf(topology_only=topology_only)

        prev = self
        cdist = 0.0 if topology_only else prev.dist
        current = prev.up
        while current is not None:
            for ch in current.children:
                if ch != prev:
                    if not ch.is_leaf():
                        fnode, fdist = ch.get_farthest_leaf(topology_only=topology_only)
                    else:
                        fnode = ch
                        fdist = 0
                    if topology_only:
                        fdist += 1.0
                    else:
                        fdist += ch.dist
                    if cdist+fdist > farthest_dist:
                        farthest_dist = cdist + fdist
                        farthest_node = fnode
            prev = current
            if topology_only:
                cdist += 1
            else:
                cdist  += prev.dist
            current = prev.up
        return farthest_node, farthest_dist

    def _get_farthest_and_closest_leaves(self, topology_only=False, is_leaf_fn=None):
        # if called from a leaf node, no necessary to compute
        if (is_leaf_fn and is_leaf_fn(self)) or self.is_leaf():
            return self, 0.0, self, 0.0

        min_dist = None
        min_node = None
        max_dist = None
        max_node = None
        d = 0.0
        for post, n in self.iter_prepostorder(is_leaf_fn=is_leaf_fn):
            if n is self:
                continue
            if post:
                d -= n.dist if not topology_only else 1.0
            else:
                if (is_leaf_fn and is_leaf_fn(n)) or n.is_leaf():
                    total_d = d + n.dist if not topology_only else d
                    if min_dist is None or total_d < min_dist:
                        min_dist = total_d
                        min_node = n
                    if max_dist is None or total_d > max_dist:
                        max_dist = total_d
                        max_node = n
                else:
                    d += n.dist if not topology_only else 1.0
        return min_node, min_dist, max_node, max_dist


    def get_farthest_leaf(self, topology_only=False, is_leaf_fn=None):
        """
        Returns node's farthest descendant node (which is always a leaf), and the
        distance to it.

        :argument False topology_only: If set to True, distance
          between nodes will be referred to the number of nodes
          between them. In other words, topological distance will be
          used instead of branch length distances.

        :return: A tuple containing the farthest leaf referred to the
          current node and the distance to it.
        """
        min_node, min_dist, max_node, max_dist = self._get_farthest_and_closest_leaves(
        topology_only=topology_only, is_leaf_fn=is_leaf_fn)
        return max_node, max_dist

    def get_closest_leaf(self, topology_only=False, is_leaf_fn=None):
        """Returns node's closest descendant leaf and the distance to
        it.

        :argument False topology_only: If set to True, distance
          between nodes will be referred to the number of nodes
          between them. In other words, topological distance will be
          used instead of branch length distances.

        :return: A tuple containing the closest leaf referred to the
          current node and the distance to it.

        """
        min_node, min_dist, max_node, max_dist = self._get_farthest_and_closest_leaves(
            topology_only=topology_only, is_leaf_fn=is_leaf_fn)

        return min_node, min_dist


    def copy(self, method="deepcopy"):
        """.. versionadded: 2.1

        Returns a copy of the current node.

        :var cpickle method: Protocol used to copy the node
        structure. The following values are accepted:

           - "newick": Tree topology, node names, branch lengths and
             branch support values will be copied by as represented in
             the newick string (copy by newick string serialisation).

           - "newick-extended": Tree topology and all node features
             will be copied based on the extended newick format
             representation. Only node features will be copied, thus
             excluding other node attributes. As this method is also
             based on newick serialisation, features will be converted
             into text strings when making the copy.

           - "cpickle": The whole node structure and its content is
             cloned based on cPickle object serialisation (slower, but
             recommended for full tree copying)

           - "deepcopy": The whole node structure and its content is
             copied based on the standard "copy" Python functionality
             (this is the slowest method but it allows to copy complex
             objects even if attributes point to lambda functions,
             etc.)

        """
        method = method.lower()
        if method=="newick":
            new_node = self.__class__(self.write(features=["name"], format_root_node=True))
        elif method=="newick-extended":
            self.write(features=[], format_root_node=True)
            new_node = self.__class__(self.write(features=[]))
        elif method == "deepcopy":
            parent = self.up
            self.up = None
            new_node = copy.deepcopy(self)
            self.up = parent
        else:
            raise TreeError("Invalid copy method")

        return new_node

    def _asciiArt(self, char1='-', show_internal=True, compact=False, attributes=None):
        """
        Returns the ASCII representation of the tree.

        Code based on the PyCogent GPL project.
        """
        if not attributes:
            attributes = ["name"]
        node_name = ', '.join(map(str, [getattr(self, v) for v in attributes if hasattr(self, v)]))

        LEN = max(3, len(node_name) if not self.children or show_internal else 3)
        PAD = ' ' * LEN
        PA = ' ' * (LEN-1)
        if not self.is_leaf():
            mids = []
            result = []
            for c in self.children:
                if len(self.children) == 1:
                    char2 = '/'
                elif c is self.children[0]:
                    char2 = '/'
                elif c is self.children[-1]:
                    char2 = '\\'
                else:
                    char2 = '-'
                (clines, mid) = c._asciiArt(char2, show_internal, compact, attributes)
                mids.append(mid+len(result))
                result.extend(clines)
                if not compact:
                    result.append('')
            if not compact:
                result.pop()
            (lo, hi, end) = (mids[0], mids[-1], len(result))
            prefixes = [PAD] * (lo+1) + [PA+'|'] * (hi-lo-1) + [PAD] * (end-hi)
            mid = int((lo + hi) / 2)
            prefixes[mid] = char1 + '-'*(LEN-2) + prefixes[mid][-1]
            result = [p+l for (p,l) in zip(prefixes, result)]
            if show_internal:
                stem = result[mid]
                result[mid] = stem[0] + node_name + stem[len(node_name)+1:]
            return (result, mid)
        else:
            return ([char1 + '-' + node_name], 0)

    def get_ascii(self, show_internal=True, compact=False, attributes=None):
        """
        Returns a string containing an ascii drawing of the tree.

        :argument show_internal: includes internal edge names.
        :argument compact: use exactly one line per tip.

        :param attributes: A list of node attributes to shown in the
            ASCII representation.

        """
        (lines, mid) = self._asciiArt(show_internal=show_internal,
                                      compact=compact, attributes=attributes)
        return '\n'+'\n'.join(lines)

    def iter_edges(self, cached_content = None):
        '''
        .. versionadded:: 2.3

        Iterate over the list of edges of a tree. Each edge is represented as a
        tuple of two elements, each containing the list of nodes separated by
        the edge.
        '''

        if not cached_content:
            cached_content = self.get_cached_content()
        all_leaves = cached_content[self]
        for n, side1 in cached_content.items():
            yield (side1, all_leaves-side1)

    def get_edges(self, cached_content = None):
        '''
        .. versionadded:: 2.3

        Returns the list of edges of a tree. Each edge is represented as a
        tuple of two elements, each containing the list of nodes separated by
        the edge.
        '''

        return [edge for edge in self.iter_edges(cached_content)]
    

def _translate_nodes(root, *nodes):
    name2node = dict([ [n, None] for n in nodes if type(n) is str])
    if name2node:
        for n in root.traverse():
            if n.name in name2node:
                if name2node[n.name] is not None:
                    raise TreeError("Ambiguous node name: "+str(n.name))
                else:
                    name2node[n.name] = n

    if None in list(name2node.values()):
        notfound = [key for key, value in name2node.items() if value is None]
        raise ValueError("Node names not found: "+str(notfound))

    valid_nodes = []
    for n in nodes:
        if type(n) is not str:
            if type(n) is not root.__class__ :
                raise TreeError("Invalid target node: "+str(n))
            else:
                valid_nodes.append(n)

    valid_nodes.extend(list(name2node.values()))
    if len(valid_nodes) == 1:
        return valid_nodes[0]
    else:
        return valid_nodes

# Alias
#: .. currentmodule:: ete3
Tree = TreeNode