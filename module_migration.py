'''
FUNCTIONS RELATING TO PROCESSING MIGRATION EVENTS
'''

import sys
import io

import pandas as pd

from module_helper import stripall
from module_tree import get_node_pairs_to_modify
from module_tree import get_attribute_filtered_tree

# load migration parameter into migration dataframe
def read_specified_mig_pattern(
        migration
        ):
    
    try: 
        
        # remove brackets, clean out tabs, and split at commas
        migration = migration[1:-1]
        migration = migration.replace("\t","")
        migration = migration.split(",")

        # replace bidirectional migrations with two unidirectional migrations, swap directions if specified backwards, and rejoin text
        migration_replaced = []
        for mig in migration:
            two_way = mig.split("<->") # A <-> B replaced with A -> B and B -> A
            reverse_dir = mig.split("<-") # A <- B replaced with B -> A
            
            if len(two_way) == 2:
                migration_replaced.append(f"{two_way[0]}->{two_way[1]}")
                migration_replaced.append(f"{two_way[1]}->{two_way[0]}")
            
            elif len(reverse_dir) == 2:
                migration_replaced.append(f"{reverse_dir[1]}->{reverse_dir[0]}")

            else:
                migration_replaced.append(mig)

            
        
        lines_text = "\n".join(migration_replaced)

        # read text into dataframe as it were a csv, with the separator '->'
        buf = io.StringIO(lines_text)
        df = pd.read_csv(buf, sep = "->", header=None, engine='python')
        df = df.applymap(lambda x: stripall(x))
        df = df.rename({0: 'source', 1: 'destination'}, axis=1)

    except:
        sys.exit("Error: migration parameter incorrectly formatted. Refer to manual for further details.")

    return df

# remaps the migration event between two populations to the populations currently accepted as species
def remap_migrate(
        tree, 
        source_name, 
        dest_name
        ):

    '''
    Migration events in the guide topology must be remapped when the species delimitation changes.
    consider the following scenario, with migration from A to B, and C to B:

    A      B     C
     \--> / <-- /
      \  /     /
       AB     /
        \    /
         \  /
         ABC

    When the topology is changed such that A and B are merged, the migration from A to B is dropped,
    and the migration from C to B is remapped, now going from C to AB

      AB <-- C
       \    /
        \  /
        ABC
    '''
    
    source_node = tree.search_nodes(name = source_name)[0]
    dest_node = tree.search_nodes(name = dest_name)[0]

    while (source_node.species == False and source_node.proposal == None):
        source_node = source_node.up

    while (dest_node.species == False and dest_node.proposal == None):
        dest_node = dest_node.up

    return source_node.name, dest_node.name


# append rows to the control file to infer migration parameters, depending on how the migration pattern was specified
def append_migrate_rows(
        tree, 
        mode, 
        ctl_file_name, 
        mig
        ):
    
    # when the migration pattern is for sister nodes, automatically generate the required migration rows
    if type(mig) == str:
        pair_list = get_node_pairs_to_modify(tree, mode)
        
        # top row corresponds to number of migration events
        txt = f'migration = {len(pair_list)*2}\n'
        
        # rest of the rows specify the populations between which migration occurs
        for pair in pair_list:
            txt += f'\t{pair[0].name} {pair[1].name}\n'
            txt += f'\t{pair[1].name} {pair[0].name}\n'

    # when the migration pattern is specified for certain nodes, get the resulting pattern
    else:
        tree = get_attribute_filtered_tree(tree, "population", newick=False)
        
        migration_list = []
        for index, row in mig.iterrows():
            remap = (remap_migrate(tree, row['source'], row['destination']))
            # if a migration event occurs between nodes with the same ancestor that is currently accepted as species, then that migration event is now intra-species, so it is not processed 
            if remap[0] != remap[1]:
                # if the migration event has already been added (this can be due to remapping)
                if remap not in migration_list:
                    migration_list.append(remap)

        # top row corresponds to number of migration events
        txt = f'migration = {len(migration_list)}\n'
        
        # rest of the rows specify the populations between which migration occurs
        for pair in migration_list:
            txt += f'\t{pair[0]} {pair[1]}\n'

    # write resulting lines to the end of the control file
    with open(ctl_file_name, "a") as myfile:
        myfile.write(txt)



# get the source of migration events to a node, and the destination of migration events from the node as lists
def get_node_migration_events(
        node, 
        migration_df
        ):

    # source      of migration events to   the node
    mig_source = migration_df[migration_df["destination"] == str(node.name)]["source"].to_list()
    # destination of migration events from the node
    mig_dest   = migration_df[migration_df["source"] == str(node.name)]["destination"].to_list()

    return mig_source, mig_dest


# check that all migration events a given node pair participate in only involve eachother
def check_migration_reciprocal(
        node_1, 
        node_2,
        migration_df
        ):

    '''
    If two node pairs A) dont participate in any migration events, or B) only participate in migration events involving
    eachother, the numerical formula for calculating gdi values 'gdi_numeric' may be used.
    '''

    # if no migration events occur, then there are no non-reciprocal migration events
    if migration_df is None:
        return True

    m_src_1, m_dest_1 = get_node_migration_events(node_1, migration_df)
    m_src_2, m_dest_2 = get_node_migration_events(node_2, migration_df)

    # if all the migration sources and destinations are shared
    if all(a.count(node_1.name) == len(a) for a in [m_src_2, m_dest_2]) and all(a.count(node_2.name) == len(a) for a in [m_src_1, m_dest_1]):
        

        return True

    else:

        return False
