#!/usr/bin/env python
# -*-coding: utf8 -*-

# Copyright 2019 by Balázs Ligeti (obalasz@gmail.com). All rights reserved.


"""

.. moduleauthor:: Balázs Ligeti <obalasz@gmail.com>

"""


class Clade(object):
    """
    Defining the taxonomy node. Each Node has an id
    """

    def __init__(self, tax_id, parent_taxid, tax_main_name, rank):
        """
        Setting up a clade
        :param tax_id: NCBI taxonomy id, number
        :param parent_taxid: The id of the parent node (root has no parent)
        :param tax_main_name: The primary name of the taxon
        :param rank: The rank of the taxon i.e. species, there is a special rank, namely strain
        that is defined as the leaves of the taxonomy tree, if no rank defined for the given node
        """
        self.tax_id = tax_id
        self.parent = None
        self.children = []
        self.tax_main_name = tax_main_name
        self.rank = rank
        self.assigned_reads = []
        self.parent_taxid = parent_taxid
        self.cum_count = 0
        self.direct_count = 0
        self.p_bernoulli_hit = 0
        self.p_value = 2
        self.abundance = 0
        self.Path_set = None
        self.height = None

    def __str__(self):
        own_str = "|id=" + str(self.tax_id) + ", name=" + self.tax_main_name + ", parent_id=" \
                  + str(self.parent_taxid) + "|" + ", rank=" + self.rank
        return own_str

    def __repr__(self):
        own_str = "|id=" + str(self.tax_id) + ", name=" + self.tax_main_name + ", parent_id=" \
                  + str(self.parent_taxid) + "|" + ", rank=" + self.rank
        return own_str

    def _get_children_tax_id_list(self):
        """
        Query the children of the given taxon.
        :return:
        """
        print("Get children list:")
        return ", ".join([str(c.tax_id) for c in self.children if c is not None])

    def get_simple_record_data_abundance_report(self):
        """
        Generates a dictionary, that contains the main information about the taxa
        :return:
        """
        record = {"taxon_name": self.tax_main_name,
                  "taxon_id": self.tax_id,
                  "rank": self.rank,
                  "direct_hits": self.direct_count,
                  "cumulative_hit": self.cum_count}
        return record

    def reset_hit_counts(self):
        """
        Setting the counts to the default values
        """
        self.cum_count = 0
        self.direct_count = 0


class PhageClade(Clade):
    """This is a special node for storing Phage (virus of bacteria) node.
    It has a special attribute, that links to its host
    """

    def __init__(self, tax_id, parent_taxid, tax_main_name, rank, hosts=None):
        """

        :param tax_id: NCBI taxonomy id, number
        :param parent_taxid: The id of the parent node (root has no parent)
        :param tax_main_name: The primary name of the taxon
        :param rank: The rank of the taxon i.e. species, there is a special rank, namely strain
        that is defined as the leaves of the taxonomy tree, if no rank defined for the given node
        :param hosts: List of taxonomy ids of bacteria that could be targeted by the virus
        :return:
        """
        self.hosts = hosts
        super(PhageClade, self).__init__(tax_id, parent_taxid, tax_main_name, rank)


class Taxonomy(object):
    """Specialized class for taxonomy. Tailered to NCBI taxonomy database"""

    def __init__(self, working_dir, taxonomy_subset_name=None, isloading_taxonomy=True):
        """
        Create a taxonomy object.
        :param working_dir: The directory where the source files where NCBI files are deposited.
        :param taxonomy_subset_name: The name of subset of the taxonomy (it is the prefix of
        names.dmp)
        """
        print("Creating a taxonomy")
        import os
        self.working_dir = working_dir
        if taxonomy_subset_name is None:
            self.file_taxon_names = os.path.join(self.working_dir, "names.dmp")
            self.file_taxon2names = os.path.join(self.working_dir, "tax2names.tsv")
            self.file_id2parent = os.path.join(self.working_dir, "nodes.dmp")
        else:
            self.file_taxon_names = os.path.join(self.working_dir,
                                                 taxonomy_subset_name + "_names.dmp")
            self.file_taxon2names = os.path.join(self.working_dir,
                                                 taxonomy_subset_name + "_tax2names.tsv")
            self.file_id2parent = os.path.join(self.working_dir, taxonomy_subset_name + "_nodes.dmp")

        self.root = None
        self.nodes = []  # Contains all the taxons
        self.id2node = {}
        self.id2names = {}
        self.id2all_names = {}
        self.taxon_name2tax_id = {}
        self.rank2nodes = {}
        self.tax_id2rank = {}
        self.tax_id2path = {}
        self.maximum_depth = -100
        if isloading_taxonomy:
            self.load_taxonomy()

    def _loading_test_clades(self):
        """
        Building a test taxonomy to asses the proper functionality of the tree
        :return:
        """
        cr1 = Clade(2, 1, 'cr1', 'no rank')  # Root a parent
        cr2 = Clade(3, 1, 'cr2', 'no rank')  # Root a parent
        c21 = Clade(4, 2, 'c21', 'no rank')  # parent: id=2
        self.root.children.extend([cr1, cr2])
        cr2.children.append(c21)
        self.nodes.extend([self.root, cr1, cr2, c21])

    def _loading_taxon_id2taxon_name_map(self):
        print("Loading NCBI tax_id2tax_name map!")
        with open(self.file_taxon_names) as fi:
            for line in fi:
                data = [attr.strip() for attr in line.split("|")]
                act_tax_id = int(data[0])
                if act_tax_id in self.id2all_names:
                    self.id2all_names[act_tax_id].append(data[1].lower())
                else:
                    self.id2all_names[act_tax_id] = [data[1].lower()]
                if "scientific name" in line:
                    line = [attr.strip() for attr in line.split("|")]
                    self.id2names[act_tax_id] = data[1].lower()


    def _reading_clades(self):
        """Reading the clade information from the NCBI taxonomy dump file"""
        print("Loading clades!")
        with open(self.file_id2parent) as f:
            for line in f:
                line = [i.strip() for i in line.split("|")]
                new_clade = Clade(int(line[0]), int(line[1]), self.id2names[int(line[0])], line[2])
                self.id2node[int(line[0])] = new_clade

    def _creating_taxonomy_tree(self):
        """
        Building the tree linking the nodes together
        :return:
        """
        for node in self.id2node.values():
            node.parent = self.id2node[node.parent_taxid]
            self.id2node[node.parent_taxid].children.append(node)
        self.root = self.id2node[1]
        self.id2node[1].parent = None
        self.id2node[1].parent_taxid = None
        self.root.children.remove(self.root)

    def _finding_strains(self):
        """
        Finding the strains, that are defined as the leaves of the taxonomy tree, if no rank is
        defined
        """
        strain_nodes = []
        print("Finding strain nodes!")
        for node in self.id2node.values():
            if not node.children and node.rank == 'no rank':
                strain_nodes.append(node)
        return strain_nodes

    def _finding_strains_having_species_ascendant(self):
        """
        Finding the strains, that are defined as the leaves of the taxonomy tree, if no rank is
        defined and the node has species ascendant
        """
        strain_nodes = []
        print("Finding strain nodes!")
        for node in self.id2node.values():
            if len(node.children) == 0:
                ascendant_nodes = self.get_ascendant_nodes(node.tax_id, "species")
                ascendant_nodes.remove(node)  # remove the element it self
                if len(ascendant_nodes) > 1:
                    strain_nodes.append(node)
        return strain_nodes

    def _creating_rank2clades_map(self, add_strains):
        """
        Building a mapping between the taxonomy level (genus, species) and clades/taxa
        :param add_strains: It is true if a special rank (not defined in the NCBI), should be added
        :return:
        """
        nodes_reads = [node for node in self.id2node.values()]
        levels = set([node.rank for node in nodes_reads])
        rank2nodes = {level: [] for level in levels}
        for node in nodes_reads:
            rank2nodes[node.rank].append(node)
        if add_strains:
            print("Adding strains!")
            rank2nodes['strain'] = self._finding_strains()
            print("Number of strains in the current dbs %d" % len(rank2nodes['strain']))
            for node in rank2nodes["strain"]:
                self.tax_id2rank[node.tax_id] = "strain"
                node.rank = "strain"

        self.rank2nodes = rank2nodes
        for rank, nodes in rank2nodes.items():
            for node in nodes:
                self.tax_id2rank[node.tax_id] = rank

    def _get_descendant_nodes_recursion(self, taxon, descendants):
        """
        Find all descendant taxa of a given taxon and returns with the list of reference of
        descdentants.
        :param taxon: The query clade
        :param descendants: list of descendant (it should be and empty list)
        :return:
        """
        for child in taxon.children:
            descendants.append(child)
            self._get_descendant_nodes_recursion(child, descendants)

    def _build_path2root(self):
        """
        Building taxon paths to root
        :return:
        """
        self.tax_id2path = {}
        print("Building path to root.")
        for tax_id in self.id2node:
            path = [taxon_node.tax_id for taxon_node in self.get_ascendant_nodes(tax_id)]
            self.tax_id2path[tax_id] = path

    def _assign_direct_hits(self, taxon_id2count):
        """
        Assign the hit counts to the taxa.
        If the taxon is not in the database, then it skip the information.
        It also modifies the state of the clades
        :param taxon_id2count: Contains the direct mapping between the taxon hits
        :return:
        """
        for taxon_id, count in taxon_id2count.items():
            if taxon_id in self.id2node:
                self.id2node[taxon_id].direct_count = count

    def _get_indexed_subtree_taxa_ids(self, taxa):
        """
        Find all the nodes that are ascendants of taxa and make an index of them
        :param taxa: list of taxa spanning the subtree
        :return: a dictionary of taxa ids of the subtree
        """
        non_null_counts_taxon_ids = self.get_ascendant_node_ids_of_set_of_nodes(taxa)
        non_null_counts_taxon_ids = dict(zip(non_null_counts_taxon_ids, [0 for i in non_null_counts_taxon_ids]))
        return non_null_counts_taxon_ids

    def _get_cumulative_hits(self, non_null_counts_taxon_ids):
        """
        Calculates the commulative sums for all the taxa
        :param non_null_counts_taxon_ids: Contains the taxa of the subtree
        :return: None
        """
        self._summing_subtree_hits(self.root, non_null_counts_taxon_ids)

    def _summing_subtree_hits(self, taxon, non_null_counts_taxon_ids):
        """
        Aggreagates all the direct hits assigned to nodes belonging to the subtree of a node
        :param taxon: reference to the taxon node. It is NOT an id, but the reference
        :type taxon: Clade
        :param non_null_counts_taxon_ids: Ids of nodes that potentatially having non-null count value. It gives only the interesting brances
        :return: with the cumulative count
        """
        cum_count = taxon.direct_count
        if len(taxon.children) > 0:
            for child in taxon.children:
                # ignoring non-interesting branches
                if child.tax_id in non_null_counts_taxon_ids:
                    cum_count += self._summing_subtree_hits(child, non_null_counts_taxon_ids)
        taxon.cum_count = cum_count
        return cum_count

    def _generate_report_table(self, non_null_counts_taxon_ids):
        """
        Generates a report for all taxa having at least 1 read assigned
        :return: A dictionary that contains the report data with the following format:
        """
        print("Generating report!")
        report = {}
        for taxon_id in non_null_counts_taxon_ids:
            path = self.get_ascendant_nodes(taxon_id)
            path.reverse()
            parent_ids = "|".join([str(parent.tax_id) for parent in path])
            parent_names = "|".join([str(parent.tax_main_name) for parent in path])
            record = self.id2node[taxon_id].get_simple_record_data_abundance_report()
            record["path_ids"] = parent_ids
            record["path"] = parent_names
            report[taxon_id] = record
        print("finished")
        return report

    def _generate_report_table_only_hits(self, non_null_counts_taxon_ids):
        """
        Generates a report for all taxa having at least 1 read assigned. (Only the cumulative hit is reported)
        :return: A dictionary that contains the report data with the following format:
        """
        report = {}
        for taxon_id in non_null_counts_taxon_ids:
            record_full = self.id2node[taxon_id].get_simple_record_data_abundance_report()
            report[taxon_id] = record_full["cumulative_hit"]
        return report

    def _removing_counts(self):
        """
        Setting a clear state for the tree
        :return:
        """
        for clade in self.id2node.values():
            clade.reset_hit_counts()


            ####################################################################################################
            ################################# PUBLIC FUNCTIONS  ################################################
            ####################################################################################################

    def get_descendant_node_ids(self, taxon_id):
        """
        Find all descendant taxa of a given taxon defined by its id and returns with the list of
        taxon ids of descendats
        :param taxon_id:  The query clade
        :return:
        """
        taxon = self.id2node[taxon_id]
        descendants = [taxon]
        self._get_descendant_nodes_recursion(taxon, descendants)
        descendant_tax_ids = [taxon.tax_id for taxon in descendants]
        return descendant_tax_ids

    def get_descendant_nodes(self, taxon):
        """
        Find all descendant taxa of a given taxon and returns with the list of reference of
        descdendants.
        :param taxon: The query clade
        :return: List of references to the descendant nodes
        """
        descendants = []
        return self._get_descendant_nodes_recursion(taxon, descendants)

    def get_ascendant_node_ids_of_set_of_nodes(self, taxa_ids):
        """
        Get the ids of all ascendant of a given set of taxon (defined by its ids) in the taxonomy
        :param taxa_ids: set of taxon
        :return: List of ids of nodes being ascendant of nodes in the set of taxa
        """
        taxa = []
        for taxon in taxa_ids:
            taxa.extend([taxon_node.tax_id for taxon_node in self.get_ascendant_nodes(taxon) if taxon in self.id2node])
        taxa = list(set(taxa))
        return taxa

    def get_ascendant_nodes(self, taxon_id, rank=None):
        """
        Get the ascendant nodes of a certain taxon. The ascendant nodes contains the query itself as
        first node ind list
        :param taxon_id:
        :param rank: The taxonomy level until the ascendants are listed. I.e. it is species,
        then the function lists all the taxon that are under the taxonomy level species.
        :return: List of ascendant nodes, the query clade is the first element in that list
        """
        ascendants = []
        act_node = self.id2node[taxon_id]
        while act_node is not None:
            ascendants.append(act_node)
            if rank is not None and act_node.rank is not None and act_node.rank == rank:
                return ascendants
            act_node = act_node.parent
        return ascendants

    def load_taxonomy(self, add_strains=True):
        """
        Loading the NCBI taxonomy (default)
        :param add_strains: it true if the leaves of taxonomy tree should handled as strains
        :return:
        """
        self._loading_taxon_id2taxon_name_map()
        self._reading_clades()
        self._creating_taxonomy_tree()
        self._creating_rank2clades_map(add_strains)

    def get_taxonomy_subset(self, output_dir, taxonomy_name, taxa, is_add_ascendants=True):
        """
        Filters the original taxonomy to a subset given in the taxa list and dump the filtered
        files into the give folder.
        :param output_dir: the directory where the taxonomy is dumped
        :param taxonomy_name: Name of the new taxonomy
        :param taxa: The list of taxon being filtered to
        :param is_add_ascendants: Forcing the tree property
        :return:
        """
        from os.path import join
        print("Filtering the original taxonomy!")
        file_new_nodes = join(output_dir, taxonomy_name + "_nodes.dmp")
        file_new_names = join(output_dir, taxonomy_name + "_names.dmp")
        fout_names = open(file_new_names, "w")
        fout_nodes = open(file_new_nodes, "w")
        if is_add_ascendants:
            new_taxa = {taxon: 0 for taxon in self.get_ascendant_node_ids_of_set_of_nodes(taxa)}
        else:
            new_taxa = {taxon: 0 for taxon in taxa}

        [fout_names.write(line) for line in open(self.file_taxon_names) if int(line.split(
            "|")[0].strip()) in new_taxa]
        [fout_nodes.write(line) for line in open(self.file_id2parent) if int(line.split(
            "|")[0].strip()) in new_taxa]
        fout_names.close()
        fout_nodes.close()

    def find_common_ancestor(self, taxa_ids):
        """
        Find the lowest common ancestor

        :param taxa_ids: to be tested taxon ids
        :return: common ancestor taxon id
        """
        if not self.tax_id2path:
            self._build_path2root()
        paths = []
        lowest_depth = -1
        for taxon_id in taxa_ids:
            if taxon_id in self.tax_id2path:
                act_path = self.tax_id2path[taxon_id]
                paths.append(act_path)
                path_length = len(act_path)
                if path_length < lowest_depth or lowest_depth < 0:
                    lowest_depth = path_length
        if not paths:
            return None
        elif len(paths) == 1:
            return paths[0][0]
        previous_common_ancestor = paths[0][-1]
        for i in range(-1, -1 * lowest_depth - 1, -1):
            act_common_ancestor = paths[0][i]
            for path in paths:
                if path[i] != act_common_ancestor:
                    return previous_common_ancestor
            previous_common_ancestor = act_common_ancestor
        return previous_common_ancestor

    def generate_cumulative_report(self, taxon_id2count, is_only_cum_hits=True):
        """
        Generates a table, that contains the taxon id and its cumulative hits.
        :param taxon_id2count: Contains the direct mapping between the taxon hits
        :param is_only_cum_hits: Report only the cumulative hits, but not the paths, parents, etc.
        :return:
        """
        try:
            taxa_ids_not_in_taxonomy = list( set(taxon_id2count.keys()) - set(self.id2node.keys()))
            taxa_ids_in_taxonomy = list(set(self.id2node.keys()) & set(taxon_id2count.keys()))
            # print("len original taxon: ", len(taxon_id2count))
            # print("taxa_ids_not_in_taxonomy: ", taxa_ids_not_in_taxonomy)
            # print("len(taxa_ids): ", len(taxa_ids_in_taxonomy))

            self._assign_direct_hits(taxon_id2count)
            non_null_counts_taxon_ids = self._get_indexed_subtree_taxa_ids(taxa_ids_in_taxonomy)
            self._get_cumulative_hits(non_null_counts_taxon_ids)
            if is_only_cum_hits:
                report = self._generate_report_table_only_hits(non_null_counts_taxon_ids)
            else:
                report = self._generate_report_table(non_null_counts_taxon_ids)
        finally:
            self._removing_counts()

        return report

    def get_pandas_dataframe(self, is_cumulative_paths=False):
        """
        It builds a pandas dataframe for the data.
        It contains the informations defined the clade records get_simple_record_data_abundance_report() function
        :param is_cumulative_paths: The datafram should contains the path as a string or not.
        :return:
        """
        import pandas as pd
        # using lists would be much better ....
        records = {taxon_id: self.id2node[taxon_id].get_simple_record_data_abundance_report() for taxon_id in self.id2node}
        if is_cumulative_paths:
            pass
#            adding path information to the data structure:
        data_table = pd.DataFrame.from_dict(records, orient='index')
        return data_table

    def get_taxon_dataframe(self):
        """
        It generates a taxon table, in which each record describe a taxon by it's NCBI taxon id, name (scientific, preffered), rank and parent taxon.
        By our definition the parent of the root is the root (self loop).
        :return: A pandas dataframe of records

        """
        import pandas as pd
        print('Creating a taxon dataframe!')
        columns = ['taxon_id', 'name', 'rank', 'parent_taxon_id']
        records = []
        print('Going through the hashed nodes!')
        for node in self.id2node.values():
            if node.parent is not None:
                records.append((node.tax_id, node.tax_main_name, node.rank, node.parent.tax_id))
            else:
                records.append((node.tax_id, node.tax_main_name, node.rank, 1))
        print('Finished')
        print('Creating DF')
        data_table = pd.DataFrame(records)
        data_table.columns = columns
        data_table.set_index('taxon_id', inplace=True, drop=False)
        return data_table


    def get_taxon_upper_taxa_dataframe(self):
        """
         It generates a taxon -> ascendants taxa table, in which each record describe a taxon's ascendant node.
         Each record contains the query taxon and it's ascendant id, the rank of the ascendant node
        :return: A pandas dataframe of records
        """
        import pandas as pd
        print("Generating the taxon -> uppertaxon table")
        columns = ['taxon_id', 'asc_rank', 'asc_taxon_id']
        print('Building table for all nodes\'s path to the root node!')
        self._build_path2root()
        print('Finished!')
        records = []
        print('Going through the node_id -> path -> node_id hash!')
        for taxon_id, path in self.tax_id2path.items():
            for asc_taxon_id in path:
                if self.id2node[asc_taxon_id].parent is not None:
                    records.append((taxon_id,self.id2node[asc_taxon_id].rank, asc_taxon_id))
        print('Finished')

        data_table = pd.DataFrame(records)
        data_table.columns = columns
        return data_table



