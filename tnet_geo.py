#!/usr/bin/python3
"""
Copyright (C) 2020 Saurav Dhar (saurav.dhar@uconn.edu),
Ion Mandoiu (ion.mandoiu@uconn.edu), and Mukul S. Bansal
(mukul.bansal@uconn.edu).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# Library imports
from Bio import Phylo
from collections import Counter
from csv import DictReader
import numpy as np
import copy, sys, os
import argparse, json

# Global variables
args = None
strain_country = {}
node_name_date = {}
leaf_strain = {}
hosts = []
score = {}
left_score = {}
right_score = {}
solution_count = {}
transmission_edges = []

def initialize_tree(input_file):
	input_tree = Phylo.read(input_file, 'newick')
	input_tree.rooted = True

	if not input_tree.is_bifurcating():
		raise IndexError("Input tree is not bifurcating.")

	return input_tree

def read_metadata_initialize_global_variables():
	with open(args.metadata, 'r') as read_obj:
		csv_dict_reader = DictReader(read_obj)
		for row in csv_dict_reader:
			node_name_date[row['strain']] = row['date']
			if row['country']:
				strain_country[row['strain']] = row['country'].replace(' ', '')

def initialize_leaf_nodes(rooted_tree):
	global hosts

	for terminal in rooted_tree.get_terminals():
		leaf_strain[terminal] = terminal.name
		if args.metadata:
			terminal.name = strain_country[terminal.name]
		else:
			terminal.name = terminal.name.split('_')[0]

		if terminal.name not in hosts:
			hosts.append(terminal.name)

	for terminal in rooted_tree.get_terminals():
		temp = []
		count = []
		for host in hosts:
			if host == terminal.name:
				temp.append(0)
				count.append(1)
			else:
				temp.append(9999999999)
				count.append(0)

		score[terminal] = temp
		solution_count[terminal] = count

def initialize_score_count(node):
	l_score = score[node.clades[0]].copy()
	r_score = score[node.clades[1]].copy()
	l_count = solution_count[node.clades[0]].copy()
	r_count = solution_count[node.clades[1]].copy()

	length = len(l_score)
	temp_score = []
	temp_left = []
	temp_right = []
	temp_count = []

	for i in range(length):
		l_score[i] -= 1
		left_count = 0
		min_left = min(l_score)
		for j in range(length):
			if l_score[j] == min_left:
				left_count += l_count[j]

		r_score[i] -= 1
		right_count = 0
		min_right = min(r_score)
		for j in range(length):
			if r_score[j] == min_right:
				right_count += r_count[j]

		temp_score.append(min_left + min_right + 2)
		temp_left.append(min_left)
		temp_right.append(min_right)
		temp_count.append(left_count * right_count)
		l_score[i] += 1
		r_score[i] += 1

	score[node] = temp_score
	left_score[node] = temp_left
	right_score[node] = temp_right
	solution_count[node] = temp_count

def initialize_internal_nodes(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'postorder'):
		initialize_score_count(nonterminal)

def get_host_from_count(count):
	if args.maxprob:
		max_count = max(count)
		for i in range(len(count)):
			if count[i] != max_count:
				count[i] = 0

	probs = [float(i)/sum(count) for i in count]
	ch = np.random.choice(len(probs), p = probs)
	return hosts[ch]

def choose_root_host(root_node):
	probs = []
	min_score = min(score[root_node])
	for i in range(len(score[root_node])):
		if score[root_node][i] == min_score:
			probs.append(solution_count[root_node][i])
		else:
			probs.append(0)

	return get_host_from_count(probs)

def choose_internal_node_host(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
		index = hosts.index(nonterminal.name)

		if not nonterminal.clades[0].is_terminal():
			l_score = score[nonterminal.clades[0]].copy()
			l_count = solution_count[nonterminal.clades[0]].copy()
			l_score[index] -= 1
			for i in range(len(l_score)):
				if l_score[i] != left_score[nonterminal][index]:
					l_count[i] = 0

			nonterminal.clades[0].name = get_host_from_count(l_count)

		if not nonterminal.clades[1].is_terminal():
			r_score = score[nonterminal.clades[1]].copy()
			r_count = solution_count[nonterminal.clades[1]].copy()
			r_score[index] -= 1
			for i in range(len(r_score)):
				if r_score[i] != right_score[nonterminal][index]:
					r_count[i] = 0

			nonterminal.clades[1].name = get_host_from_count(r_count)

def choose_internal_node_host_with_bias(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
		index = hosts.index(nonterminal.name)

		if not nonterminal.clades[0].is_terminal():
			l_score = score[nonterminal.clades[0]].copy()
			l_count = solution_count[nonterminal.clades[0]].copy()
			if min(l_score) == l_score[index]:
				nonterminal.clades[0].name = hosts[index]
			elif min(l_score) + 1 == l_score[index]:
				nonterminal.clades[0].name = hosts[index]
			else:
				l_score[index] -= 1
				for i in range(len(l_score)):
					if l_score[i] != left_score[nonterminal][index]:
						l_count[i] = 0

				nonterminal.clades[0].name = get_host_from_count(l_count)

		if not nonterminal.clades[1].is_terminal():
			r_score = score[nonterminal.clades[1]].copy()
			r_count = solution_count[nonterminal.clades[1]].copy()
			if min(r_score) == r_score[index]:
				nonterminal.clades[1].name = hosts[index]
			elif min(r_score) + 1 == r_score[index]:
				nonterminal.clades[1].name = hosts[index]
			else:
				r_score[index] -= 1
				for i in range(len(r_score)):
					if r_score[i] != right_score[nonterminal][index]:
						r_count[i] = 0

				nonterminal.clades[1].name = get_host_from_count(r_count)

def get_labeled_trees_json_data(rooted_tree):
	# for hosding the labeled trees
	labeled_trees = []
	# for holding the json data
	data = {}
	data['Country of exposure'] = {}
	data['Transmission edges'] = {}
	data['Dated edges'] = []

	# store the dates for all the nodes
	node_date = {}
	exposure = {}
	transmission_edges = []
	for node in rooted_tree.get_terminals():
		node_date[node] = node_name_date[leaf_strain[node]]
		exposure[leaf_strain[node]] = []

	for node in rooted_tree.get_nonterminals():
		node_date[node] = node_name_date[node.name]

	sample_times = 1 if not args.times else args.times
	for i in range(sample_times):
		rooted_tree.root.name = choose_root_host(rooted_tree.root)

		if args.biasedsampling:
			choose_internal_node_host_with_bias(rooted_tree)
		else:
			choose_internal_node_host(rooted_tree)

		labeled_trees.append(copy.deepcopy(rooted_tree))

		# filling the json data
		temp_edges = []
		for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
			for clade in nonterminal.clades:
				if nonterminal.name != clade.name:
					transmission_edge = nonterminal.name + '->' + clade.name
					data['Dated edges'].append([transmission_edge, node_date[nonterminal]])
					temp_edges.append(transmission_edge)
				if clade.is_terminal():
					exposure[leaf_strain[clade]].append(nonterminal.name)

		temp_edges = list(set(temp_edges))
		transmission_edges.extend(temp_edges)

	data['Transmission edges'] = dict(Counter(transmission_edges))
	for strain, countries in exposure.items():
		country_count = dict(Counter(countries))
		data['Country of exposure'][strain] = {'count': country_count, 'country': max(country_count, key=country_count.get)}

	return labeled_trees, data

def create_json_extradata(json_data):
	output_json = args.OUTPUT_FILE + '.json'
	with open(output_json, 'w') as outfile:
		json.dump(json_data, outfile)

def parse_arguments():
	parser = argparse.ArgumentParser(description = 'Process TNet-Geo arguments.')
	parser.add_argument('INPUT_TREE_FILE', action = 'store', type = str, help = 'input tree file name')
	parser.add_argument('OUTPUT_FILE', action = 'store', type = str, help = 'output file name')
	parser.add_argument('-md', '--metadata', default = None, type = str, help = 'csv table with meta data for the tree nodes')
	parser.add_argument('-sd', '--seed', default = None, type = int, help = 'random number generator seed')
	parser.add_argument('-bs', '--biasedsampling', default = False, action = 'store_true', help = 'sample optimal solutions with back transmission bias')
	parser.add_argument('-t', '--times', default = None, type = int, help = 'sample TNet multiple times')
	parser.add_argument('-mx', '--maxprob', default = False, action = 'store_true', help = 'compute highest-probability solution')
	parser.add_argument('-ex', '--extradata', default = False, action = 'store_true', help = 'output a json file with extra output data for further analysis')
	parser.add_argument('-v', '--version', action = 'version', version = 'You are using %(prog)s 1.0')
	return parser.parse_args()

def main():
	# read argparse
	global args
	args = parse_arguments()

	# initialize input_tree, hosts, score, solution_count
	input_tree = initialize_tree(args.INPUT_TREE_FILE)
	if args.metadata:
		read_metadata_initialize_global_variables()

	initialize_leaf_nodes(input_tree)
	initialize_internal_nodes(input_tree)

	# label internal nodes
	np.random.seed(args.seed)
	labeled_trees, json_data = get_labeled_trees_json_data(input_tree)

	# create output files
	Phylo.write(labeled_trees, args.OUTPUT_FILE, 'newick')
	if args.extradata:
		create_json_extradata(json_data)

if __name__ == "__main__":
	sys.setrecursionlimit(100000)
	main()
