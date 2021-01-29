#!/usr/bin/python3

# Library imports
import argparse, json
from collections import defaultdict
from datetime import date

# global variables
args = None

def create_dated_edges_groups_from_json():
	edge_date = json.load(open(args.INPUT_JSON_FILE))["Dated edges"]
	active_months = []
	for i in range(len(edge_date)):
		parts = edge_date[i][1].split('-')
		dt = date(int(parts[0]), int(parts[1]), int(parts[2]))
		edge_date[i][1] = dt
		active_months.append(date(dt.year, dt.month, 1))

	active_months = list(set(active_months))
	active_months.sort()

	edge_time_count = {}

	for edge, dt in edge_date:
		if edge in edge_time_count:
			edge_time_count[edge][date(dt.year, dt.month, 1)] += 1
		else:
			edge_time_count[edge] = defaultdict(int)
			edge_time_count[edge][date(dt.year, dt.month, 1)] = 1

	result = open(args.edgedategroup, 'w+')
	result.write('edges/dates,{}\n'.format(','.join(str(am.year) + '-' + str(am.month) for am in active_months)))
	for edge, time_count in edge_time_count.items():
		result.write('{},{}\n'.format(edge, ','.join(str(edge_time_count[edge][date(am.year, am.month, 1)]) for am in active_months)))


def parse_arguments():
	parser = argparse.ArgumentParser(description = 'Process TNet-Geo analyse json arguments.')
	parser.add_argument('INPUT_JSON_FILE', action = 'store', type = str, help = 'input json file name')
	parser.add_argument('-eg', '--edgedategroup', default = None, type = str, help = 'csv table with edges grouped into months')
	return parser.parse_args()

def main():
	# read argparse
	global args
	args = parse_arguments()

	# analyze and output file
	if args.edgedategroup:
		create_dated_edges_groups_from_json()

if __name__ == "__main__":
	main()