// code by Eli Levy Karin

#ifndef ___CHOP_TABLE_H
#define ___CHOP_TABLE_H	

#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include <map>
#include <fstream>
#include <string>
#include "simulation.h"
using namespace std;

class chop_table
{
public:

	chop_table(const simulation & simulation_obj, bool should_use_fixed_segment_length) :
		_simulation_obj(simulation_obj), _should_use_fixed_segment_length(should_use_fixed_segment_length)
	{
		fill_N_map();
		fill_L_R_B_maps();
	} //constructor

	map<pair<size_t, size_t>, size_t> get_N_chops_map() const { return _N_chops_map; }
	map<pair<size_t, size_t>, size_t> get_L_chops_map() const { return _L_chops_map; }
	map<pair<size_t, size_t>, size_t> get_R_chops_map() const { return _R_chops_map; }
	map<pair<size_t, size_t>, size_t> get_B_chops_map() const { return _B_chops_map; }
	
	size_t get_total_number_of_N_chops() const { return _total_number_of_N_chops; }
	size_t get_total_number_of_L_chops() const { return _total_number_of_L_chops; }
	size_t get_total_number_of_R_chops() const { return _total_number_of_R_chops; }
	size_t get_total_number_of_B_chops() const { return _total_number_of_B_chops; }

	void write_tables (string table_file_path) const;
	void write_true_trimmed_alignment(string alignment_file_path, size_t number_positions_to_write) const;


private:
	simulation _simulation_obj;
	bool _should_use_fixed_segment_length;
	size_t _max_segment_size = 200; // fixed

	map<pair<size_t, size_t>, size_t> _N_chops_map;
	map<pair<size_t, size_t>, size_t> _L_chops_map;
	map<pair<size_t, size_t>, size_t> _R_chops_map;
	map<pair<size_t, size_t>, size_t> _B_chops_map;

	vector<double> _anc_length_cum_probs; // will be computed according to a geometric distribution.

	size_t _total_number_of_N_chops;
	size_t _total_number_of_L_chops;
	size_t _total_number_of_R_chops;
	size_t _total_number_of_B_chops;

	void fill_N_map();
	void fill_L_R_B_maps();

	void compute_cumu_dist_of_ancestral_lengths();
	size_t sample_segment_length();

};

#endif

