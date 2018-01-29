// code by Eli Levy Karin

#ifndef ___SIMULATION_H
#define ___SIMULATION_H

#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <map>
#include "rates.h"
using namespace std;

class simulation
{
public:
	simulation(const vector<size_t> & full_ancestral_state, const double branch_length, const rates & indel_rates, const size_t flank_to_trim) :
		_full_ancestral_state(full_ancestral_state), _indel_rates(indel_rates), _branch_length(branch_length), _flank_to_trim(flank_to_trim)
	{

	} //constructor

	vector<vector<size_t>> get_true_alignment() const;
	vector<vector<size_t>> get_trimmed_true_alignment() const;
	double get_basic_gamma() const { return _indel_rates.get_basic_gamma(); }
	void simulate();


private:
	const vector<size_t> & _full_ancestral_state; // '1' denotes an ancestral character, '2' denotes a descendant character
	const rates & _indel_rates;
	const double _branch_length;
	const size_t _flank_to_trim;

	vector<size_t> _full_descendant_state;
	vector<size_t> _descendant_character_ids;

	map<size_t, vector<size_t>> _event_tracker; // debug - track events

	vector<vector<size_t>> _true_alignment; // '0' denotes a gap, '1' an ancestral character, '2' denotes a descendant character
	vector<vector<size_t>> _trimmed_true_alignment; // '0' denotes a gap, '1' an ancestral character, '2' denotes a descendant character, flanks removed

	void reconstruct_alignment();
	void trim_reconstructed_alignment();

	void get_pos_and_adjust_size_del_event(size_t & rand_left_most_position, size_t & end_deletion_pos, size_t & event_size);

	void get_pos_ins_event(size_t & rand_left_start_position);

};

#endif
