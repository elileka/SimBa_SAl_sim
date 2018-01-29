// code by Eli Levy Karin

#include "chop_table.h"

void chop_table::write_tables (string table_file_path) const
{
	ofstream myfile;
	myfile.open(table_file_path);

	myfile << "type_chop\ti\tj\tnum_times_chop_observed\ttotal_num_chops_of_this_type\tchop_probability" << endl;

	for (auto const& chop : _N_chops_map)
	{
		
		pair<size_t, size_t> curr_ij = chop.first;
		size_t count = chop.second;
		double chop_prob = (double)count / _total_number_of_N_chops;
		myfile << "N\t" << curr_ij.first << "\t" << curr_ij.second << "\t" << count << "\t" << _total_number_of_N_chops << "\t" << chop_prob << endl;
	}

	for (auto const& chop : _L_chops_map)
	{

		pair<size_t, size_t> curr_ij = chop.first;
		size_t count = chop.second;
		double chop_prob = (double)count / _total_number_of_L_chops;
		myfile << "L\t" << curr_ij.first << "\t" << curr_ij.second << "\t" << count << "\t" << _total_number_of_L_chops << "\t" << chop_prob << endl;
	}

	for (auto const& chop : _R_chops_map)
	{

		pair<size_t, size_t> curr_ij = chop.first;
		size_t count = chop.second;
		double chop_prob = (double)count / _total_number_of_R_chops;
		myfile << "R\t" << curr_ij.first << "\t" << curr_ij.second << "\t" << count << "\t" << _total_number_of_R_chops << "\t" << chop_prob << endl;
	}

	for (auto const& chop : _B_chops_map)
	{

		pair<size_t, size_t> curr_ij = chop.first;
		size_t count = chop.second;
		double chop_prob = (double)count / _total_number_of_B_chops;
		myfile << "B\t" << curr_ij.first << "\t" << curr_ij.second << "\t" << count << "\t" << _total_number_of_B_chops << "\t" << chop_prob << endl;
	}

	myfile.close();
}

void chop_table::write_true_trimmed_alignment(string alignment_file_path, size_t number_positions_to_write) const
{
	vector<vector<size_t>> timmed_true_simulated_alignment = _simulation_obj.get_trimmed_true_alignment();
	size_t length_of_trimmed_alignment = timmed_true_simulated_alignment.size();
	if (number_positions_to_write > length_of_trimmed_alignment)
	{
		number_positions_to_write = length_of_trimmed_alignment;
	}

	ofstream myfile;
	myfile.open(alignment_file_path);

	myfile << ">anc" << endl;

	for (size_t curr_alignment_ind = 0; curr_alignment_ind < number_positions_to_write; curr_alignment_ind++)
	{
		size_t anc_char = timmed_true_simulated_alignment[curr_alignment_ind][0];
		myfile << anc_char;
	}
	myfile << endl;

	myfile << ">des" << endl;
	for (size_t curr_alignment_ind = 0; curr_alignment_ind < number_positions_to_write; curr_alignment_ind++)
	{
		size_t des_char = timmed_true_simulated_alignment[curr_alignment_ind][1];
		myfile << des_char;
	}
	myfile << endl;

	myfile.close();
}

void chop_table::fill_N_map()
{
	vector<vector<size_t>> timmed_true_simulated_alignment = _simulation_obj.get_trimmed_true_alignment();

	size_t curr_chop_i = 0; // number deleted from ancestor 
	size_t curr_chop_j = 0; // number inserted in descendant

	_total_number_of_N_chops = 0;

	bool left_flag = true;

	size_t alignment_ind = 0;

	// '0' denotes a gap, '1' an ancestral character, '2' denotes a descendant character
	for (size_t alignment_ind = 0; alignment_ind < timmed_true_simulated_alignment.size(); alignment_ind++)
	{
		if ((timmed_true_simulated_alignment[alignment_ind][0] != 0) && (timmed_true_simulated_alignment[alignment_ind][1] != 0)) // match
		{
			// if this is the first match from the left - disregard it
			if (left_flag)
			{
				left_flag = false;
			}
			else // this is not the first time we encounter a match - it means an N chop ended
			{
				_total_number_of_N_chops++;
				pair<size_t, size_t> curr_ij = make_pair(curr_chop_i, curr_chop_j);
				if (_N_chops_map.find(curr_ij) == _N_chops_map.end())
				{
					// not in map
					_N_chops_map[curr_ij] = 1;
				}
				else
				{
					// exists in map
					_N_chops_map[curr_ij] = _N_chops_map[curr_ij] + 1;
				}
			}
			// a chop ended - restart the counts
			curr_chop_i = 0; // number deleted from ancestor 
			curr_chop_j = 0; // number inserted in descendant
		}
		else if (timmed_true_simulated_alignment[alignment_ind][0] != 0) // only at the ancestor
		{
			curr_chop_i++;
		}
		else // only at the descendant
		{
			curr_chop_j++;
		}
	}
}

void chop_table::fill_L_R_B_maps()
{
	// the trimmed alignment is cut to segments according to a geometric distribution of the ancestral length. 
	// The edges of each segment are used to estimate the left, right and both chops

	vector<vector<size_t>> timmed_true_simulated_alignment = _simulation_obj.get_trimmed_true_alignment();
	size_t length_of_trimmed_alignment = timmed_true_simulated_alignment.size();
	_total_number_of_L_chops = 0;
	_total_number_of_R_chops = 0;
	_total_number_of_B_chops = 0;

	// fill in the cumulative distribution of ancestral lengths
	compute_cumu_dist_of_ancestral_lengths();

	size_t curr_alignment_ind = 0;
	while (curr_alignment_ind < length_of_trimmed_alignment)
	{
		size_t curr_segment_length = sample_segment_length();
		size_t curr_boundary = curr_alignment_ind + curr_segment_length;
		if (curr_boundary > length_of_trimmed_alignment)
		{
			curr_boundary = length_of_trimmed_alignment;
		}

		bool left_flag = true;
		bool not_event_a_single_match = true; // if we see a match this becomes 'false'
		pair<size_t, size_t> left_chop;
		pair<size_t, size_t> right_chop;

		size_t curr_chop_i = 0; // number deleted from ancestor 
		size_t curr_chop_j = 0; // number inserted in descendant

		// '0' denotes a gap, '1' an ancestral character, '2' denotes a descendant character
		while (curr_alignment_ind < curr_boundary)
		{
			if ((timmed_true_simulated_alignment[curr_alignment_ind][0] != 0) && (timmed_true_simulated_alignment[curr_alignment_ind][1] != 0)) // match
			{
				not_event_a_single_match = false; // we saw a match!
				// handle left chop
				if (left_flag)
				{
					left_chop = make_pair(curr_chop_i, curr_chop_j);

					_total_number_of_L_chops++;
					if (_L_chops_map.find(left_chop) == _L_chops_map.end())
					{
						// not in map
						_L_chops_map[left_chop] = 1;
					}
					else
					{
						// exists in map
						_L_chops_map[left_chop] = _L_chops_map[left_chop] + 1;
					}
				}

				// restart
				left_flag = false; // we only look at the start of the segment
				curr_chop_i = 0; // deleted from ancestor 
				curr_chop_j = 0; // inserted in descendant
			}
			else if (timmed_true_simulated_alignment[curr_alignment_ind][0] != 0) // only at the ancestor
			{
				curr_chop_i++;
			}
			else // only at the descendant
			{
				curr_chop_j++;
			}

			// handle right chop - will be overridden until the end of the inner loop
			right_chop = make_pair(curr_chop_i, curr_chop_j);

			curr_alignment_ind++;
		}

		// finished going over the segment
		if (not_event_a_single_match) // B chop!!!
		{
			_total_number_of_B_chops++;
			if (_B_chops_map.find(right_chop) == _B_chops_map.end())
			{
				// not in map
				_B_chops_map[right_chop] = 1;
			}
			else
			{
				// exists in map
				_B_chops_map[right_chop] = _B_chops_map[right_chop] + 1;
			}
		}
		else // R chop!!! - we saw a match on the left
		{
			_total_number_of_R_chops++;
			if (_R_chops_map.find(right_chop) == _R_chops_map.end())
			{
				// not in map
				_R_chops_map[right_chop] = 1;
			}
			else
			{
				// exists in map
				_R_chops_map[right_chop] = _R_chops_map[right_chop] + 1;
			}
		}
	}
}

void chop_table::compute_cumu_dist_of_ancestral_lengths()
{
	double basic_gamma = _simulation_obj.get_basic_gamma();

	double sum_of_small_probs = 0.0;
	for (size_t i = 1; i < _max_segment_size; i++)
	{
		double curr_prob = pow(basic_gamma, i) * (1 - basic_gamma);
		sum_of_small_probs = sum_of_small_probs + curr_prob;
		_anc_length_cum_probs.push_back(sum_of_small_probs); // cumulative probability
	}
	_anc_length_cum_probs.push_back(1); // complete to 1
}

size_t chop_table::sample_segment_length()
{
	
	if (_should_use_fixed_segment_length)
	{
		return (_max_segment_size);
	}

	random_device rd;  // will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
	uniform_real_distribution<double> distribution(0.0, 1);
	double rand_in_range = distribution(gen);

	for (size_t i = 0; i < _anc_length_cum_probs.size(); i++)
	{
		size_t curr_length = (i + 1);
		if (rand_in_range <= _anc_length_cum_probs[i])
		{
			return (curr_length);
		}
	}

	// should not get here:
	return (0);
}
