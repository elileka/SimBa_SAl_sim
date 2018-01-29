// code by Eli Levy Karin

#include "simulation.h"

vector<vector<size_t>> simulation::get_true_alignment() const
{
	return _true_alignment;
}

vector<vector<size_t>> simulation::get_trimmed_true_alignment() const
{
	return _trimmed_true_alignment;
}

void simulation::simulate()
{
	// simulation starts with descendant as ancestor
	_full_descendant_state = _full_ancestral_state; // this creates a new copy
	
	// initialize the character ids
	for (size_t i = 0; i < _full_descendant_state.size(); i++)
	{
		_descendant_character_ids.push_back(i);
	}
	size_t global_insertion_counter = _full_descendant_state.size();

	size_t curr_descendant_length = _full_descendant_state.size();

	double total_time = 0.0;

	size_t event_id = 0; // debug
	size_t proxy_for_num_del_chars = 0; // debug
	size_t proxy_for_num_ins_chars = 0; // debug
	size_t num_pos_inserted_then_deleted = 0;

	while (1)
	{
		// first - take care of edge case (curr_descendant_length is at 0 length)
		if (curr_descendant_length == 0)
		{
			// we will avoid a sink by smapling an insertion event for the "immortal" link:
			double time_to_next_insertion_jump = _indel_rates.sample_insertion_jump_time();
			total_time = total_time + time_to_next_insertion_jump;
			// if the next jump exceeds the branch - we are done
			if (total_time > _branch_length)
			{
				// we remain with the last state - nothing changes
				break;
			}
			event_id++; // debug

			size_t event_size = 0;
			_indel_rates.sample_insertion_size(event_size);
			for (size_t i = 1; i <= event_size; i++)
			{
				_full_descendant_state.push_back(2);  // 2 denotes a descendent character
				_descendant_character_ids.push_back(global_insertion_counter);
				global_insertion_counter++;
			}
		}
		// done dealing with the edge case - below is the actual code

		// sample the time to the next jump
		double time_to_next_jump = _indel_rates.sample_jump_time(curr_descendant_length);
		
		total_time = total_time + time_to_next_jump;
		// if the next jump exceeds the branch - we are done
		if (total_time > _branch_length)
		{
			// we remain with the last state - nothing changes
			break;
		}
		else
		{
			event_id++; // debug

			// sample event according to relative rates
			size_t event_type = 0;
			size_t event_size = 0;

			_indel_rates.sample_event_type(event_type, event_size);
			size_t rand_left_start_position; // will be sampled with a helper function

			// events happen to the right
			if (event_type == 0) // deletion
			{
				proxy_for_num_del_chars += event_size; // debug

				size_t end_deletion_pos;

				// call helper function to sample:
				get_pos_and_adjust_size_del_event(rand_left_start_position, end_deletion_pos, event_size);
				
				// for debug - track the event:
				for (size_t i = rand_left_start_position; i <= end_deletion_pos; i++)
				{
					size_t curr_char_id = _descendant_character_ids[i];
					if (_event_tracker.find(curr_char_id) == _event_tracker.end())
					{
						// not in map
						vector<size_t> pos_events;
						pos_events.push_back(event_id);
						_event_tracker[curr_char_id] = pos_events;
					}
					else
					{
						// exists in map
						_event_tracker[curr_char_id].push_back(event_id);
						num_pos_inserted_then_deleted++;
					}
				}
				// end track the event

				// remove the characters and their ids:
				_full_descendant_state.erase(_full_descendant_state.begin() + rand_left_start_position, _full_descendant_state.begin() + end_deletion_pos);
				_descendant_character_ids.erase(_descendant_character_ids.begin() + rand_left_start_position, _descendant_character_ids.begin() + end_deletion_pos);
			}
			if (event_type == 1) // insertion
			{
				proxy_for_num_ins_chars += event_size; // debug

				// call helper function to sample:
				get_pos_ins_event(rand_left_start_position);
				for (size_t i = 1; i <= event_size; i++)
				{
					_full_descendant_state.insert(_full_descendant_state.begin() + rand_left_start_position + i, 2); // 2 denotes a descendent character
					_descendant_character_ids.insert(_descendant_character_ids.begin() + rand_left_start_position + i, global_insertion_counter);
					global_insertion_counter++;
				}

				// for debug - track the event:
				for (size_t i = (rand_left_start_position + 1); i <= event_size; i++)
				{
					size_t curr_char_id = _descendant_character_ids[i];
					// inserted position is guaranteed not to be in map
					vector<size_t> pos_events;
					pos_events.push_back(event_id);
					_event_tracker[curr_char_id] = pos_events;
				}
				// end track the event
			}
		}
	}

	cout << "simulation had " << event_id << " events." << endl;
	cout << "proxy_for_num_del_chars is " << proxy_for_num_del_chars << endl;
	cout << "proxy_for_num_ins_chars is " << proxy_for_num_ins_chars << endl;
	cout << "number of positions inserted then deleted is " << num_pos_inserted_then_deleted << endl;

	reconstruct_alignment();
}

void simulation::reconstruct_alignment()
{
	// this method is called at the end of the simulation
	size_t num_ancestral_chars = _full_ancestral_state.size();
	size_t num_descendant_chars = _full_descendant_state.size(); // after the simulation this is full

	size_t anc_ind = 0; // this is also equal to the character id
	size_t des_ind = 0; // this might not be equal to the character id

	while ((anc_ind < num_ancestral_chars) && (des_ind < num_descendant_chars)) // while anc and des vectors are not empty
	{
		size_t des_id = _descendant_character_ids[des_ind];
		if (des_id >= num_ancestral_chars) // the character is a new des character
		{
			vector<size_t> curr_pos;
			curr_pos.push_back(0); // 0 denotes a gap character
			curr_pos.push_back(_full_descendant_state[des_ind]);
			_true_alignment.push_back(curr_pos);

			des_ind++;
		}
		else
		{
			if (des_id == anc_ind) // homology
			{
				vector<size_t> curr_pos;
				curr_pos.push_back(_full_ancestral_state[anc_ind]);
				curr_pos.push_back(_full_descendant_state[des_ind]);
				_true_alignment.push_back(curr_pos);

				anc_ind++;
				des_ind++;
			}
			else // des_id > anc_ind --> a deleted ancestral character
			{
				vector<size_t> curr_pos;
				curr_pos.push_back(_full_ancestral_state[anc_ind]);
				curr_pos.push_back(0);
				_true_alignment.push_back(curr_pos);

				anc_ind++;
			}
		}
	}

	// it could be that one of the vectors is not empty
	while (anc_ind < num_ancestral_chars) // some ancestral characters are left
	{
		vector<size_t> curr_pos;
		curr_pos.push_back(_full_ancestral_state[anc_ind]);
		curr_pos.push_back(0); // 0 denotes a gap character
		_true_alignment.push_back(curr_pos);

		anc_ind++;
	}

	while (des_ind < num_descendant_chars)
	{
		vector<size_t> curr_pos;
		curr_pos.push_back(0); // 0 denotes a gap character
		curr_pos.push_back(_full_descendant_state[des_ind]);
		_true_alignment.push_back(curr_pos);

		des_ind++;
	}

	trim_reconstructed_alignment();
}

void simulation::trim_reconstructed_alignment()
{
	// if we get here, we already have the reconstructed true alignment

	if (_true_alignment.size() < 2 * _flank_to_trim + 1)
	{
		cerr << "true alignment is shorter than the flanks to remove - there is a problem" << endl;
		exit(1);
	}
	
	for (size_t i = _flank_to_trim; i < (_true_alignment.size() - _flank_to_trim); i++)
	{
		_trimmed_true_alignment.push_back(_true_alignment[i]);
	}
}


void simulation::get_pos_and_adjust_size_del_event(size_t & rand_left_start_position, size_t & end_deletion_pos, size_t & event_size)
{
	size_t curr_descendant_length = _full_descendant_state.size();
	random_device rd;  // will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
	uniform_int_distribution<size_t> dis(0, (curr_descendant_length - 1));
	rand_left_start_position = dis(gen);

	// trim event if it deletes characters beyond the right border:
	end_deletion_pos = rand_left_start_position + event_size;
	if (end_deletion_pos > curr_descendant_length)
	{
		end_deletion_pos = curr_descendant_length;
	}
}


void simulation::get_pos_ins_event(size_t & rand_left_start_position)
{
	size_t curr_descendant_length = _full_descendant_state.size();
	random_device rd;  // will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
	uniform_int_distribution<size_t> dis(0, (curr_descendant_length - 1));
	rand_left_start_position = dis(gen);
}
