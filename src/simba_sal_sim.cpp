// code by Eli Levy Karin

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <map>
#include <cmath>
#include "rates.h"
#include "simulation.h"
#include "chop_table.h"

using namespace std;

void check_reversibility_eq_30(map<pair<size_t, size_t>, size_t> chops_map, size_t total_number_of_chops, double basic_gamma, ofstream & myfilestream)
{
	// check reversibility:
	// for all items in curr map, we should get:
	// gamma^(i+1) * curr_ij = gamma^(j+1) * curr_ji
	myfilestream << "a total of " << total_number_of_chops << " were examined" << endl;

	for (auto const& chop : chops_map)
	{
		pair<size_t, size_t> curr_ij = chop.first;

		size_t i = curr_ij.first;
		size_t j = curr_ij.second;

		double prob_ij = (double)chops_map[curr_ij] / total_number_of_chops;
		double left_side = pow(basic_gamma, (i + 1)) * prob_ij;

		pair<size_t, size_t> curr_ji = make_pair(j, i);
		if (chops_map.find(curr_ji) == chops_map.end())
		{
			// curr_ji not in map
			myfilestream << "chop_" << j << "," << i << " not in map; prob_" << i << "," << j << " = " << prob_ij << "; left_side = " << left_side << endl;
		}
		else // curr_ji exists in map
		{
			double prob_ji = (double)chops_map[curr_ji] / total_number_of_chops;
			double right_side = pow(basic_gamma, (j + 1)) * prob_ji;

			myfilestream << "both in map; prob_" << i << "," << j << " = " << prob_ij << "; left_side = " << left_side << "; prob_" << j << "," << i << " = " << prob_ji << "; right_side = " << right_side << endl;
		}

	}

}

void check_reversibility_eq_31(map<pair<size_t, size_t>, size_t> chops_map_L, size_t total_number_of_L_chops, map<pair<size_t, size_t>, size_t> chops_map_B, size_t total_number_of_B_chops, double basic_gamma, ofstream & myfilestream)
{
	// check reversibility:
	// for all items in curr map, we should get:
	// sum: Lij * gamma^(i+1) + sum: Bij * gamma^i * (1 - gamma) = 1
	double sum_L = 0.0;
	double sum_of_L_probs = 0.0;

	for (auto const& chop : chops_map_L)
	{
		pair<size_t, size_t> curr_ij = chop.first;
		size_t i = curr_ij.first;

		double prob_ij = (double)chops_map_L[curr_ij] / total_number_of_L_chops;
		sum_L = sum_L + (prob_ij * pow(basic_gamma, (i + 1)));

		sum_of_L_probs = sum_of_L_probs + prob_ij;
	}

	double sum_B = 0.0;
	double sum_of_B_probs = 0.0;

	for (auto const& chop : chops_map_B)
	{
		pair<size_t, size_t> curr_ij = chop.first;
		size_t i = curr_ij.first;

		double prob_ij = (double)chops_map_B[curr_ij] / total_number_of_B_chops;
		sum_B = sum_B + (prob_ij * pow(basic_gamma, i) * (1 - basic_gamma));

		sum_of_B_probs = sum_of_B_probs + prob_ij;
	}

	double sum = sum_L + sum_B;

	myfilestream << "sum: Lij * gamma^(i+1) + sum: Bij * gamma^i * (1 - gamma) = " << sum << endl;
	myfilestream << "sum: Lij probs = " << sum_of_L_probs << endl;
	myfilestream << "sum: Bij probs = " << sum_of_B_probs << endl;
}

void check_reversibility_eq_32(map<pair<size_t, size_t>, size_t> chops_map_N, map<pair<size_t, size_t>, size_t> chops_map_R, size_t total_number_of_N_chops, size_t total_number_of_R_chops, double basic_gamma, ofstream & myfilestream)
{
	// check reversibility:
	// for all items in curr map, we should get:
	// sum: Nij * gamma^(i+1) + sum: Rij * gamma^i * (1 - gamma) = 1
	double sum_N = 0.0;
	double sum_of_N_probs = 0.0;

	for (auto const& chop : chops_map_N)
	{
		pair<size_t, size_t> curr_ij = chop.first;
		size_t i = curr_ij.first;

		double prob_ij = (double)chops_map_N[curr_ij] / total_number_of_N_chops;
		sum_N = sum_N + (prob_ij * pow(basic_gamma, (i + 1)));

		sum_of_N_probs = sum_of_N_probs + prob_ij;
	}

	double sum_R = 0.0;
	double sum_of_R_probs = 0.0;

	for (auto const& chop : chops_map_R)
	{
		pair<size_t, size_t> curr_ij = chop.first;
		size_t i = curr_ij.first;

		double prob_ij = (double)chops_map_R[curr_ij] / total_number_of_R_chops;
		sum_R = sum_R + (prob_ij * pow(basic_gamma, i) * (1 - basic_gamma));

		sum_of_R_probs = sum_of_R_probs + prob_ij;
	}

	double sum = sum_N + sum_R;

	myfilestream << "sum: Nij * gamma^(i+1) + sum: Rij * gamma^i * (1 - gamma) = " << sum << endl;
	myfilestream << "sum: Nij probs = " << sum_of_N_probs << endl;
	myfilestream << "sum: Rij probs = " << sum_of_R_probs << endl;
}

void check_reversibility_eq_33_N(map<pair<size_t, size_t>, size_t> chops_map_N, size_t total_number_of_N_chops, double sum_mus, ofstream & myfilestream)
{
	// check reversibility:
	// for all items in curr map, we should get:
	// sum: N0i = exp (- sum(mu_k))
	double sum = 0.0;

	for (auto const& chop : chops_map_N)
	{
		pair<size_t, size_t> curr_ij = chop.first;
		size_t i = curr_ij.first;
		if (i != 0)
		{
			continue;
		}

		double prob_0i = (double)chops_map_N[curr_ij] / total_number_of_N_chops;
		sum = sum + prob_0i;
	}

	double right_side = exp(-sum_mus);

	myfilestream << "sum: N0i = " << sum << "; exp (- sum(mu_k)) = " << right_side << endl;
}

void check_reversibility_eq_33_L(map<pair<size_t, size_t>, size_t> chops_map_L, size_t total_number_of_L_chops, double sum_k_mus, ofstream & myfilestream)
{
	// check reversibility:
	// for all items in curr map, we should get:
	// sum: L0i = exp (- sum(k * mu_k))
	double sum = 0.0;

	for (auto const& chop : chops_map_L)
	{
		pair<size_t, size_t> curr_ij = chop.first;
		size_t i = curr_ij.first;

		if (i != 0)
		{
			continue;
		}

		double prob_0i = (double)chops_map_L[curr_ij] / total_number_of_L_chops;
		sum = sum + prob_0i;
	}

	double right_side = exp(-sum_k_mus);

	myfilestream << "sum: L0i = " << sum << "; exp (- sum(k * mu_k)) = " << right_side << endl;
}

string getCmdOption(int num_args, const char* argv[], const std::string& option)
{
	string val;
	for (int i = 0; i < num_args; ++i)
	{
		string arg = argv[i];
		size_t start_position_of_param_name = arg.find(option); // found is the start position of the match
		if (start_position_of_param_name == 0)
		{
			size_t start_pos_of_value = start_position_of_param_name + option.size();
			val = arg.substr(start_pos_of_value);
			return val;
		}
	}
	return val;
}

int main(int argc, const char * argv[])
{
	if (argc < 2)
	{
		cerr << "The program takes parameter values using '=' pairs. The order does not matter." << endl;
		cerr << "The parameter names are: " << endl;
		cerr << "type_model (long-indel / power-law inspired: LI/PL, defaults to 'LI')" << endl;
		cerr << "r_param (controls tendency for longer indel under the long indel model, range: [0,1))" << endl;
		cerr << "basic_mu (controls the total per-position deletion rate under the long indel model, range: [0,1))" << endl;
		cerr << "length_of_t (length of branch connecting the ancestor to descendant)" << endl;
		cerr << "out_files_prefix (where to write the output, parameter values will be added to out_files_prefix)" << endl;
		cerr << "number_postiotions_to_write (from the true trimmed PWA, defaults to 5000, won't exceed the length of the trimmed PWA)" << endl;
		cerr << "A (controls tendency for longer indel under the power-law inspired model, needed only if 'PL', range: (1,))" << endl;
		cerr << "IR (controls the ratio of indel to substitution events, needed only if 'PL', range: [0,1])" << endl;

		cerr << "Example usage: " << endl << argv[0] << " r_param=0.3 basic_mu=0.1 length_of_t=1 out_files_prefix=PATH/SimBa-SAl_sim/" << endl;
		return 1;
	}

	// collection of user parameter values (some may be missing)
	string user_r_param = getCmdOption(argc, argv, "r_param=");
	string user_basic_mu = getCmdOption(argc, argv, "basic_mu=");
	string user_A = getCmdOption(argc, argv, "A=");
	string user_IR = getCmdOption(argc, argv, "IR=");
	string user_length_of_t = getCmdOption(argc, argv, "length_of_t=");
	string user_out_files_prefix = getCmdOption(argc, argv, "out_files_prefix=");
	string user_number_postiotions_to_write = getCmdOption(argc, argv, "number_postiotions_to_write=");
	string user_type_model = getCmdOption(argc, argv, "type_model=");

	// parameter values as strings:
	string r_param_str = user_r_param;
	string basic_mu_str = user_basic_mu;
	string A_str = user_A;
	string IR_str = user_IR;
	string length_of_t_str = user_length_of_t;
	string out_files_prefix = user_out_files_prefix;

	// special treatment of parameters with default values:
	string type_model;
	if (user_type_model == "PL")
	{
		type_model = "PL";
		if ((A_str == "") || (IR_str == ""))
		{
			cerr << "PL model chosen but A and/or IR are not defined" << endl;
			return 1;
		}
	}
	else
	{
		type_model = "LI";
		if ((r_param_str == "") || (basic_mu_str == ""))
		{
			cerr << "LI model chosen but r_param and/or basic_mu are not defined" << endl;
			return 1;
		}
	}
	string number_postiotions_to_write_str;
	if (user_number_postiotions_to_write != "")
	{
		number_postiotions_to_write_str = user_number_postiotions_to_write;
	}
	else
	{
		number_postiotions_to_write_str = "5000";
	}


	// some parameters we fix for all runs - can be changed in the future, if needed:
	size_t length_of_anc = 1000000;
	double basic_gamma = 0.99; // basic_gamma = lambda_1 / mu_1 (estimated from unaligned seqs)
	size_t max_indel_length = 50;
	size_t flank_to_trim = 200;  // this will be removed from either end to avoid edge problems

	bool is_long_indel = true; // long indel model (by N formulae)
	double r_param; // long-indel param - [0,1). If 0 - TKF91
	double basic_mu; // long-indel param - equal to the sum of k*mu_k for all k
	double A; // power law inspired shape param
	double IR; // per position insertion = deletion rate (relative to sub.) 

	if (type_model == "PL")
	{
		is_long_indel = false; // power law
		A = atof(A_str.c_str());
		IR = atof(IR_str.c_str());
	}
	else
	{
		r_param = atof(r_param_str.c_str()); // [0,1). If 0 - TKF91
		basic_mu = atof(basic_mu_str.c_str()); // equal to the sum of k*mu_k for all k
	}

	double length_of_t = atof(length_of_t_str.c_str()); // branch connecting the anc and des, measured in expected number of substitution per position
	size_t number_postiotions_to_write = (size_t)atoi(number_postiotions_to_write_str.c_str()); // the whole alignment is very long - avoid big files...

	// prepare file names prefix
	stringstream ss_pref;
	ss_pref << out_files_prefix;

	// prepare file names suffix
	stringstream ss_suff;
	if (is_long_indel)
	{
		ss_suff << "_r_param_" << r_param << "_basic_mu_" << basic_mu << "_basic_gamma_" << basic_gamma << "_max_indel_length_" << max_indel_length << "_length_of_anc_" << length_of_anc << "_t_" << length_of_t << ".txt";
	}
	else
	{
		ss_suff << "_A_param_" << A << "_IR_" << IR << "_max_indel_length_" << max_indel_length << "_length_of_anc_" << length_of_anc << "_t_" << length_of_t << ".txt";
	}
	
	// prepare tables file name
	stringstream ss_table;
	ss_table << ss_pref.str() << "chop_tables" << ss_suff.str();
	std::string chop_tables_file_string = ss_table.str();
	cout << "chop tables will be written to: " << chop_tables_file_string << endl;

	// prepare partial alignment file name
	stringstream ss_true_trimmed_alignment;
	ss_true_trimmed_alignment << ss_pref.str() << "true_alignment_fasta_" << number_postiotions_to_write << "_positions" << ss_suff.str();
	std::string ss_true_trimmed_alignment_file_string = ss_true_trimmed_alignment.str();

	// prepare sanity checks file name
	stringstream ss_report;
	ss_report << ss_pref.str() << "report_consistency_checks_reversibility" << ss_suff.str();
	std::string report_file_string = ss_report.str();

	// set ancestral vector
	vector<size_t> anc_state;
	for (size_t i = 0; i < length_of_anc; i++)
	{
		anc_state.push_back(1);
	}

	// let's get started!
	cout << "starting to simulate " << ss_suff.str() << endl;
	
	if (is_long_indel)
	{
		// construct rates object
		rates rates_for_long_indel(basic_mu, basic_gamma, r_param, max_indel_length);

		simulation test_sim(anc_state, length_of_t, rates_for_long_indel, flank_to_trim);
		test_sim.simulate();

		cout << "starting to compute maps" << endl;
		bool should_use_fixed_segment_length = false;
		chop_table chop_table_obj(test_sim, should_use_fixed_segment_length);

		// write tables to file
		chop_table_obj.write_tables(chop_tables_file_string);

		// write out part of the true alignment
		chop_table_obj.write_true_trimmed_alignment(ss_true_trimmed_alignment_file_string, number_postiotions_to_write);

		// perform consistency checks according to the long indel model paper, appendix B
		map<pair<size_t, size_t>, size_t> N_chops_map = chop_table_obj.get_N_chops_map();
		size_t total_number_of_N_chops = chop_table_obj.get_total_number_of_N_chops();

		map<pair<size_t, size_t>, size_t> L_chops_map = chop_table_obj.get_L_chops_map();
		size_t total_number_of_L_chops = chop_table_obj.get_total_number_of_L_chops();

		map<pair<size_t, size_t>, size_t> R_chops_map = chop_table_obj.get_R_chops_map();
		size_t total_number_of_R_chops = chop_table_obj.get_total_number_of_R_chops();

		map<pair<size_t, size_t>, size_t> B_chops_map = chop_table_obj.get_B_chops_map();
		size_t total_number_of_B_chops = chop_table_obj.get_total_number_of_B_chops();

		// write the consistency checks to file:
		ofstream myfile;
		myfile.open(report_file_string);

		myfile << "checking N (eq 30)" << endl;
		check_reversibility_eq_30(N_chops_map, total_number_of_N_chops, basic_gamma, myfile);
		myfile << endl;
		myfile << "checking L (eq 30)" << endl;
		check_reversibility_eq_30(L_chops_map, total_number_of_L_chops, basic_gamma, myfile);
		myfile << endl;
		myfile << "checking R (eq 30)" << endl;
		check_reversibility_eq_30(R_chops_map, total_number_of_R_chops, basic_gamma, myfile);
		myfile << endl;
		myfile << "checking B (eq 30)" << endl;
		check_reversibility_eq_30(B_chops_map, total_number_of_B_chops, basic_gamma, myfile);
		myfile << endl;
		myfile << "checking L + B (eq 31)" << endl;
		myfile << endl;
		check_reversibility_eq_31(L_chops_map, total_number_of_L_chops, B_chops_map, total_number_of_B_chops, basic_gamma, myfile);
		myfile << endl;
		myfile << "checking N + R (eq 32)" << endl;
		myfile << endl;
		check_reversibility_eq_32(N_chops_map, R_chops_map, total_number_of_N_chops, total_number_of_R_chops, basic_gamma, myfile);
		myfile << endl;
		
		double sum_k_mus = 0.0;
		double sum_mus = 0.0;
		vector<double> mus = rates_for_long_indel.get_mus();
		double sum_k_lambdas = 0.0;
		vector<double> lambdas = rates_for_long_indel.get_lambdas();
		myfile << "here are the deletion rates:" << endl;
		for (size_t ind = 0; ind < mus.size(); ind++)
		{
			size_t k = ind + 1;
			sum_k_mus = sum_k_mus + (k * mus[ind]);
			sum_mus = sum_mus + mus[ind];
			myfile << "mu_" << k << " = " << mus[ind] << ", ";
		}
		myfile << endl << "here are the insertion rates:" << endl;
		for (size_t ind = 0; ind < lambdas.size(); ind++)
		{
			size_t k = ind + 1;
			sum_k_lambdas = sum_k_lambdas + (k * lambdas[ind]);
			myfile << "lambda_" << k << " = " << lambdas[ind] << ", ";
		}
		myfile << endl << "sum_k_mus is: " << sum_k_mus << endl;
		myfile << "sum_k_lambdas is: " << sum_k_lambdas << endl;

		myfile << endl;
		myfile << "checking N vs. sum mus (eq 33)" << endl;
		myfile << endl;
		check_reversibility_eq_33_N(N_chops_map, total_number_of_N_chops, sum_mus, myfile);
		myfile << endl;
		myfile << "checking L vs. sum k*mus (eq 33)" << endl;
		myfile << endl;
		check_reversibility_eq_33_L(L_chops_map, total_number_of_L_chops, sum_k_mus, myfile);
		myfile << endl;

		myfile.close();
	}

	else // INDELible power law
	{
		// construct rates object
		rates rates_for_INDELible_PL(IR, A, max_indel_length);

		simulation test_sim(anc_state, length_of_t, rates_for_INDELible_PL, flank_to_trim);
		test_sim.simulate();

		cout << "starting to compute maps" << endl;
		bool should_use_fixed_segment_length = true;
		chop_table chop_table_obj(test_sim, should_use_fixed_segment_length);

		// write tables to file
		chop_table_obj.write_tables(chop_tables_file_string);

		// write out part of the true alignment
		chop_table_obj.write_true_trimmed_alignment(ss_true_trimmed_alignment_file_string, number_postiotions_to_write);

		ofstream myfile;
		myfile.open(report_file_string);
		double sum_k_mus = 0.0;
		vector<double> mus = rates_for_INDELible_PL.get_mus();
		double sum_k_lambdas = 0.0;
		vector<double> lambdas = rates_for_INDELible_PL.get_lambdas();
		myfile << "here are the deletion rates:" << endl;
		for (size_t ind = 0; ind < mus.size(); ind++)
		{
			size_t k = ind + 1;
			sum_k_mus = sum_k_mus + (k * mus[ind]);
			myfile << "mu_" << k << " = " << mus[ind] << ", ";
		}
		myfile << endl << "here are the insertion rates:" << endl;
		for (size_t ind = 0; ind < lambdas.size(); ind++)
		{
			size_t k = ind + 1;
			sum_k_lambdas = sum_k_lambdas + (k * lambdas[ind]);
			myfile << "lambda_" << k << " = " << lambdas[ind] << ", ";
		}
		myfile << endl << "sum_k_mus is: " << sum_k_mus << endl;
		myfile << "sum_k_lambdas is: " << sum_k_lambdas << endl;

		myfile.close();


	}

	return 0;

}


