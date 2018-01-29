// code by Eli Levy Karin

#ifndef ___RATES_H
#define ___RATES_H

#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <iostream>
using namespace std;

class rates
{
public:
	
	rates(double basic_mu, double basic_gamma, double r_param, size_t max_indel_size) :
		_basic_mu(basic_mu), _basic_gamma(basic_gamma), _r_param(r_param), _max_indel_size(max_indel_size)
	{
		if ((_r_param < 0) || (_r_param >= 1.0))
		{
			cerr << "r should be between 0 and 1 (not including 1)" << endl;
			exit(1);
		}
		if (_basic_mu < 0)
		{
			cerr << "mu should be non negative" << endl;
			exit(1);
		}
		_type_rate = 0; // long indel
		cout << "constructing the long indel model" << endl;
		cout << "basic_mu: " << _basic_mu << endl;
		cout << "basic_gamma: " << _basic_gamma << endl;
		cout << "r_param: " << _r_param << endl;
		cout << "max_indel_size: " << _max_indel_size << endl;

		set_rate_vectors();
		set_rate_cumu_vecs();
	} //constructor long indel

	rates(double IR, double A, size_t max_indel_size) :
		_IR(IR), _A(A), _max_indel_size(max_indel_size)
	{
		if ((_IR < 0.0) || ((_IR > 1.0)))
		{
			cerr << "IR should be between 0 and 1" << endl;
			exit(1);
		}
		if (_A <= 1.0)
		{
			cerr << "A (power law shape parameter) should be greater than 1" << endl;
			exit(1);
		}
		_type_rate = 1; // power law
		cout << "constructing the power law model" << endl;
		cout << "basic_A " << _A << endl;
		cout << "basic_IR: " << _IR << endl;
		cout << "max_indel_size: " << _max_indel_size << endl;

		set_rate_vectors();
		set_rate_cumu_vecs();
	} //constructor power law

	vector<double> get_mus() const { return _mus; }
	vector<double> get_lambdas() const { return _lambdas; }
	double get_total_rate_to_change() const { return _total_rate_to_change; }
	double get_basic_gamma() const { return _basic_gamma; }
	size_t get_max_indel_size() const { return _max_indel_size; }
	void sample_event_type(size_t & type, size_t & size) const;
	double sample_jump_time(size_t num_positions) const;

	double sample_insertion_jump_time() const; // can be used if sequence goes to 0 (immortal link)
	void sample_insertion_size(size_t & size) const; // can be used if sequence goes to 0 (immortal link)


private:
	size_t _type_rate; // if 0 - long indel (by formulae for N), if 1 - by INDELible power law

	double _basic_mu; // long indel
	double _basic_gamma; // long indel
	double _r_param; // long indel

	size_t _max_indel_size;

	double _IR; // per site insertion rate = per site deletion rate
	double _A; // power law shape parameter
	double _riemann_zeta_A = 0.0;
	
	vector<double> _mus;
	vector<double> _lambdas;
	vector<double> _cumu_mus;
	vector<double> _cumu_lambdas;
	double _total_rate_to_change; // sum of all indel rates per position
	double _total_insertion_rate; // sum of all insertion rates per position
	double _total_deletion_rate; // sum of all deletion rates per position

	void set_rate_vectors();
	void set_rate_cumu_vecs();

};

#endif
