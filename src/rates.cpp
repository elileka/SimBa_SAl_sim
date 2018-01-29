// code by Eli Levy Karin

#include "rates.h"

void rates::sample_event_type(size_t & type, size_t & size) const
{
	type = 0; // denotes deletion

	random_device rd;  // will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
	uniform_real_distribution<double> distribution(0.0, _total_rate_to_change);
	double rand_in_range = distribution(gen);

	// we can search nicer (lion in the desert) but not now
	for (size_t i = 0; i < _cumu_mus.size(); i++)
	{
		size_t curr_size = (i + 1);
		if (rand_in_range <= _cumu_mus[i])
		{
			size = curr_size;
			return;
		}
	}

	// set indel type to insertion
	type = 1; // denotes insertion
	for (size_t i = 0; i < _cumu_lambdas.size(); i++)
	{
		size_t curr_size = (i + 1);
		if (rand_in_range <= _cumu_lambdas[i])
		{
			size = curr_size;
			return;
		}
	}
}


double rates::sample_jump_time(size_t num_positions) const
{
	double rate_to_change_over_all_positions = _total_rate_to_change * num_positions;
	// sample the time to the next jump
	random_device rd;  // Will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
	uniform_real_distribution<double> distribution(0.0, 1.0);
	double uniform_helper = distribution(gen);
	double time_to_jump = -(1.0 / rate_to_change_over_all_positions) * log(uniform_helper);
	return(time_to_jump);
}

double rates::sample_insertion_jump_time() const
{
	// sample the time to the next jump of an insertion from the immortal link: _total_insertion_rate
	random_device rd;  // Will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
	uniform_real_distribution<double> distribution(0.0, 1.0);
	double uniform_helper = distribution(gen);
	double time_to_jump = -(1.0 / _total_insertion_rate) * log(uniform_helper);
	return(time_to_jump);
}

void rates::sample_insertion_size(size_t & size) const
{
	random_device rd;  // will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
	uniform_real_distribution<double> distribution(_total_deletion_rate, _total_rate_to_change);
	double rand_in_range = distribution(gen);

	for (size_t i = 0; i < _cumu_lambdas.size(); i++)
	{
		size_t curr_size = (i + 1);
		if (rand_in_range <= _cumu_lambdas[i])
		{
			size = curr_size;
			return;
		}
	}
}

void rates::set_rate_cumu_vecs()
{
	// rate vectors are already initialized
	_cumu_mus.push_back(_mus[0]);

	for (size_t i = 1; i < _mus.size(); i++)
	{
		_cumu_mus.push_back(_cumu_mus.back() + _mus[i]);
	}

	double total_mus = _cumu_mus.back();

	_cumu_lambdas.push_back(total_mus + _lambdas[0]);
	for (size_t i = 1; i < _lambdas.size(); i++)
	{
		_cumu_lambdas.push_back(_cumu_lambdas.back() + _lambdas[i]);
	}

}

void rates::set_rate_vectors()
{
	// sums of all rates per position:
	_total_rate_to_change = 0.0;
	_total_deletion_rate = 0.0;
	_total_insertion_rate = 0.0;
	

	if (_type_rate == 0) // long indel model - N rate function
	{
		double epsilon = 0.0001;
		if (abs(_r_param - 0) < epsilon) // TKF91 model
		{
			double mu_1 = _basic_mu;
			double lambda_1 = mu_1 * _basic_gamma; // here we don't account for the type of inserted characters

			_total_rate_to_change += (mu_1 + lambda_1);
			_total_deletion_rate += mu_1;
			_total_insertion_rate += lambda_1;
			
			_mus.push_back(mu_1);
			_lambdas.push_back(lambda_1);
		}
		else
		{
			for (size_t k = 1; k <= _max_indel_size; k++)
			{
				double curr_mu = _basic_mu * pow((1.0 - _r_param), 2) * pow((_r_param), (k - 1));
				double curr_lambda = curr_mu * pow(_basic_gamma,k); // here we don't account for the type of inserted characters

				_total_rate_to_change += (curr_mu + curr_lambda);
				_total_deletion_rate += curr_mu;
				_total_insertion_rate += curr_lambda;
				
				_mus.push_back(curr_mu);
				_lambdas.push_back(curr_lambda);
			}

			// add the rate for longer indels to the truncated vector:
			double remainder_mu = _basic_mu * (1.0 - _r_param) * pow((_r_param), _max_indel_size);
			double remainder_lambda = _basic_gamma * _basic_mu * pow((1.0 - _r_param), 2) * pow((_r_param * _basic_gamma), _max_indel_size) / (1.0 - (_r_param * _basic_gamma));

			_mus[(_max_indel_size - 1)] += remainder_mu;
			_lambdas[(_max_indel_size - 1)] += remainder_lambda;

			_total_rate_to_change += (remainder_mu + remainder_lambda);
			_total_deletion_rate += remainder_mu;
			_total_insertion_rate += remainder_lambda;
		}
	}

	if (_type_rate == 1) // power law inspired model
	{
		//_total_rate_to_change = 2 * _IR;
		//_total_deletion_rate = _IR;
		//_total_insertion_rate = _IR;
		
		// we generate the spread of the rates according to a truncated power law distribution
		vector<double> event_pow_probs;
		double sum_unnormalized_probs = 0.0;
		for (size_t k = 1; k <= _max_indel_size; k++)
		{
			double curr_unnormalized_prob = pow(k,(- _A));
			sum_unnormalized_probs = sum_unnormalized_probs + curr_unnormalized_prob;
			event_pow_probs.push_back(curr_unnormalized_prob); // before noramlization
		}
		for (size_t i = 0; i < event_pow_probs.size(); i++)
		{
			event_pow_probs[i] = event_pow_probs[i] / sum_unnormalized_probs; // normalize
		}

		// now prepare the rate vectors
		for (size_t k = 1; k <= _max_indel_size; k++)
		{
			double curr_mu = event_pow_probs[k - 1] * _IR / k;
			double curr_lambda = curr_mu; // equal deletion and insertion

			_total_rate_to_change += (curr_mu + curr_lambda);
			_total_deletion_rate += curr_mu;
			_total_insertion_rate += curr_lambda;

			_mus.push_back(curr_mu);
			_lambdas.push_back(curr_lambda);
		}
	}

}
