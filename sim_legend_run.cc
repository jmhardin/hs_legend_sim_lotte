#include <random>
#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cmath>

#define STAR_FLOOR 26
#define MAX_GAMES 1000

double calc_binom(int in, int im)
{
	double rval = 1;
	int n = in;
	int m = im;
	if (m > n/2)
	{
		m = n - m;
	}

	for (int i = 0; i < m; i++)
	{
		rval *= (n - i);
		rval /= (m - i);
	}
	return rval;
}
double fastpow(double b, int e)
{//Positive integers only
	double rval = 1;
	double cpow = b;
	int expon = e;	
	
	while (expon > 0)
	{
		if (expon % 2 != 0)
		{
			rval *= cpow;
		}
		cpow *= cpow;
		expon /= 2;
	}
	if (e > 0)
	{
		return rval;
	}
	else
	{
		return 1/rval;
	}

}
double calc_pdf(
	double win_rate,\
	int stars,\
	int games)
{
	if (games < stars || (games - stars) %2 != 0)
	{//not possible, therefore doesn't make sense
		return 0;
	}
	double rval = fastpow(win_rate*(1-win_rate),(games-stars)/2)*fastpow(win_rate,stars)*calc_binom(games,(games-stars)/2)*stars/games;
	rval = std::max(rval,0.0);
	return rval;
}
double calc_pdf_bayes(
	double win_rate,\
	int stars,\
	int games,
	double win_rate_error)
{
	int divisions = 1000;
	double half_width = 3*win_rate_error;
	double dwr = 2*half_width/divisions;

	double rval = 0;
	for (double wr = win_rate - half_width; wr < win_rate + half_width; wr += dwr)
	{
		rval += dwr*exp(-(win_rate - wr)*(win_rate - wr)/(2*win_rate_error*win_rate_error))*calc_pdf(wr,stars,games);
	}
	rval /= (sqrt(2*3.14159)*win_rate_error);
	return rval;
}


unsigned int sim_run(
	double win_rate,\
	int stars_needed,\
	std::mt19937_64 &rgen)
{
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	unsigned int rval = 0;
	int sn = stars_needed;
	while (rval < MAX_GAMES)
	{
		if (win_rate > distribution(rgen))
		{//won a game
			sn--;
			if (sn <= 0)
			{//made legend
				break;
			}
		}
		else
		{//lost a game
			if (sn < STAR_FLOOR)
			{//above rank 5, 0 stars
				sn++;
			}
		}
		rval++;
	}
	return rval;
}

void fill_calc_pdf(
	double win_rate,\
	int stars,\
	std::vector<double> &pdf_probs)
{
	pdf_probs.clear();
	for (unsigned int i = 0; i < MAX_GAMES; i++)
	{
		pdf_probs.push_back(0);
		pdf_probs[i] = calc_pdf(win_rate,stars,i);
	}
}
void fill_calc_pdf_bayes(
	double win_rate,\
	int stars,\
	double win_rate_error,\
	std::vector<double> &pdf_probs)
{
	pdf_probs.clear();
	for (unsigned int i = 0; i < MAX_GAMES; i++)
	{
		pdf_probs.push_back(0);
		pdf_probs[i] = calc_pdf_bayes(win_rate,stars,i,win_rate_error);
	}
}
void fill_sim_pdf(
	double win_rate,\
	int stars,\
	int num_trials,\
	std::vector<double> &pdf_probs,\
	std::mt19937_64 &rgen)
{
	pdf_probs.clear();
	for (unsigned int i = 0; i < MAX_GAMES; i++)
	{
		pdf_probs.push_back(0);
	}
	for (unsigned int i = 0; i < num_trials; i++)
	{
		unsigned int g = sim_run(win_rate,stars,rgen);
		pdf_probs[g] += 1.0;
	}
	for (unsigned int i = 0; i < MAX_GAMES; i++)
	{
		pdf_probs[i] /= num_trials;
	}
}
void fill_sim_pdf_bayes(
	double win_rate,\
	int stars,\
	int num_trials,\
	double win_rate_error,\
	std::vector<double> &pdf_probs,\
	std::mt19937_64 &rgen)
{
	std::normal_distribution<double> wr_dist(win_rate,win_rate_error);
	pdf_probs.clear();
	for (unsigned int i = 0; i < MAX_GAMES+1; i++)
	{
		pdf_probs.push_back(0);
	}
	for (unsigned int i = 0; i < num_trials; i++)
	{
		double tmp_win_rate = wr_dist(rgen);
		unsigned int g = sim_run(tmp_win_rate,stars,rgen);
		pdf_probs[g] += 1.0;
	}
	for (unsigned int i = 0; i < MAX_GAMES; i++)
	{
		pdf_probs[i] /= num_trials;
	}
}
	

int main( int argc, char** argv)
{
	std::random_device dev;
	std::mt19937_64 rgen(dev());
	bool print_cdf = true;
	double wr = .53;
	double wr_err = .01;
	int stars = 20;
	int trials = 100000;

	bool sucessful_options = true;
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i],"-n") == 0)
		{
			i++;
			trials = atoi(argv[i]);
		}
		else if (strcmp(argv[i],"-wr") == 0)
		{
			i++;
			wr = atof(argv[i]);
		}
		else if (strcmp(argv[i],"-wr_err") == 0)
		{
			i++;
			wr_err = atof(argv[i]);
		}
		else if (strcmp(argv[i],"-stars") == 0)
		{
			i++;
			stars = atoi(argv[i]);
		}
		else if (strcmp(argv[i],"-pdf") == 0)
		{
			print_cdf = false;
		}
		else
		{
			sucessful_options = false;
			break;
		}
	}
	if (sucessful_options == false)
	{
		printf("USAGE:\n");
		printf("-n <trials>\n");
		printf("-wr <win rate>\n");
		printf("-wr_err <gaussian error on win rate>\n");
		printf("-stars <needed stars>\n");
		printf("-pdf (Prints the pdf instead of the cdf)\n");
		return -1;
	}

	std::vector<double> sim_pdf;
	std::vector<double> sim_pdf_bayes;
	std::vector<double> calc_pdf;
	std::vector<double> calc_pdf_bayes;

	std::vector<double> sim_cdf;
	std::vector<double> sim_cdf_bayes;
	std::vector<double> calc_cdf;
	std::vector<double> calc_cdf_bayes;

	fill_sim_pdf(\
		wr,\
		stars,\
		trials,\
		sim_pdf,\
		rgen);
	fill_sim_pdf_bayes(\
		wr,\
		stars,\
		trials,\
		wr_err,\
		sim_pdf_bayes,\
		rgen);
	
	fill_calc_pdf(
		wr,\
		stars,\
		calc_pdf);
	fill_calc_pdf_bayes(
		wr,\
		stars,\
		wr_err,\
		calc_pdf_bayes);
	for (unsigned int i = 0; i < MAX_GAMES; i++)
	{	
		sim_cdf.push_back(0);
		sim_cdf_bayes.push_back(0);
		calc_cdf.push_back(0);
		calc_cdf_bayes.push_back(0);
	}

	
	sim_cdf[0] = sim_pdf[0];
	sim_cdf_bayes[0] = sim_pdf_bayes[0];
	calc_cdf[0] = calc_pdf[0];
	calc_cdf_bayes[0] = calc_pdf_bayes[0];

	for (unsigned int i = 1; i < sim_cdf.size(); i++)
	{
		sim_cdf[i] = sim_cdf[i-1]+sim_pdf[i];
		sim_cdf_bayes[i] = sim_cdf_bayes[i-1]+sim_pdf_bayes[i];
		calc_cdf[i] = calc_cdf[i-1]+calc_pdf[i];
		calc_cdf_bayes[i] = calc_cdf_bayes[i-1]+calc_pdf_bayes[i];
	}

	for (unsigned int i = 0; i < MAX_GAMES; i++)
	{
		if (print_cdf == true)
		{
			printf("%8d %12.04e %12.04e %12.04e %12.04e\n",i,sim_cdf[i],sim_cdf_bayes[i],calc_cdf[i],calc_cdf_bayes[i]);
		}
		else
		{
			printf("%8d %12.04e %12.04e %12.04e %12.04e\n",i,sim_pdf[i],sim_pdf_bayes[i],calc_pdf[i],calc_pdf_bayes[i]);
		}
	}
	return 0;
}
