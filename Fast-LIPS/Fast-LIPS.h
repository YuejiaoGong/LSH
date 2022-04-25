#ifndef _RPSO_H
#define _RPSO_H

#include <vector>
using std::vector;

#define MAXDIM   20	// max dimension
#define MAXSIZE  1000    // max swarm size

struct particle_type
{
	double v[MAXDIM];	// velocity
	double x[MAXDIM];	// current position
	double fitness;	
};

void initialize();
void process();

/* auxiliary functions */
int random_int(int low, int high);
double  rand_uni(double low, double high);
void convert_to_vec(vector<vector<double> >& s);
double rand_gau(double mu, double thegma);

#endif