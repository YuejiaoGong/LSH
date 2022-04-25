#ifndef BNDE_H
#define BNDE_H

#include <map>
using namespace std;

CEC2013 *pFunc;		

#define NR      50
#define MAXDIM  20
#define MAXFES  400000
#define POPSIZE 1000

/*parameters of the algorithm*/
int popsize;  
int incsize;
double F = 0.5;
double CR = 0.9;

struct Individual
{
	double x[MAXDIM];
	double u[MAXDIM];
	double fit;
	double ufit;
	double CR;
	double F;
};
Individual pop[POPSIZE];
Individual centroid;
int threshold = 30;

struct Projection
{
	double a1[MAXDIM], a2[MAXDIM];  // random Gaussian vector
	double b1, b2;          // offset
	double r1, r2;          // segment width
	double p1[POPSIZE], p2[POPSIZE]; // projection
	int hashid[POPSIZE];
	map<int, vector<int> > bucket;
};

Projection proj[30];
int num_projs;
map<int, vector<int> >::iterator iter;

double gbestval;            //use the found global best fitness to estimate the global optima

#define betterthan >=
#define worsethan  <
const double accuracy_level[] = {0.1, 0.01, 0.001, 0.0001, 0.00001};//five levels of accuracy 
int dim;
int maxfes,fes;
int gen;
double xmax[MAXDIM];
double xmin[MAXDIM];
vector<vector<double>> solu; 
vector<vector<double>> seed; //solu & seed are used in the final output of the algorithm

/*performance measure*/
double peak_ratio[5];
double succ_rate[5];
double run_time;

inline double rand_uni(double low,double high);
inline double rand_gau(double mu, double thegma);
inline void boundsctl(double & x_id, int d);
inline double eu_dis(double p1[], double p2[], int dd);
void outputsolu();
double caleval(double pos[]);
int cmp ( const void *a , const void *b );
double rand_cauchy(double alpha, double beta);

#endif