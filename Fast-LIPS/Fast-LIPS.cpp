#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <map>
#include "Fast-LIPS.h"
#include "cec2013.h"
using namespace std;

#define TIMES	   50
#define betterthan >=

struct particle_type gBest;				// global best
struct particle_type lBest;				// local best
struct particle_type particle[MAXSIZE]; // partical swarm
struct particle_type pBest[MAXSIZE];	// personal best
struct particle_type centroid;

int swarm_size;				   // swarm size
int dim;			           // dim of the problem

double chi;					   // constriction coefficient
double phi1, phi2;			   // acceeleration coefficient
double vmax[MAXDIM];		   // max velocity of each particle
double lbound[MAXDIM];	       // lower bound of each dim
double ubound[MAXDIM];	       // upper bound of each dim
int    max_eval;			   // max number of evaluations
int    num_eval;			   // current number of evaluations
int    max_iter;			   // max number of iterations
int    num_iter;			   // current number of iterations
int    neighbor[MAXSIZE][MAXSIZE];
int    nsize = 5;			   // maximum number of particles _informed_ by a given one, i.e. number of neighbors


CEC2013 *pFunc;

double calc_fitness(double* x)
{
	return pFunc->evaluate(x);
}

struct Projection
{
	double a1[MAXDIM], a2[MAXDIM];  // random Gaussian vector
	double b1, b2;          // offset
	double r1, r2;          // segment width
	double p1[MAXSIZE], p2[MAXSIZE]; // projection
	int hashid[MAXSIZE];
	map<int, vector<int> > bucket;
};
Projection proj[200];
int num_projs;     // number of projection lines
map<int, vector<int> >::iterator iter;

void generate_projs()
{
	int i, j, k, r1, r2, cnt, id;
	double min1, max1, width;
	double min2, max2;

	for(i=0; i<num_projs; i++){
		
		for(j=0; j<dim; j++){
			proj[i].a1[j] = rand_gau(0, 1);
			proj[i].a2[j] = rand_gau(0, 1);
		}

		for(iter=proj[i].bucket.begin(); iter!=proj[i].bucket.end(); iter++)
			iter->second.clear();
		proj[i].bucket.clear();
		min1 = min2 = 1e300;
		max1 = max2 = -1e300;
		for(j=0; j<swarm_size; j++){
			proj[i].p1[j] = 0;
			proj[i].p2[j] = 0;
			for(k=0; k<dim; k++){
				proj[i].p1[j] += proj[i].a1[k]*(pBest[j].x[k]-centroid.x[k]);
				proj[i].p2[j] += proj[i].a2[k]*(pBest[j].x[k]-centroid.x[k]);
			}
			if(proj[i].p1[j]>max1) max1 = proj[i].p1[j];
			if(proj[i].p1[j]<min1) min1 = proj[i].p1[j];

			if(proj[i].p2[j]>max2) max2 = proj[i].p2[j];
			if(proj[i].p2[j]<min2) min2 = proj[i].p2[j];
		}

		double d = 20;  // number of buckets
		proj[i].r1 = (max1-min1)/(d);
		proj[i].r2 = (max2-min2)/(d);
		proj[i].b1 = rand_uni(0, proj[i].r1);
		proj[i].b2 = rand_uni(0, proj[i].r2);
		for(j=0; j<swarm_size; j++){
			id = ceil((proj[i].p1[j]+proj[i].b1)/proj[i].r1) * 50 + ceil((proj[i].p2[j]+proj[i].b2)/proj[i].r2);
			proj[i].hashid[j]  = id;
			proj[i].bucket[id].push_back(j);
		}
	}
}

void define_neighbors()
{
	int i, j, r1, r2, cnt;
	int begin, end;
	int np, id, bs;

	generate_projs();

	for(i=0; i<swarm_size; i++){
		for(j=0; j<nsize; j++){
			cnt = 0;
			//do
			{
			np  =  rand()%num_projs;
			id = proj[np].hashid[i];
			bs = proj[np].bucket[id].size();
			cnt++;
			}
			//while(bs>0.05*swarm_size && bs>1 && cnt<100);
			neighbor[i][j] = proj[np].bucket[id][rand()%bs];
		}
	}
}

void initialize()
{
	/* initialize parameters */
	chi = 0.729843788;
	phi1 = phi2 = 2.05;

	num_eval = 0;
	num_iter = 0;

	for(int j=0; j<dim; j++)
	{
		vmax[j] = 0.5*(ubound[j]-lbound[j]);
	}

	/* initialize velocity and position */
	for(int i=0; i<swarm_size; i++)
	{
		for(int j=0; j<dim; j++)
		{
			particle[i].v[j] = rand_uni(-vmax[j], vmax[j]);
			particle[i].x[j] = rand_uni(lbound[j], ubound[j]);
			centroid.x[j] = (lbound[j]+ubound[j])/2;
		}
		particle[i].fitness = calc_fitness(particle[i].x);
		num_eval++;
	}
	

	/* initialize pBest and gBest */
	int best = 0;
	for(int i=0; i<swarm_size; i++)
	{
		pBest[i] = particle[i];
		if(pBest[i].fitness betterthan pBest[best].fitness)
			best = i;
	}
	gBest = pBest[best];
}


void swap(int i, int j)
{
	struct particle_type t1 = particle[i];
	struct particle_type t2 = pBest[i];
	particle[i] = particle[j];
	pBest[i] = pBest[j];
	particle[j] = t1;
	pBest[j] = t2;
}

void process()
{
	int multiple;

	int nsize = 2;
	if(num_iter<=0.25*max_iter)
		nsize=2;
    else if(num_iter<=0.5*max_iter)
        nsize=3;
    else if(num_iter<=0.75*max_iter)
        nsize=4;
    else
        nsize=5;

	num_projs = 10; // number of projections (hash functions)

	if(num_iter%(2*num_projs)==0)
		define_neighbors();

	
	for(int i=0; i<swarm_size; i++)
	{
		// find the lBest

		double phi[MAXDIM] = {0.0};
		struct particle_type P;
		int indices[5];

		for(int j=0; j<dim; j++)
			P.x[j] = 0;
		for(int nn=0; nn<nsize; nn++)
		{
			double r;
			for(int j=0; j<dim; j++)
			{
				r = rand_uni(0, 4.1/nsize);
				P.x[j] += r*pBest[neighbor[i][nn]].x[j];
				phi[j] += r;
			}
		}
		for(int j=0; j<dim; j++)
		{
			P.x[j] = P.x[j]/phi[j];
		}


		// update velocity and position 
		for(int j=0; j<dim; j++)
		{
			particle[i].v[j] = chi*particle[i].v[j]+phi[j]*(P.x[j]-particle[i].x[j]);


			if(fabs(particle[i].v[j])>vmax[j])
			{
				particle[i].v[j] = rand_uni(-vmax[j]/2, vmax[j]/2);
			}

			particle[i].x[j] += particle[i].v[j]; 

			/* Attention: bound handling technique has a great impact on the performance of f1, because its optimum happen to be on the boundary */
			/* adjust if the position has been out of range */
			if(particle[i].x[j]<lbound[j])
			{
				particle[i].x[j] = lbound[j];
				particle[i].v[j] = 0;
			}
			if(particle[i].x[j]>ubound[j])
			{
				particle[i].x[j] = ubound[j];
				particle[i].v[j] = 0;
			}
		}

		// calculate fitness
		particle[i].fitness = calc_fitness(particle[i].x);
		num_eval++;
	}

	/* use synchronous update according to the author's source code */
	for(int i=0; i<swarm_size; i++)
	{
		// update pBest and gBest
		if(particle[i].fitness betterthan pBest[i].fitness)
		{
			pBest[i] = particle[i];
		}
		if(pBest[i].fitness betterthan gBest.fitness)
		{
			gBest  = pBest[i];
		}
	}
	num_iter++;
}

int main()
{
	srand(time(0));
	FILE* fp;
	int num_opt;     // number of global optima
	int num_peak[5]; // number of peaks that have been found
	double peak_ratio[5];
	double succ_rate[5];  
	double accuracy_level[5] = {0.1, 0.01, 0.001, 0.0001, 0.00001};

	FILE* total_pr = fopen("total_pr.txt", "a");
	FILE* total_sr = fopen("total_sr.txt", "a");
    fprintf(total_pr,"func\taccruacy level\t1.0E-01\t1.0E-02\t1.0E-03\t1.0E-04\t1.0E-05\n");
	fclose(total_pr);
	fprintf(total_sr,"func\taccruacy level\t1.0E-01\t1.0E-02\t1.0E-03\t1.0E-04\t1.0E-05\n");
	fclose(total_sr);

	for(int fn=1; fn<=20; fn++)
	{
		pFunc = new CEC2013(fn);
		dim   = pFunc->get_dimension();
		swarm_size  = 100; 
		max_eval	= pFunc->get_maxfes();
		max_iter    = max_eval/swarm_size;
		for(int d=0; d<dim; d++)
		{
			lbound[d] = pFunc->get_lbound(d);
			ubound[d] = pFunc->get_ubound(d);
		}
		num_opt = pFunc->get_no_goptima();

		printf("F%d is running...\n", fn);
		printf("%d\n", num_opt);
				
		total_pr = fopen("total_pr.txt", "a");
		fprintf(total_pr, "%d\t", fn);
		total_sr = fopen("total_sr.txt", "a");
		fprintf(total_sr, "%d\t", fn);

		FILE* stat_pr = fopen("stat_pr.txt", "a");
		FILE* stat_np = fopen("stat_np.txt", "a");

		memset(peak_ratio, 0, sizeof(peak_ratio));
		memset(succ_rate, 0, sizeof(succ_rate));

		for(int i=0; i<TIMES; i++)
		{
			initialize();
			while(num_eval<max_eval)
			{
				process();
			}

			vector<vector<double> > res, seed;
			convert_to_vec(res);
			for(int j=0; j<5; j++)	//five levels of accuracy
			{ 
				num_peak[j] = how_many_goptima(res,seed,pFunc,accuracy_level[j],pFunc->get_rho());
				printf("%d\t", num_peak[j]);
				fprintf(stat_np, "%d\t", num_peak[j]);
				fprintf(stat_pr, "%lf\t", (double)num_peak[j]/num_opt);
				peak_ratio[j] += num_peak[j];
				if(num_peak[j] >= num_opt)
					succ_rate[j]++;
			}
			printf("\n");
			fprintf(stat_np, "\n");
			fprintf(stat_pr, "\n");
		}
		fprintf(stat_np, "\n");
		fprintf(stat_pr, "\n");

		fprintf(total_pr, "peak_ratio\t");
		for(int j=0; j<5; j++){
			peak_ratio[j] /= (TIMES*num_opt + 0.0);
			succ_rate[j] /= (TIMES + 0.0);
			printf("level %d : %f\t%f\n", j+1, peak_ratio[j],succ_rate[j]);
			fprintf(total_pr, "%f\t", peak_ratio[j]);
		}
		fprintf(total_pr, "\n");
		fprintf(total_sr, "\tsuccess_rate\t");
		for(int j=0; j<5; j++)fprintf(total_sr, "%f\t", succ_rate[j]);
		fprintf(total_sr, "\n");
		
		fclose(total_pr);
		fclose(total_sr);
		fclose(stat_pr);
		fclose(stat_np);
	}
	getchar();
	return 0;
}

/* auxiliary functions */
inline double rand_uni(double low, double high)//generate uniformly random numbers
{
	return (double(rand())/RAND_MAX)*(high-low)+low;
}

void convert_to_vec(vector<vector<double> >& s)
{
	vector<double> e;
	s.clear();
	for(int i=0; i<swarm_size; i++)
	{
		e.clear();
		for(int j=0; j<dim; j++)
		{
			e.push_back(pBest[i].x[j]);
		}
		s.push_back(e);
	}
}

static int phase = 0;
double rand_gau(double mu, double thegma)//generate Gaussian distributed random numbers
{
	static double V1, V2, S; 
	double X;    
	if ( phase == 0 ) {
		do {
			double U1 = (double)rand() /(double)RAND_MAX;
			double U2 = (double)rand() /(double)RAND_MAX;
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;
	return mu+X*thegma;
}