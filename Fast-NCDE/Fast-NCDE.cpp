#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cstdio>
#include<vector>
#include<ctime>
using std::vector;

#include"cec2013.h"
#include"Fast-NCDE.h"

void initialization()
{
	int i, j;
	fes = 0;
	gbestval = -1e300;
	popsize = 100;
	incsize = 0;

	for(i=0; i<popsize; i++){
		for(j=0; j<dim; j++){
			pop[i].x[j] = rand_uni(xmin[j], xmax[j]);
			centroid.x[j] = (xmin[j]+xmax[j])/2.0;
		}
		pop[i].fit = caleval(pop[i].x);
		if(pop[i].fit betterthan gbestval)
			gbestval = pop[i].fit;

	}
	gen = 0;
	num_projs = 10;   // number of projection lines (hash functions)
}


void crossover_mutation(int i)
{
	int r1, r2, r3, jrand;

	int np, id, bs;
	int cnt = 0;
	do
	{
		np  =  rand()%num_projs;
		id = proj[np].hashid[i];
		bs = proj[np].bucket[id].size();
		cnt++;
		if(cnt>num_projs)
			break;
	}
	while(bs<3);

	if(cnt>num_projs){
		do{r1=rand()%popsize;}while(r1==i);
		do{r2=rand()%popsize;}while(r2==r1 || r2==i);
		do{r3=rand()%popsize;}while(r3==r2 || r3==r1 || r3==i);
	}
	else{
		r1 = proj[np].bucket[id][rand()%bs];

		do
		{
			np  =  rand()%num_projs;
			id = proj[np].hashid[i];
			bs = proj[np].bucket[id].size();
		}
		while(bs<3);
		r2 = proj[np].bucket[id][rand()%bs];

		do
		{
			np  =  rand()%num_projs;
			id = proj[np].hashid[i];
			bs = proj[np].bucket[id].size();
		}
		while(bs<3);
		r3 = proj[np].bucket[id][rand()%bs];
	}

	jrand=rand()%dim;
	for(int j=0; j<dim; j++){
		if(rand_uni(0,1)<=CR || j==jrand){
			pop[i].u[j]=pop[r1].x[j]+F*(pop[r2].x[j]-pop[r3].x[j]);
			boundsctl(pop[i].u[j],j);
		}
		else 
			pop[i].u[j]=pop[i].x[j];
	}	
}

int find_nearest(int i)
{
	int nei;
	double mindis = 1e300, dis;
	for(int j=0; j<popsize; j++){
		dis = eu_dis(pop[i].u, pop[j].x, dim);
		if(dis < mindis){
			mindis = dis;
			nei = j;
		}
	}
	return nei;
}

int find_nearest_in_bucket(int i)
{
	int nei = -1;
	double mindis = 1e300, dis;
	int np = rand()%num_projs;
	double p1=0, p2=0;
	int k, id;
	for(k=0; k<dim; k++){
		p1 += proj[np].a1[k]*(pop[i].u[k]-centroid.x[k]);
		p2 += proj[np].a2[k]*(pop[i].u[k]-centroid.x[k]);
	}
	id = ceil((p1+proj[np].b1)/proj[np].r1) * 50 + ceil((p2+proj[np].b2)/proj[np].r2);
	vector<int> bucket = proj[np].bucket[id];
	int bs = bucket.size();
	for(int j=0; j<bs; j++){
		dis = eu_dis(pop[i].u, pop[bucket[j]].x, dim);
		if(dis < mindis){
			mindis = dis;
			nei = bucket[j];
		}
	}
	if(nei == -1){
		nei = 0;
	}
	return nei;
}

void selection()
{
	int i, j, k, n, id;
	double min, max;

	// generate random projections
	if(gen%(2*num_projs)==0)
	{
		double min1, max1;
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
			for(j=0; j<popsize; j++){
				proj[i].p1[j] = 0;
				proj[i].p2[j] = 0;
				for(k=0; k<dim; k++){
					proj[i].p1[j] += proj[i].a1[k]*(pop[j].x[k]-centroid.x[k]);
					proj[i].p2[j] += proj[i].a2[k]*(pop[j].x[k]-centroid.x[k]);
				}
				if(proj[i].p1[j]>max1) max1 = proj[i].p1[j];
				if(proj[i].p1[j]<min1) min1 = proj[i].p1[j];

				if(proj[i].p2[j]>max2) max2 = proj[i].p2[j];
				if(proj[i].p2[j]<min2) min2 = proj[i].p2[j];
			}

			double d = 5; // number of buckets
			//proj[i].r = (max-min)/(20+rand()%15);
			proj[i].r1 = (max1-min1)/(d);
			proj[i].r2 = (max2-min2)/(d);
			proj[i].b1 = rand_uni(0, proj[i].r1);
			proj[i].b2 = rand_uni(0, proj[i].r2);
			for(j=0; j<popsize; j++){
				id = ceil((proj[i].p1[j]+proj[i].b1)/proj[i].r1) * 50 + ceil((proj[i].p2[j]+proj[i].b2)/proj[i].r2);
				proj[i].hashid[j]  = id;
				proj[i].bucket[id].push_back(j);
			}
		}

	}

	for(i=0; i<popsize; i++){
		crossover_mutation(i);
	//}

	//for(i=0; i<popsize; i++){
		
		pop[i].ufit = caleval(pop[i].u);

		n = find_nearest_in_bucket(i); 
		if(pop[i].ufit betterthan pop[n].fit){
			for(j=0; j<dim; j++)
				pop[n].x[j] = pop[i].u[j];
			pop[n].fit = pop[i].ufit;
		}

	}

	gen++;
}

void process()
{
	initialization();
	while(fes < maxfes)
	{	
		selection();	
	}
	outputsolu();
}

void main()
{	
	int func_id, d;
	int nkp;    //the number of known global optima
	int nfp[5]; //the number of found global optima (with different accuracy levels)
	double start_time, end_time;

	srand((unsigned)time(NULL));
	FILE* total_pr = fopen("total_pr.txt", "a");
	FILE* total_sr = fopen("total_sr.txt", "a");
    fprintf(total_pr,"func\taccruacy level\t1.0E-01\t1.0E-02\t1.0E-03\t1.0E-04\t1.0E-05\n");
	fclose(total_pr);
	fprintf(total_sr,"func\taccruacy level\t1.0E-01\t1.0E-02\t1.0E-03\t1.0E-04\t1.0E-05\n");
	fclose(total_sr);
	
	for(func_id = 1; func_id <= 20; func_id ++){
		/*initialize benchmark instance*/
		pFunc = new CEC2013(func_id);
		dim = pFunc->get_dimension();	
		maxfes = pFunc->get_maxfes();
		for(d=0;d<dim;++d) {
			xmax[d]=pFunc->get_ubound(d);
			xmin[d]=pFunc->get_lbound(d);
		}
		printf("F%d is running...\n",func_id);
		total_pr = fopen("total_pr.txt", "a");
		fprintf(total_pr, "%d\t", func_id);
		total_sr = fopen("total_sr.txt", "a");
		fprintf(total_sr, "%d\t", func_id);

		FILE* stat_pr = fopen("stat_pr.txt", "a");
		FILE* stat_np = fopen("stat_np.txt", "a");

		/*run algorithm & output result*/
		nkp = pFunc->get_no_goptima();
		printf("%d\n", nkp);//////////
		memset(peak_ratio, 0, sizeof(peak_ratio));
		memset(succ_rate, 0, sizeof(succ_rate));
		run_time = 0;
		for(int i=0; i<NR; i++){
			start_time = clock();
			process();
			end_time = clock();
			run_time = (end_time-start_time)/double(CLOCKS_PER_SEC);
			printf("run time: %lf\n", run_time);

			for(int j=0; j<5; j++){ //five levels of accuracy
				nfp[j] = how_many_goptima(solu,seed,pFunc,accuracy_level[j],pFunc->get_rho());
				printf("%d\t", nfp[j]);///////////
				fprintf(stat_np, "%d\t", nfp[j]);
				fprintf(stat_pr, "%lf\t", (double)nfp[j]/nkp);
				peak_ratio[j] += nfp[j];
				if(nfp[j] >= nkp)succ_rate[j]++;
			}
			printf("\n");/////////
			fprintf(stat_np, "\n");
			fprintf(stat_pr, "\n");
		}
		fprintf(stat_np, "\n");
		fprintf(stat_pr, "\n");

		fprintf(total_pr, "peak_ratio\t");
		for(int j=0; j<5; j++){
			peak_ratio[j] /= (NR*nkp + 0.0);
			succ_rate[j] /= (NR + 0.0);
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
	system("Pause");
}

/*           auxiliary functions          */
inline double rand_uni(double low,double high)//generate uniformly random numbers
{
	return (double(rand())/RAND_MAX)*(high-low)+low;
}

static int phase = 0;
inline double rand_gau(double mu, double thegma)//generate Gaussian distributed random numbers
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

inline void boundsctl(double & x_id, int d)//bound control, to restrict the variable in the range
{
	if(x_id < xmin[d])
		x_id = xmin[d];
	else if(x_id > xmax[d])
		x_id = xmax[d];
}

inline double eu_dis(double p1[], double p2[], int dd)
{
	double dis = 0;
	for(int d=0; d<dd; d++)
		dis += (p1[d]-p2[d])*(p1[d]-p2[d]);
	return sqrt(dis);
}

void outputsolu()
{
	vector<double> s;
	solu.clear();
	for(int i=0; i<popsize; i++){
		s.clear();
		for(int j=0; j< dim; j++) s.push_back(pop[i].x[j]);
		solu.push_back(s);
	}
}

double caleval(double pos[])
{
	fes++;
	double fit = pFunc->evaluate(pos);
	if(fit betterthan gbestval) gbestval = fit;
	return fit;
}

int cmp ( const void *a , const void *b ) 
{ 
	struct Individual * p = (struct Individual *) a;
	struct Individual * q = (struct Individual *) b;
	return (p->fit - q->fit < 0) ? 1 : -1; ; 
} 

void reinitialize_ind(int i)
{
	for(int j=0; j<dim; j++)
	{
		pop[i].x[j] = rand_uni(xmin[j], xmax[j]);
	}
	pop[i].fit = caleval(pop[i].x);
}

double rand_cauchy(double alpha, double beta)
{
	return alpha+beta*tan(3.1415926*(double(rand())/RAND_MAX-0.5));
}