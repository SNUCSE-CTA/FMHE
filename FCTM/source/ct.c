/*
 * SMART: string matching algorithms research tool.
 * Copyright (C) 2012  Simone Faro and Thierry Lecroq
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 * 
 * contact the authors at: faro@dmi.unict.it, thierry.lecroq@univ-rouen.fr
 * download the tool at: http://www.dmi.unict.it/~faro/smart/
 *
 * This is an implementation of the Knuth Morris Pratt algorithm
 * in D. E. Knuth and J. H. Morris and V. R. Pratt. 
 * Fast pattern matching in strings. SIAM J. Comput., vol.6, n.1, pp.323--350, (1977).
 */

#include "include/define.h"
#include "include/timer.h"
//#include <sys/types.h>
//#include <sys/ipc.h>
//#include <sys/shm.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <random>
#include <functional>

#define BEGIN_FILTERING		{timer_start(fil_timer);}
#define BEGIN_PREPROCESSING	{timer_start(pre_timer);}
#define BEGIN_SEARCHING		{timer_start(sea_timer);}
#define BEGIN_VERIFY		{timer_start(ver_timer);}
#define END_FILTERING		{timer_stop(fil_timer);}
#define END_PREPROCESSING	{timer_stop(pre_timer);}
#define END_SEARCHING		{timer_stop(sea_timer);}
#define END_VERIFY			{timer_stop(ver_timer);}

/* global variables used for computing preprocessing and searching times */
clock_t start, end;
TIMER *fil_timer, *pre_timer, *sea_timer, *ver_timer;


void getPD(dtype *s, int n, int *PD);
void getFailure(dtype *x, int m, int *pd, int *failure);
int searchct(dtype *x, int m, dtype *y, int n) ;

int *candidates;

int main(int argc, char *argv[])
{
	std::default_random_engine generator;
	dtype *rawp, *rawt;
	int m, n;
	m=std::stoi(argv[1]);
	n=TSIZE;
	std::uniform_int_distribution<int> distribution(0, n-m);
	auto myrand = std::bind(distribution, generator);
	if(argc<2){
		printf("no arg\n");
		return 0;
	}
	std::string filename = TFILE;

	fil_timer=(TIMER *) malloc(sizeof(TIMER));
	sea_timer=(TIMER *) malloc(sizeof(TIMER));
	pre_timer=(TIMER *) malloc(sizeof(TIMER));
	ver_timer=(TIMER *) malloc(sizeof(TIMER));
	FILE *f;
	if((f = fopen(filename.c_str(), "r"))==NULL){
		printf("Cannot open the file\n");
		return 0;
	}
	int i;
//	rawp = (dtype *) malloc(sizeof(dtype) * m);
	rawt = (dtype *) malloc(sizeof(dtype) * n);
	if(isINT){
//		for(i=0;i<m;i++){
//			fscanf(f, "%d", &rawp[i]);
//		}
		if(isLARGE){
			for(i=0;i<n;i++){
				fscanf(f, "%d", &rawt[i]);
			}
		}
		else{
			for(i=0;i<n;i++){
				fscanf(f, "%hhd", &rawt[i]);
			}
		}
	}
	else{
		if(isLARGE){
//			for(i=0;i<m;i++){
//				fscanf(f, "%f", &rawp[i]);
//			}
			for(i=0;i<n;i++){
				fscanf(f, "%lf", &rawt[i]);
			}
		}
		else{
//			for(i=0;i<m;i++){
//				fscanf(f, "%lf", &rawp[i]);
//			}
			for(i=0;i<n;i++){
				fscanf(f, "%f", &rawt[i]);
			}
		}
	}
	int total_occ=0;
	double totaltime =0;
//	candidates=(int *) malloc(n*sizeof(int));
	int occ;
	for(int num=0;num<RNUM;num++){
		rawp=rawt+myrand();
		BEGIN_SEARCHING
		occ = searchct(rawp, m, rawt, n);
		END_SEARCHING
		totaltime += timer_elapsed(sea_timer);
		total_occ+=occ;
//		printf("found %d occurrences\n", occ);
//		for(int i=0;i<occ;i++){
//			for(int j=0;j<m;j++){
//				printf("%d ", rawt[candidates[i]+j]);
//			}
//			printf("\n");
//		}
//		printf("\n");
	}
	fprintf(stderr, "%d ", total_occ);
	printf("%.2lf", totaltime);

//	free(candidates);
	free(rawt);
	free(fil_timer);
	free(sea_timer);
	free(pre_timer);
	free(ver_timer);
	return 0;
}

void getPD(dtype *s, int n, int *PD){
	int st=-1;
	int *sti = (int *) malloc(n*sizeof(int));
	for(int i=0;i<n;i++){
		int index;
		while(st!=-1){
			index=sti[st];
			if(s[index]<=s[i]) break;
			st--;
		}
		if(st==-1) PD[i]=0;
		else PD[i]=i-index;
		sti[++st]=i;
	}
	free(sti);
	return;
}

void getFailure(dtype *x, int m, int *pd, int *failure){
	int len=0;
	failure[0]=0;
	int i;
	int npd;
	for(i=1;i<m;i++){
		npd=pd[i];
//		while(len!=0){
		while(1){
			if(npd>len) npd=0;
			if(npd==pd[len]) break;
			else{
				len=failure[len-1];
			}
		}
		len++;
		failure[i]=len;
	}
	return;
}

int searchct(dtype *x, int m, dtype *y, int n) {
	int count=0;
	int first=0;
	int ele=0;
	int q = 1;
	while(q<m) q=(q<<1);
	int qM1 = q-1;
	int *dqi = (int *) malloc(q*sizeof(int));

	int *pd = (int *) malloc(m*sizeof(int));
	int *failure = (int *) malloc(m*sizeof(int));
	getPD(x, m, pd);
	getFailure(x, m, pd, failure);
	int len=0;
	int i;
	int mM1=m-1;
	for(i=0;i<n;i++){
		int index;
		int npd;
		while(ele!=0){
			index=dqi[(first+ele-1)&qM1];
			if(y[index]<=y[i]) break;
			ele--;
		}
		npd = (ele>0)?i-index:0;
		while(true){
			if(npd>len) npd=0;
			if(npd == pd[len]) break;
			else{
				len = failure[len-1];
				while(ele!=0){
					if(dqi[first&qM1] < i-len){
						first++;
						ele--;
					}
					else break;
				}
			}
		}
		len++;
		if(len == m){
			OUTPUT2(i-m+1);
			len=failure[len-1];
		}

		while(ele!=0){
			if(dqi[first&qM1] <= i-len){
				first++;
				ele--;
			}
			else break;
		}
		dqi[(first+ele)&qM1]=i;
		ele++;
	}
	
	free(pd);
	free(failure);
	free(dqi);
   return count;
}


