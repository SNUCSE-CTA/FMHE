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
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <random>
#include <functional>
#include <emmintrin.h>
#include <immintrin.h>
#include <xmmintrin.h>

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
	rawt = (dtype *) malloc(sizeof(dtype) * (n+m));
	if(isINT){
		if(isLARGE){
//		for(i=0;i<m;i++){
//			fscanf(f, "%d", &rawp[i]);
//		}
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
	double totaltime =0;
//	candidates=(int *) malloc(n*sizeof(int));
	int occ, total_occ=0;
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
			//printf("%d ", candidates[i]);
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

void getRPD(dtype *s, int n, int *PD){
	int st=-1;
	dtype *stv = (dtype *) malloc(n*sizeof(dtype));
	int *sti = (int *) malloc(n*sizeof(int));
	for(int i=n-1;i>=0;i--){
		dtype value;
		int index;
		while(st!=-1){
			value=stv[st];
			index=sti[st];
			if(value<s[i]) break;
			st--;
		}
		if(st==-1) PD[i]=0;
		else PD[i]=i-index;
		stv[++st]=s[i];
		sti[st]=i;
	}
	free(stv);
	free(sti);
	return;
}

void getPD(dtype *s, int n, int *PD, int *CD){
	int st=-1;
	int *sti = (int *) malloc(n*sizeof(int));
	for(int i=0;i<n;i++){
		int index;
		CD[i]=0;
		while(st!=-1){
			index=sti[st];
			if(s[index]<=s[i]) break;
			st--;
			CD[i]=i-index;
		}
		if(st==-1) PD[i]=0;
		else PD[i]=i-index;
		sti[++st]=i;
	}
	free(sti);
	return;
}

void getFailure(dtype *x, int m, int *pd, int *failure){
	int len=-1;
	failure[0]=-1;
	int i;
	int npd;
	for(i=1;i<=m;i++){
		npd=pd[i-1];
		while(len>0){
//		while(1){
			if(npd>len) npd=0;
			if(npd==pd[len]) break;
			else{
				len=failure[len];
			}
		}
		len++;
		failure[i]=len;
	}
	return;
}

void getParent(dtype *s, int n, int *PD){
	int st=-1;
	int sti[XSIZE];
	for(int i=0;i<n;i++){
		int cd=-1;
		int index=-1;
		while(st!=-1){
			index=sti[st];
			if(s[index]<=s[i]) break;
			st--;
			cd = index;
		}
		if(st==-1){
			PD[i]=i;
		}
		else PD[i]=index;
		if(cd>-1) PD[cd]=i;
		sti[++st]=i;
	}
	return;
}
int searchct(dtype *x, int m, dtype *y, int n) {
	if(m>16) return 0;
	int count=0;
	int *pd = (int *) malloc(m*sizeof(int));
	int *failure = (int *) malloc((m+1)*sizeof(int));
	int *cd = (int *) malloc(m*sizeof(int));
	int *gp = (int *) malloc(m*sizeof(int));
	getPD(x, m, pd, cd);
	getFailure(x, m, pd, failure);
//	int rpd[XSIZE], parent[XSIZE];
//	getRPD(x,m,rpd);
	getParent(x,m,gp);
	int shift1[XSIZE], shift2[XSIZE];
	int gpd1[XSIZE], gpd2[XSIZE];
	int n1=0, n2=0;
	for(int i=0;i<m;i++){
		if(gp[i]<i){
			gpd1[n1]=i-gp[i];
			shift1[n1]=i;
			n1++;
		}
		else if(gp[i]>i){
			gpd2[n2]=gp[i]-i;
			shift2[n2]=i;
			n2++;
		}
	}
	int tem;
	for(int i=0;i<n1;i++){
		for(int j=0;j<n1-i-1;j++){
			if(gpd1[j]>gpd1[j+1]){
				tem=gpd1[j];
				gpd1[j]=gpd1[j+1];
				gpd1[j+1]=tem;
				tem=shift1[j];
				shift1[j]=shift1[j+1];
				shift1[j+1]=tem;
			}
		}
	}
	for(int i=0;i<n2;i++){
		for(int j=0;j<n2-i-1;j++){
			if(gpd2[j]>gpd2[j+1]){
				tem=gpd2[j];
				gpd2[j]=gpd2[j+1];
				gpd2[j+1]=tem;
				tem=shift2[j];
				shift2[j]=shift2[j+1];
				shift2[j+1]=tem;
			}
		}
	}
	int m17=17-m;
	int i=0;
	int qprime;
	int q;
	__m128i a,b;
	for(i=0;i<=n-16;i+=m17){
		dtype *s = y+i;
		q=(1<<(m17))-1;
		a = _mm_loadu_si128((__m128i *) s);
		if(n1>0){
			b = _mm_loadu_si128((__m128i *) (s-gpd1[0]));
			qprime = _mm_movemask_epi8(~_mm_cmpgt_epi8(b,a));
			q=q&(qprime>>shift1[0]);
		}
		for(int j=1;j<n1;j++){
			if(gpd1[j-1]!=gpd1[j]){
				b = _mm_loadu_si128((__m128i *) (s-gpd1[j]));
				qprime = _mm_movemask_epi8(~_mm_cmpgt_epi8(b,a));
			}
			q=q&(qprime>>shift1[j]);
		}
		if(n2>0){
			b = _mm_loadu_si128((__m128i *) (s+gpd2[0]));
			qprime = _mm_movemask_epi8(_mm_cmpgt_epi8(a,b));
			q=q&(qprime>>shift2[0]);
		}
		for(int j=1;j<n2;j++){
			if(gpd2[j-1]!=gpd2[j]){
				b = _mm_loadu_si128((__m128i *) (s+gpd2[j]));
				qprime = _mm_movemask_epi8(_mm_cmpgt_epi8(a,b));
			}
			q=q&(qprime>>shift2[j]);
		}
		count+=__builtin_popcount(q);
	}
	if(i<n){
		q=(1<<(n-m+1-i))-1;
		__m128i a = _mm_loadu_si128((__m128i *) (y+i));
		for(int j=0;j<n1;j++){
			__m128i b = _mm_loadu_si128((__m128i *) (y+i-gpd1[j]));
			q=q&((_mm_movemask_epi8(~_mm_cmpgt_epi8(b,a)))>>shift1[j]);
		}
		for(int j=0;j<n2;j++){
			__m128i b = _mm_loadu_si128((__m128i *) (y+i+gpd2[j]));
			q=q&((_mm_movemask_epi8(_mm_cmpgt_epi8(a,b)))>>shift2[j]);
		}
		count+=__builtin_popcount(q);
	}
			
	
	free(pd);
	free(cd);
	free(failure);
	free(gp);
   return count;
}
