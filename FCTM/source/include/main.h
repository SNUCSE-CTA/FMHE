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
 */

#include "timer.h"
#include "ct.h"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
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
dtype t[10001024];
dtype *rawp;
int m, n;
int parent[XSIZE];
/*
int verification(dtype *x, int m, dtype *y, int n, int *candidates, int ele, int *res){
	int count=0;
	int *pd = (int *) malloc(m*sizeof(int));
	int *rpd = (int *) malloc(m*sizeof(int));
	getPD(x, m, pd);
	getRPD(x, m, rpd);
	for(int i=0;i<m;i++){
		if(pd[i]==0) pd[i]=i-rpd[i];
		else if(rpd[i]==0) pd[i]=i-pd[i];
		else if(x[i-pd[i]]>x[i-rpd[i]]) pd[i]=i-pd[i];
		else pd[i]=i-rpd[i];
	}
	while(ele>0&&candidates[ele-1]>n-m) ele--; //2개씩 묶을 경우 초과가능?
	for(int e=0;e<ele;e++){
		dtype *s = y + candidates[e];
		int i;
		for(i=0;i<m;i++){
			//if(s[i]<s[pd[i]]||(s[i]==s[pd[i]]&&i<pd[i])) break;
			if(s[i]<s[pd[i]]) break;
			if(s[i]==s[pd[i]]&&i<pd[i]) break;
		}
		if(i==m) res[count++]=candidates[e];
	}
	free(pd);
	free(rpd);
	return count;
}

int verification2(dtype *x, int m, dtype *y, int n, int *candidates, int ele, int *res){
    int count=0;
    int *pd=(int *) malloc(m*sizeof(int));
    getPD(x, m, pd);
    int st=-1;
    dtype *stv = (dtype *) malloc(m*sizeof(dtype));
    int *sti = (int *) malloc(m*sizeof(int));
	while(ele>0&&candidates[ele-1]>n-m) ele--;
	//if(ele>n-m) ele--;
    for(int e=0;e<ele;e++){
        dtype *s = y+candidates[e];
        st=-1;
        int i;
        for(i=0;i<m;i++){
            dtype value;
            int index;
            int npd;
            while(st!=-1){
                value=stv[st];
                index=sti[st];
                if(value<=s[i]) break;
                st--;
            }
            if(st==-1) npd=0;
            else npd=i-index;
            if(npd!=pd[i]) break;
            stv[++st]=s[i];
            sti[st]=i;
        }
        if(i==m){
            res[count++]=candidates[e];
        }
    }
	free(pd);
	free(stv);
	free(sti);
    return count;
}

*/
int search(unsigned char *p, int m, dtype* t, int n);

int *candidates;

int main(int argc, char *argv[])
{
	std::default_random_engine generator;
    unsigned char p[XSIZE];
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
	for(i=0;i<n;i++){
#if isINT
	#if isLARGE
		fscanf(f, "%d", &t[i]);
	#else
		fscanf(f, "%hhd", &t[i]);
	#endif
#else
	#if isLARGE
		fscanf(f, "%lf", &t[i]);
	#else
		fscanf(f, "%f", &t[i]);
	#endif
#endif
	}
	double totaltime=0, verifytime=0, filtertime=0, searchtime=0;
	int occ, count;
	int total_occ=0;
	for(int num=0;num<RNUM;num++){
		rawp=t+myrand();
		BEGIN_FILTERING
		for(i=0;i<m-1;i++){
			p[i] = (rawp[i]>rawp[i+1]);
		}
		for(i=0;i<m;i++){
			t[n+i]=rawp[i];
		}
		getParent(rawp, m, parent);
		occ = search(p, m-1, t, n);
		total_occ+=occ;
		END_FILTERING
//		printf("found %d candidates\n", occ);
//		printf("occ: ");
//		for(i=0;i<occ;i++){
//			printf("%d ", candidates[i]);
//		}
//		printf("\n");
		BEGIN_VERIFY
//		int count = verification(rawp, m, t, n, candidates, occ, res);
		END_VERIFY
		filtertime += timer_elapsed(fil_timer);
//		verifytime += timer_elapsed(ver_timer);
//		searchtime += timer_elapsed(sea_timer);
//		total_occ+=occ;
///		printf("found %d occurrences\n", count);
//		BEGIN_SEARCHING
//		count = searchct(rawp, m, t, n);
//		END_SEARCHING
//
//		printf("found %d occurrences\n", occ);
//		printf("%d ", occ);
//		printf("occ: ");
//		for(i=0;i<occ;i++){
//			printf("%d ", candidates[i]);
//		}
//		printf("\n");

//		printf("CT: %lf\n", timer_elapsed(sea_timer));

//		if(count == occ){
//			printf("success!\n");
//		}
//		else printf("failed!\n");
	}
	totaltime = filtertime;
//	fprintf(stderr, "%lf ", searchtime);
	fprintf(stderr, "%d ", total_occ);
	printf("%.2lf", totaltime);
	free(fil_timer);
	free(sea_timer);
	free(pre_timer);
	free(ver_timer);
	return 0;
}
