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
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
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

int verification(dtype *x, int m, dtype *y, int n, int *candidates, int ele, int *res){
    int count=0;
    int *pd=(int *) malloc(m*sizeof(int));
    getPD(x, m, pd);
    int st=-1;
    dtype *stv = (dtype *) malloc(m*sizeof(dtype));
    int *sti = (int *) malloc(m*sizeof(int));
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


int search(unsigned char* p, int m, dtype* t, int n);

int *candidates;

int main(int argc, char *argv[])
{
	std::default_random_engine generator;
    unsigned char *p, *t;
	dtype *rawp, *rawt;
	int m, n;
	m=std::stoi(argv[1]);
	n=10000000;
	std::uniform_int_distribution<int> distribution(0, n-m);
	auto myrand = std::bind(distribution, generator);
	if(argc<2){
		printf("no arg\n");
		return 0;
	}
	std::string filename;
	if(isINT) filename = "data/intd";
	else filename = "data/floatd";
	fil_timer=(TIMER *) malloc(sizeof(TIMER));
	sea_timer=(TIMER *) malloc(sizeof(TIMER));
	pre_timer=(TIMER *) malloc(sizeof(TIMER));
	ver_timer=(TIMER *) malloc(sizeof(TIMER));
	FILE *f;
	if((f = fopen(filename.c_str(), "r"))==NULL){
		printf("Cannot open the file\n");
		return 0;
	}
//	if((f = fopen(TFILE, "r")) == NULL){
//		printf("Cannot open the file\n");
//		return 0;
//	}
	int i;
//	rawp = (dtype *) malloc(sizeof(dtype) * m);
	rawt = (dtype *) malloc(sizeof(dtype) * (n+m));
	p = (unsigned char *) malloc(sizeof(unsigned char) * m);
	if(isINT){
//		for(i=0;i<m;i++){
//			fscanf(f, "%d", &rawp[i]);
//		}
		for(i=0;i<n;i++){
			fscanf(f, "%d", &rawt[i]);
		}
	}
	else{
		if(isFLOAT){
//			for(i=0;i<m;i++){
//				fscanf(f, "%f", &rawp[i]);
//			}
			for(i=0;i<n;i++){
				fscanf(f, "%f", &rawt[i]);
			}
		}
		else{
//			for(i=0;i<m;i++){
//				fscanf(f, "%lf", &rawp[i]);
//			}
			for(i=0;i<n;i++){
				fscanf(f, "%lf", &rawt[i]);
			}
		}
	}
	double totaltime=0;
	candidates= (int *) malloc(sizeof(int) * n);
	int *res = (int *) malloc(sizeof(int) * n);
	for(int num=0;num<RNUM;num++){
		rawp=rawt+myrand();
		BEGIN_FILTERING
		for(i=0;i<m-1;i++){
			p[i] = (rawp[i]<=rawp[i+1]);
		}
		dtype diff=rawt[n-1]-rawp[0];
		for(i=1;i<m;i++){
			rawt[n-1+i]=rawp[i]+diff;
		}

		int occ = search(p, m-1, rawt, n-1);
		END_FILTERING
//		printf("found %d candidates\n", occ);
//		printf("occ: ");
//		for(i=0;i<occ;i++){
//			printf("%d ", candidates[i]);
//		}
//		printf("\n");
	
		BEGIN_VERIFY
		int count = verification(rawp, m, rawt, n, candidates, occ, res);
		END_VERIFY
		totaltime += timer_elapsed(fil_timer)+timer_elapsed(ver_timer);

//		printf("res: ");
//		for(i=0;i<count;i++){
//			printf("%d ", res[i]);
//		}
//		printf("\n");
//		printf("found %d occurrences\n", count);
//		printf("FILTER: %lf, VERIFY: %lf, SUM: %lf\n", timer_elapsed(fil_timer), timer_elapsed(ver_timer), timer_elapsed(fil_timer) + timer_elapsed(ver_timer));
//		BEGIN_SEARCHING
//		occ = searchct(rawp, m, rawt, n);
//		END_SEARCHING

//		printf("found %d occurrences\n", occ);
//		printf("CT: %lf\n", timer_elapsed(sea_timer));

//		if(count == occ){
//			for(i=0;i<occ;i++){
//				if(candidates[i]!=res[i]) break;
//			}
//			if(i==occ) printf("test success!\n");
//			else printf("test failed!\n");
//		}
//		else printf("test failed!\n");
	}

	printf("%lf", totaltime);
	free(candidates);
	free(res);
	free(p);
	free(t);
	free(rawt);
	free(fil_timer);
	free(sea_timer);
	free(pre_timer);
	free(ver_timer);
	return 0;
}
