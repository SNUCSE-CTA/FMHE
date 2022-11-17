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
#ifndef __CT_H
#define __CT_H

#include "define.h"
#include <cstdio>
#include <cstdlib>

//dtype pd[XSIZE];
extern int parent[XSIZE];
int id1[XSIZE], id2[XSIZE], gp1[XSIZE], gp2[XSIZE], gpm1, gpm2;//1은 부모가 앞에 있음.
extern int m, n;
extern dtype *rawp;

inline void verification2(dtype *x, dtype *y, int c, int *parent, int &count){
	dtype *s = y + c;
	int i;
	for(i=0;i<gpm1;i++){
		if(s[id1[i]]<s[gp1[i]]) return;
	}
	for(i=0;i<gpm2;i++){
		if(s[id2[i]]<=s[gp2[i]]) return;
	}
	if(c<n-m+1) count++;
}

inline void verification(dtype *x, dtype *y, int c, int *parent, int &count){
	dtype *s = y + c;
	int i;
/*	for(i=0;i<m;i++){
		if(s[i]<s[parent[i]]) break;
		if(s[i]==s[parent[i]]&&i<parent[i]) break;
	}
	if(i==m){
		if(c < n-m+1) count++;
	}
	*/
	for(i=0;i<m;i++){
		if(s[i]<s[parent[i]]) return;
		if(s[i]==s[parent[i]]&&i<parent[i]) return;
	}
	if(c<n-m+1) count++;
}

int sti[XSIZE];

void getPD(dtype *s, int n, int *PD, int *CD){
	int st=-1;
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
	return;
}

void getParent(dtype *s, int n, int *PD){
	int st=-1;
	for(int i=0;i<n;i++){
		int cd = -1;
		int index;
		while(st!=-1){
			index=sti[st];
			if(s[index]<=s[i]) break;
			st--;
			cd = index;
		}
		if(st==-1) PD[i]=i;
		else PD[i]=index;
		if(cd>-1) PD[cd]=i;
		sti[++st]=i;
	}
	gpm1=gpm2=0;
	for(int i=0;i<n;i++){
		if(PD[i]<i-1){
			id1[gpm1]=i;
			gp1[gpm1++]=PD[i];
		}
		else if(PD[i]>i+1){
			id2[gpm2]=i;
			gp2[gpm2++]=PD[i];
		}
	}
/*	for(int i=0;i<n;i++){
		printf("%d ", i);
	}
	printf("\n");
	for(int i=0;i<n;i++){
		printf("%d ", PD[i]);
	}
	printf("\n");
	for(int i=0;i<gpm;i++){
		printf("%d ", gp1[i]);
	}
	printf("\n");
	for(int i=0;i<gpm;i++){
		printf("%d ", gp2[i]);
	}
	printf("\n");
*/	
	return;
}


void getFailure(dtype *x, int m, int *pd, int *failure){
	int len=0;
	failure[0]=0;
	int i;
	int npd;
	for(i=1;i<m;i++){
		npd=pd[i];
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

int searchct(dtype *p, int m, dtype *y, int n) {
	int count=0;
	int first=0;
	int ele=0;
	int q=1;
	while(q<m) q=(q<<1);
	int qM1=q-1;
	int *dqi = (int *) malloc(q*sizeof(int));

	int *pd = (int *) malloc(m*sizeof(int));
	int *failure = (int *) malloc(m*sizeof(int));
	int cd[XSIZE];
	getPD(p, m, pd, cd);
	getFailure(p, m, pd, failure);
	int len=0;
	int i;
	int last;
	for(i=0;i<n;i++){
		int index;
		int npd;
		while(ele!=0){
			last=(first+ele-1)&qM1;
			index=dqi[last];
			if(y[index]<=y[i]) break;
			ele--;
		}
		npd = (ele>0)?i-index:0;
		while(1){
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
			OUTPUT(i-m+1);
			len=failure[len-1];
		}

		while(ele!=0){
			if(dqi[first&qM1] <= i-len){
				first++;
				ele--;
			}
			else break;
		}
		last =(first+ele)&qM1;
		dqi[last]=i;
		ele++;
	}
	
	free(pd);
	free(failure);
	free(dqi);
   return count;
}

#endif
