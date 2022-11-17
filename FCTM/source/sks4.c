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
#include "include/main.h"
#include <emmintrin.h>
#include <immintrin.h>
#include <xmmintrin.h>

inline int simd_comp(dtype *y, int i){
	return ((y[i]>y[i+1]))+((y[i+1]>y[i+2])<<1)+((y[i+2]>y[i+3])<<2)+((y[i+3]>y[i+4])<<3);
}


const int q=4;
int delta[1<<q];
int si[XSIZE];
int nex[XSIZE];

int search(unsigned char *x, int m, dtype *y, int n) {
   /* Preprocessing */
	int count=0;
	if(m<q) return -1;
	si[0]=0;
	memset(delta, -1, sizeof(int)*(1<<q));
	memset(nex, -1, sizeof(int)*XSIZE);
	for(int i=0;i<q;i++){
		si[0]+=(int(x[i])<<i);
	}
   for(int i=0;i<m-q+1;i++){
	   nex[i] = delta[si[i]];
	   delta[si[i]] = i;
	   if(i<m-q) si[i+1] = (si[i]>>1)+(int(x[i+q])<<(q-1));
   }
  /* Searching */
   int tem;
   int mMq1=m-q+1;
   int nMq1=n-q+1;
   for(int i=m-q;i<nMq1;i+=mMq1){
	   tem=delta[simd_comp(y, i)];
	   while(tem>-1){
		   int j=mMq1;
		   for(j=0;j<mMq1;j+=q){
			   if(simd_comp(y, i-tem+j)!=si[j]) break;
		   }
		   if(j>=mMq1){
			   OUTPUT(i-tem);
		   }
		   tem=nex[tem];
	   }
   }
   return count;
}


