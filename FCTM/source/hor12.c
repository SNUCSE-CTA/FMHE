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
	return ((y[i]>y[i+1]))+((y[i+1]>y[i+2])<<1)+((y[i+2]>y[i+3])<<2)+((y[i+3]>y[i+4])<<3)+
			((y[i+4]>y[i+5])<<4)+((y[i+5]>y[i+6])<<5)+((y[i+6]>y[i+7])<<6)+((y[i+7]>y[i+8])<<7)+
			((y[i+8]>y[i+9])<<8)+((y[i+9]>y[i+10])<<9)+((y[i+10]>y[i+11])<<10)+((y[i+11]>y[i+12])<<11);
}

const int q=12;
int delta[1<<q];
int si[XSIZE];

int search(unsigned char *x, int m, dtype *y, int n) {
   /* Preprocessing */
	const int mMq=m-q;
	if(m<q) return -1;
   int count=0;
   int k=0;
   for(int j=0;j<q;j++){
	   for(int i=0;i<(1<<(q-j));i++){
			delta[i+k]=m-j;
		}
	   k=(k>>1)+(int(x[j])<<(q-1));
	}
   si[0]=(x[0]<<0)+(x[1]<<1)+(x[2]<<2)+(x[3]<<3)+(x[4]<<4)+(x[5]<<5)+(x[6]<<6)+(x[7]<<7)+(int(x[8])<<8)+(int(x[9])<<9)+(int(x[10])<<10)+(int(x[11])<<11);
   for(int i=0;i<mMq+1;i++){
	   delta[si[i]]=mMq-i;
	   if(i!=mMq) si[i+1]=(si[i]>>1)+(int(x[i+q])<<(q-1));
   }
//   unsigned int mask=0;
//  for(int i=m-1;i>=0;i--){
//	   mask=(mask<<1)+x[i];
//  }
  /* Searching */
   int i=mMq;
   int j;
   while (true) {
	   do {
		   k=delta[simd_comp(y, i)];
		   i+=k;
	   } while (k>0);
	   for(j=0;j<mMq;j+=q){
		   if(simd_comp(y, i-mMq+j)!=si[j]) break;
		  // if(simd_comp(y, i-mMq+j)!=((mask>>(j))&0xff)) break;
		}
	   if(j>=mMq){
		   if(i>n-q) break;
		   OUTPUT(i-mMq);
	   }
	   i++;
   }
   return count;
}


