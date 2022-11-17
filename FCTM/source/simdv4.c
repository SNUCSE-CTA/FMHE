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
//	dtype *dy = y+i;
//	__m128i am = _mm_set_epi8(dy[15],dy[14],dy[13],dy[12],dy[11],dy[10],dy[9],dy[8],dy[7],dy[6],dy[5],dy[4],dy[3],dy[2],dy[1], dy[0]);
//	__m128i bm = _mm_set_epi8(dy[16],dy[15],dy[14],dy[13],dy[12],dy[11],dy[10],dy[9],dy[8],dy[7],dy[6],dy[5],dy[4],dy[3],dy[2],dy[1]);
	__m128i am = _mm_loadu_si128((__m128i *)(y+i));
	__m128i bm = _mm_loadu_si128((__m128i *)(y+i+1));
	__m128i rm = _mm_cmpgt_epi8(am, bm);
	return _mm_movemask_epi8(rm);

//	return (int(y[i]>y[i+1])<<0)+(int(y[i+1]>y[i+2])<<1)+(int(y[i+2]>y[i+3])<<2)+(int(y[i+3]>y[i+4])<<3)+
//			(int(y[i+4]>y[i+5])<<4)+(int(y[i+5]>y[i+6])<<5)+(int(y[i+6]>y[i+7])<<6)+(int(y[i+7]>y[i+8])<<7)+
//			(int(y[i+8]>y[i+9])<<8)+(int(y[i+9]>y[i+10])<<9)+(int(y[i+10]>y[i+11])<<10)+(int(y[i+11]>y[i+12])<<11)+
//			(int(y[i+12]>y[i+13])<<12)+(int(y[i+13]>y[i+14])<<13)+(int(y[i+14]>y[i+15])<<14)+(int(y[i+15]>y[i+16])<<15);
}

int search(unsigned char *x, int m, dtype *y, int n) {
   /* Preprocessing */
	const int q=16;
	const int mMq=m-q;
	if(m<q-1) return -1;
   int delta[1<<q];
   int count=0;
   int k=0;
   for(int j=0;j<q;j++){
	   for(int i=0;i<(1<<(q-j));i++){
			delta[i+k]=m-j;
	   }
	   k=(k>>1)+(int(x[j])<<(q-1));
	}
   int si[XSIZE];
   si[0]=0;
   for(int j=0;j<q;j++){
	   si[0]+=(int(x[j])<<j);
   }
   for(int i=0;i<m-q+1;i++){
	   delta[si[i]]=m-i-q;
	   if(i<m-q) si[i+1]=(si[i]>>1)+(int(x[i+q])<<(q-1));
   }
  /* Searching */
   int i=mMq;
   int j;
   while (i<n) {
	   do {
		   k=delta[simd_comp(y, i)];
		   i+=k;
	   } while (k>0);
	   for(j=0;j<mMq;j+=q){
		   if(simd_comp(y, i-mMq+j)!=si[j]) break;
		}
	   if(j>=mMq){
		   if(i>n-q) break;
		   OUTPUT(i-mMq);
	   }
	   i++;
   }
   return count;
}


