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

#if isINT
inline int simd_comp(dtype *y, int i){
//	__m256i am = _mm256_loadu_si256((__m256i *)(y+i));
//	__m256i bm = _mm256_loadu_si256((__m256i *)(y+i+1));
//	__m256i rm = _mm256_cmpgt_epi32(am, bm);
//	return _mm256_movemask_ps(_mm256_castsi256_ps(rm));

	__m128i am = _mm_loadu_si128((__m128i *)(y+i));
	__m128i bm = _mm_loadu_si128((__m128i *)(y+i+1));
	__m128i rm = _mm_cmpgt_epi32(am, bm);
	 am = _mm_loadu_si128((__m128i *)(y+i+4));
	 bm = _mm_loadu_si128((__m128i *)(y+i+5));
	 __m128i rm2 = _mm_cmpgt_epi32(am, bm);
	return (_mm_movemask_ps(_mm_castsi128_ps(rm))) + (_mm_movemask_ps(_mm_castsi128_ps(rm2))<<4);

//	return ((y[i]>y[i+1])<<0)+((y[i+1]>y[i+2])<<1)+((y[i+2]>y[i+3])<<2)+((y[i+3]>y[i+4])<<3)+
//			((y[i+4]>y[i+5])<<4)+((y[i+5]>y[i+6])<<5)+((y[i+6]>y[i+7])<<6)+((y[i+7]>y[i+8])<<7);
}
#else
inline int simd_comp(dtype *y, int i){
//	__m256 am = _mm256_loadu_ps(y+i);
//	__m256 bm = _mm256_loadu_ps(y+i+1);
//	__m256 rm = _mm256_cmp_ps(am, bm, 13);
//	return _mm256_movemask_ps(rm);

//	__m128 am = _mm_loadu_ps(y+i);
//	__m128 bm = _mm_loadu_ps(y+i+1);
//	__m128 rm = _mm_cmpgt_ps(am, bm);
//	 am = _mm_loadu_ps(y+i+4);
//	 bm = _mm_loadu_ps(y+i+5);
//	 __m128 rm2 = _mm_cmpgt_ps(am, bm);
//	return (_mm_movemask_ps(rm)) + (_mm_movemask_ps(rm2)<<4);

	return ((y[i]>y[i+1])<<0)+((y[i+1]>y[i+2])<<1)+((y[i+2]>y[i+3])<<2)+((y[i+3]>y[i+4])<<3)+
			((y[i+4]>y[i+5])<<4)+((y[i+5]>y[i+6])<<5)+((y[i+6]>y[i+7])<<6)+((y[i+7]>y[i+8])<<7);
}
#endif

int search(unsigned char *x, int m, dtype *y, int n) {
   /* Preprocessing */
	const int q=8;
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
   si[0]=(x[0]<<0)+(x[1]<<1)+(x[2]<<2)+(x[3]<<3)+(x[4]<<4)+(x[5]<<5)+(x[6]<<6)+(x[7]<<7);
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


