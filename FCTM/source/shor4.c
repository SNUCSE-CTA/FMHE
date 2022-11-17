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
/*
#if isINT
#define simd_comp(y, i) _mm_movemask_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(_mm_loadu_si128((__m128i *)((y)+(i))), _mm_loadu_si128((__m128i *)((y)+(i)+1)))))
#else
#define simd_comp(y, i) _mm_movemask_ps(_mm_cmpgt_ps(_mm_loadu_ps((y)+(i)), _mm_loadu_ps((y)+(i)+1)))
#endif
*/

#if isINT
#if isLARGE
inline int simd_comp(dtype *y, int i){
	__m128i am = _mm_loadu_si128((__m128i *)(y+i));
	__m128i bm = _mm_loadu_si128((__m128i *)(y+i+1));
	__m128i rm = _mm_cmpgt_epi32(am, bm);
	return _mm_movemask_ps(_mm_castsi128_ps(rm));

//	__m128i rm = _mm_cmplt_epi32(am, bm);
//	return _mm_movemask_ps(_mm_castsi128_ps(~rm));

//	return ((y[i]>y[i+1]))+((y[i+1]>y[i+2])<<1)+((y[i+2]>y[i+3])<<2)+((y[i+3]>y[i+4])<<3);
}
#else

inline int simd_comp(dtype *y, int i){
	__m128i am = _mm_loadu_si128((__m128i *)(y+i));
	__m128i bm = _mm_loadu_si128((__m128i *)(y+i+1));
	__m128i rm = _mm_cmpgt_epi8(am, bm);
	return _mm_movemask_epi8(rm) & 0xf;
}
#endif
#else
inline int simd_comp(dtype *y, int i){
	__m128 am = _mm_loadu_ps(y+i);
	__m128 bm = _mm_loadu_ps(y+i+1);
	__m128 rm = _mm_cmpgt_ps(am,bm);
	return _mm_movemask_ps(rm);

//	return ((y[i]>y[i+1]))+((y[i+1]>y[i+2])<<1)+((y[i+2]>y[i+3])<<2)+((y[i+3]>y[i+4])<<3);
}
#endif

const int q=4;
int delta[1<<q];
int si[XSIZE];

int search(unsigned char *x, int m, dtype *y, int n) {
   /* Preprocessing */
	if(m<q) return -1;
   int count=0;
   for(int i=0;i<16;i++) delta[i]=m;
   for(int i=0;i<8;i++) delta[(i<<0)+(x[0]<<3)]=m-1;
   for(int i=0;i<4;i++) delta[(i<<0)+(x[0]<<2)+(x[1]<<3)]=m-2;
   for(int i=0;i<2;i++) delta[(i<<0)+(x[0]<<1)+(x[1]<<2)+(x[2]<<3)]=m-3;
   for(int i=0;i<m-3;i++){
	   si[i]=(x[i]<<0)+(x[i+1]<<1)+(x[i+2]<<2)+(x[i+3]<<3);
	   delta[si[i]]= m-i-4;
   }
  /* Searching */
   int i=m-4;
   int j, k;
   const int mMq=m-q;
   while (true) {
	   do {
		   k=delta[simd_comp(y, i)];
		   i+=k;
	   } while (k>0);
	   for(j=0;j<mMq;j+=4){
		   if(simd_comp(y, i-mMq+j)!=si[j]) break;
	   }
	   if(j>=mMq){
		   if(i>n-4) break;
		   OUTPUT(i-mMq);
	   }
	   i++;
   }
   return count;
}


