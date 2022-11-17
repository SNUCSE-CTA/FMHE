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
#define simd_comp(y, i) (_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(_mm_loadu_si128((__m128i *)((y)+(i))), _mm_loadu_si128((__m128i *)((y)+(i)+1)))))+(_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(_mm_loadu_si128((__m128i *)((y)+(i)+4)), _mm_loadu_si128((__m128i *)((y)+(i)+5)))))<<4))
#else
#define simd_comp(y, i) (_mm_movemask_ps(_mm_cmpgt_ps(_mm_loadu_ps((y)+(i)), _mm_loadu_ps((y)+(i)+1)))+(_mm_movemask_ps(_mm_cmpgt_ps(_mm_loadu_ps((y)+(i)+4), _mm_loadu_ps((y)+(i)+5)))<<4))
#endif
*/

#if isINT
#if isLARGE
inline int simd_comp(dtype *y, int i){
	__m128i am = _mm_loadu_si128((__m128i *)(y+i));
	__m128i bm = _mm_loadu_si128((__m128i *)(y+i+1));
	__m128i rm = _mm_cmpgt_epi32(am, bm);
	__m128i am2 = _mm_loadu_si128((__m128i *)(y+i+4));
	__m128i bm2 = _mm_loadu_si128((__m128i *)(y+i+5));
	 __m128i rm2 = _mm_cmpgt_epi32(am2, bm2);
	__m128i am3 = _mm_loadu_si128((__m128i *)(y+i+8));
	__m128i bm3 = _mm_loadu_si128((__m128i *)(y+i+9));
	 __m128i rm3 = _mm_cmpgt_epi32(am3, bm3);
	__m128i am4 = _mm_loadu_si128((__m128i *)(y+i+12));
	__m128i bm4 = _mm_loadu_si128((__m128i *)(y+i+13));
	 __m128i rm4 = _mm_cmpgt_epi32(am4, bm4);
	return (_mm_movemask_ps(_mm_castsi128_ps(rm))) + (_mm_movemask_ps(_mm_castsi128_ps(rm2))<<4)+ (_mm_movemask_ps(_mm_castsi128_ps(rm3))<<8)+ (_mm_movemask_ps(_mm_castsi128_ps(rm4))<<12);
}
#else

inline int simd_comp(dtype *y, int i){
	__m128i am = _mm_loadu_si128((__m128i *)(y+i));
	__m128i bm = _mm_loadu_si128((__m128i *)(y+i+1));
	__m128i rm = _mm_cmpgt_epi8(am, bm);
	return _mm_movemask_epi8(rm);
}
#endif
#else
inline int simd_comp(dtype *y, int i){
	__m128 am = _mm_loadu_ps(y+i);
	__m128 bm = _mm_loadu_ps(y+i+1);
	__m128 rm = _mm_cmpgt_ps(am, bm);
	__m128 am2 = _mm_loadu_ps(y+i+4);
	__m128 bm2 = _mm_loadu_ps(y+i+5);
	 __m128 rm2 = _mm_cmpgt_ps(am2, bm2);
	__m128 am3 = _mm_loadu_ps(y+i+8);
	__m128 bm3 = _mm_loadu_ps(y+i+9);
	 __m128 rm3 = _mm_cmpgt_ps(am3, bm3);
	__m128 am4 = _mm_loadu_ps(y+i+12);
	__m128 bm4 = _mm_loadu_ps(y+i+13);
	 __m128 rm4 = _mm_cmpgt_ps(am4, bm4);
	return (_mm_movemask_ps(rm)) + (_mm_movemask_ps(rm2)<<4) + (_mm_movemask_ps(rm3)<<8) + (_mm_movemask_ps(rm4)<<12);
}
#endif

const int q=16;
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
   si[0]=(x[0]<<0)+(x[1]<<1)+(x[2]<<2)+(x[3]<<3)+(x[4]<<4)+(x[5]<<5)+(x[6]<<6)+(x[7]<<7)+(int(x[8])<<8)+(int(x[9])<<9)+(int(x[10])<<10)+(int(x[11])<<11)+(int(x[12])<<12)+(int(x[13])<<13)+(int(x[14])<<14)+(int(x[15])<<15);
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


