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
const int bsize=12;
#if isINT
#if isLARGE
inline int simd_comp(dtype *y, int i){
	__m128i am = _mm_loadu_si128((__m128i *)(y+i));
	__m128i bm = _mm_loadu_si128((__m128i *)(y+i+1));
	__m128i rm = _mm_cmpgt_epi32(am, bm);
	return _mm_movemask_ps(_mm_castsi128_ps(rm));
}
const int q=4;
#define simd_comp_block(y,i) (simd_comp(y,i)+(simd_comp(y,(i)+4)<<4)+(simd_comp(y,(i)+8)<<8))
#define si_block(i) (si[i]+(si[(i)+4]<<4)+(si[(i)+8]<<8))
#else
inline int simd_comp(dtype *y, int i){
	__m128i am = _mm_loadu_si128((__m128i *)(y+i));
	__m128i bm = _mm_loadu_si128((__m128i *)(y+i+1));
	__m128i rm = _mm_cmpgt_epi8(am, bm);
	return _mm_movemask_epi8(rm) & 0xfff;
}
const int q=bsize;
#define simd_comp_block(y,i) simd_comp(y,i)
#define si_block(i) si[i]
#endif
#else
inline int simd_comp(dtype *y, int i){
	__m128 am = _mm_loadu_ps(y+i);
	__m128 bm = _mm_loadu_ps(y+i+1);
	__m128 rm = _mm_cmpgt_ps(am, bm);
	return _mm_movemask_ps(rm);
}
const int q=4;
#define simd_comp_block(y,i) (simd_comp(y,i)+(simd_comp(y,(i)+4)<<4)+(simd_comp(y,(i)+8)<<8))
#define si_block(i) (si[i]+(si[(i)+4]<<4)+(si[(i)+8]<<8))
#endif

int delta[1<<bsize];
int si[XSIZE];
int nex[XSIZE];

int search(unsigned char *x, int m, dtype *y, int n) {
   /* Preprocessing */
	int count=0;
	if(m<bsize) return -1;
	si[0]=0;
	memset(delta, -1, sizeof(int)*(1<<bsize));
	memset(nex, -1, sizeof(int)*XSIZE);
	for(int i=0;i<q;i++){
		si[0]+=(int(x[i])<<i);
	}
	for(int i=0;i<m-q;i++){
		si[i+1]=(si[i]>>1)+(int(x[i+q])<<(q-1));
	}
   for(int i=0;i<m-bsize+1;i++){
	   nex[i] = delta[si_block(i)];
	   delta[si_block(i)] = i;
   }
  /* Searching */
   int tem;
   int mMb1=m-bsize+1;
   int nMb1=n-bsize+1;
   int mMq1=m-q+1;
   for(int i=m-bsize;i<nMb1;i+=mMb1){
	   tem=delta[simd_comp_block(y, i)];
	   while(tem>-1){
		   int j=0;
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


