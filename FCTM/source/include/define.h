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

#ifndef __DEF_H__
#define __DEF_H__


#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define FALSE      0
#define TRUE       1
//#define TSIZE		658795
#define TSIZE		10000000
#define XSIZE       1024			//maximal length of the pattern
#define WSIZE  	    256				//greater int value fitting in a computer word
//#define SIGMA       256				//constant alphabet size
#define SIGMA 2
//#define ALPHA       256				//constant alphabet size
//#define ASIZE		256				//constant alphabet size
#define UNDEFINED       -1
#define HALFDEFINED     -2
#define WORD	    32				//computer word size (in bit)
//#define OUTPUT(j)   {candidates[count++]=(j);}

#define OUTPUT(j) {verification(rawp,y,j,parent,count);}
#define OUTPUT2(j) {count++;}
//#define OUTPUT2(j) {candidates[count++]=j;}
#define RNUM 100

#define isINT 1 //10:signed char 11:int 20:float 21:double
#define isLARGE 0

#if isINT
	#if isLARGE
		typedef int dtype;
		//#define TFILE "data/reald"
		#define TFILE "data/intd"
	#else
		typedef signed char dtype;
		#define TFILE "data/chad"
	#endif
#else
	#if isLARGE
		typedef double dtype
		#define TFILE "data/doud"
	#else
		typedef float dtype;
		#define TFILE "data/flod"
	#endif
#endif

#define gety(y, i) ((y)[(i)]>(y)[(i)+1])
//#define gety(y, i) (((y)[(i)]<=(y)[(i)+1])<<1)+((y)[(i)+1]<=(y)[(i)+2])
//#define gety(i) (gety2((i)<<1))

//inline unsigned char gety(dtype *y, int i){
//	//i=i<<1;
//	return ((y[i]<=y[i+1])<<1)+(y[i+1]<=y[i+2]);
//}

inline bool mycmp(unsigned char *x, dtype *y, int m){
	for(int i=0;i<m;i++){
		if(x[i]!=gety(y, i)) return true;
	}
	return false;
}

inline bool mycmp(unsigned char *x, dtype *y, int m, int q){
	for(int i=0;i<m;i++){
		if(x[i]!=gety(y, q+i)) return true;
	}
	return false;
}


#endif
