#include "include/define.h"
#include "include/main.h"

#define BYTE 8
#define MLEN 128

typedef struct skip_bucket
{
   //int pos;
	struct skip_bucket *next;
   int i,pos;
} SKBUCKET;

void PrecomputePatt(unsigned char *P, int m, unsigned char Patt[BYTE][MLEN], unsigned char Mask[BYTE][MLEN], unsigned int *Last)
{
	int w,i,j,k,h, plen, plast, last;
   unsigned char M[MLEN];
   k = 8;
//   plen = (int)ceil((double)m / (double)k);
	plen = (m+k-1)/k;
   plast = m%k;
   last = plen;
   for(i=0;i<MLEN;i++) M[i]=0;
   for(h=0;h<plen;h++) M[h]=~0;
   if(plast>0) M[plen-1] <<= BYTE-plast;
	for(i=0;i<8;i++)
   {
		for(h=0;h<last;h++)
      	{
				Patt[i][h]=P[h]>>i;
				Mask[i][h]=M[h]>>i;
            if(h>0)
            {
            	Patt[i][h] |= (P[h-1]<<(k-i));
            	Mask[i][h] |= (M[h-1]<<(k-i));
            }
         }
      Last[i] = last;
   	plast++;
      if(plast==8) plast=0;
      if(plast==1)
      {
      	//plast=0;
         last++;
      }
   }
}

int search(unsigned char *x, int m, dtype *y, int n) {
   int B[256],            // Bit mask constructed for each alphabet character
   	 M[8], I[8],        // Ending and Starting positions of Patt[i]
       BC[256];           // Bad character rule
   SKBUCKET *t, *Sk[256]; // bucket table
   int i, j, k, s, d, last, count, nw, minlen, pos, len, nchar, h;
   unsigned char Patt[BYTE][MLEN], Mask[BYTE][MLEN], HASH, c, A;
   unsigned int Last[BYTE];
   k = 8;
   nchar = (int)ceil((double)n / (double)k);
   count = 0;

   /* Preprocessing of the pattern */
   PrecomputePatt(x, m, Patt, Mask, Last);
   minlen = m;
   for(i=0; i<8; i++) {
      M[i] = Last[i]-1;
      if( ((i+m)%8) > 0 ) M[i]--;
      I[i] = 0; if(i>0) I[i]++;
      if (minlen>(M[i]-I[i]+1)) minlen = (M[i]-I[i]+1);
   }
   len = minlen;
   if(minlen>32) minlen=32;
   /* construction of the bucket table and bc rule*/
  	for(c=0,i=0; i<256; i++,c++)	{BC[c] = (2*minlen)-1; Sk[i]=NULL;}
   i=0; h=0;
   for(j=0;j<m-k+1;j++)
   {
      s = Last[i]-h-1; if( ((i+m)%8) > 0 ) s--;
      BC[Patt[i][h]] = minlen-1+s;
      if( s == 0)
      {
   		t = (SKBUCKET *) malloc (sizeof(SKBUCKET));
      	t->i = i;
      	t->pos = h;
      	t->next = Sk[Patt[i][h]];
      	Sk[Patt[i][h]] = t;
      }
   	i--;
      if(i<0)
      {
      	i=7;
         h++;
      }
   }

   /* construction of the NFA */
   count = 0;
   for(i=0; i<256; i++) B[i]=0;
   for(i=0; i<8; i++) {
      A = 1;
      for (j = 0; j < minlen; j++) {
      	B[Patt[i][M[i]-j]] |= A;
         A <<= 1;
      }
   }
   /* searching */
   j = len-1;
   while (j < nchar) {
      d = (B[gety(y,j)]<<1) & B[gety(y,j-1)];
      if (d != 0) {
         pos = j;
         while (d=(d+d)&B[gety(y,j-2)]) --j;
         j += minlen-2;
         if (j == pos) {
		      //match
      		t = Sk[gety(y,j)];
      		while(t!=NULL)
      		{
        		 	i = t->i;
         		pos = t->pos;
         		h = 1;
         		while(h<Last[i]-1 && Patt[i][h]==gety(y,j-pos+h)) h++;
         		if(h==Last[i]-1)
         		{
            		if( (Patt[i][0]==(gety(y,j-pos)&Mask[i][0])) &&
                		 (Patt[i][h]==(gety(y,j-pos+h)&Mask[i][h])) )  printf("%d %d\n", j,pos);//count++;//OUTPUT(j-pos);
         		}
       			t = t->next;
      		}
            ++j;
         }
      }
      else j += BC[gety(y,j+minlen-1)];
   }
	return count;
}



