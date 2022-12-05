#ifndef _CTREE_
#define _CTREE_

#include <omp.h>

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long ulong;

#define SEPERATOR '\n'

#include "AuxData.h"
#include "oList.h"

struct tree_node{
	int pos;
	int from;
	int to;
    int sdepth;
    int rem[5];
	struct tree_node *ptr[5];
};

struct stack_node{
	int value;
	struct stack_node *next;
};

class AuxData; 

struct tree_node *create_tree(uchar *text,uint k,ulong n,ulong startp[],uint sorted[]);
struct tree_node *create_tree_modified(uchar *text,uint k,ulong n,ulong startp[],struct stack_node **stacks);

struct tree_node *create_node(int pos,int from,int to);
struct stack_node* create_stack_node(uint i,struct stack_node*ptr);

int get_index(uchar c);
void display_node(struct tree_node* p);
ulong find_pos(ulong startpos[],uint sorted[], int i,ulong k,ulong n);
ulong find_pos_modified(ulong startpos[], int i,ulong k,ulong n);
void traverse_tree_prefix(struct tree_node *ptr, uchar *T, ulong startpos[], uint sorted[], int curpos, int hash, AuxData* wm);
int traverse_tree_modified(struct tree_node *ptr,uint sorted [],int *counter,struct stack_node **stacks, int depth, struct tree_node **sptr, struct tree_node *parent, int sdepth);
void updateNodeRange(struct tree_node *ptr,int i);
void display_tree(struct tree_node* p);
uchar get_char(int c);
void remove_index(struct tree_node *ptr,int i);
void find_all_pairs(struct tree_node *ptr,uchar *T,ulong N,ulong K,ulong startpos[],uint sorted[],
			int output,int threads,int distribution_method,int min, AuxData* ad, struct tree_node**);
void find_all_pairs_modified(struct tree_node *ptr,uchar *T,ulong N,ulong k,ulong startpos[],uint sorted[],int output);
void traversal_modified(struct tree_node *ptr,int i,int **A,struct tree_node *root,uchar *T,int curpos,int curpos2,ulong startpos[],uint sorted[],ulong N,ulong k);
void traversal_helper_modified(struct tree_node* ptr2,int i2,struct tree_node *ptr1,int i1,uchar * T,int curpos,int curpos2,uint k,int N,ulong startpos[],uint sorted[],int **A);
void traversal1(struct tree_node *ptr,int **A,struct tree_node *root,uchar *T,int curpos,int curpos2,ulong startpos[],uint sorted[],ulong N,ulong k,int output);
//void traversal1_no_output(struct tree_node *ptr,int **A,struct tree_node *root,uchar *T,int curpos,int curpos2,ulong startpos[],uint sorted[],ulong N,ulong k);
void traversal2(struct tree_node* ptr2,int i2,struct tree_node *ptr1,int i1,uchar * T,int curpos,int curpos2,uint k,int N,ulong startpos[],uint sorted[],int **A,int output);
//void traversal2_no_output(struct tree_node* ptr2,int i2,struct tree_node *ptr1,int i1,uchar * T,int curpos,int curpos2,uint k,int N,ulong startpos[],uint sorted[],int **A);

/*
    CT's root : struct tree_node* ptr
    sorted : sorted
    qrm : ad->qrm
    qList : ad->qList
    Bprefix : ad->Bprefix
*/
void matching_step(uchar *T,ulong startpos[],uint sorted[],int **A,struct  tree_node* ptr,ulong k,ulong N,int output,int min, AuxData* ad, struct tree_node**);
void matching_step_par(int threads,uchar *T,ulong startpos[],uint sorted[],int **A,struct  tree_node* ptr,ulong k,ulong N,int output,int min, AuxData* ad, oList*, struct tree_node**);

void insert_output(int, unsigned int*, int , struct tree_node*, int, oList*, int**);

char decode(uchar *final,ulong pos);
#endif
