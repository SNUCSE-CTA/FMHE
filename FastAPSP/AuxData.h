#ifndef _AUXDATA_
#define _AUXDATA_

#include <vector>
#include <cmath>
#include "Ctree.h"
#include <list>
#include <vector>
#include <malloc.h>
#include "LinkedList.h"

using namespace std;

class AuxData {
public:
	unsigned int k, m, B, m2, B2;
    unsigned int *qrm;
    LinkedList *qList;
    struct tree_node ** Bprefix;
    
    unsigned int BMask;
    int* Bprefixp;
    
    unsigned int* prefixH;

    unsigned int Alphabet;
    unsigned int ShiftBits;
    ulong Size;
    unsigned int rank[128];

    AuxData();
	~AuxData();
	void Preprocessing(uchar *text, uint k, ulong startp[], int min, int _m);
	void Preprocessing_par(uchar *text, uint k, ulong startp[], int min, int threads, int _m);

private:
};

#endif
