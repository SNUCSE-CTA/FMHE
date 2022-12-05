#include "AuxData.h"

using namespace std;

AuxData::AuxData(){
}

AuxData::~AuxData() {
}

void AuxData::Preprocessing(uchar *text,uint k, ulong startp[], int min, int _m) {
    m = _m;
    Alphabet = 4;
    rank['A'] = 0; rank['C'] = 1; rank['G'] = 2; rank['T'] = 3;
    B = (int)(log(2*m*k) / log(Alphabet) );
    if(m < B) B = m;
    ShiftBits = (unsigned int)ceil(log((double)Alphabet)/ log(2.0));

    cerr << "B : " << B << endl;
    cerr << "m : " << m << endl;
    BMask = 1;
    for(int i = 0; i < B*ShiftBits; i++) BMask <<= 1;
    BMask--;

    Size = (ulong) pow((double)Alphabet, B);
    cerr << "2km : " << Size << endl; 
    qrm = new unsigned int[Size];
    if( min < m ){
        qList = new LinkedList [Size];
    }
    if( min < B ){
        Bprefixp = new int [B];
        Bprefixp[0] = 0; int p = 1;
        for(int i = 1; i < B; i++){
            p *= Alphabet;
            Bprefixp[i] = Bprefixp[i-1] + p;
        }

        Bprefix = new struct tree_node * [Size];
        for(int j = 0; j < Size; j++){
            Bprefix[j] = NULL;
        }

    }
    for (size_t i = 0; i < Size; ++i) {
        qrm[i] = B - 1;
    }
    uint* hash = new uint[k];
    unsigned int* T = new unsigned int[Size/32];    // We use bit vector $T$
    for(int i = 0; i < Size/32; i++) T[i] = 0;
    for(int i = 0; i < k; i++){
        if( startp[i+1] - startp[i] >= B ){
            uint h = 0;
            for(int r = B; r >= 1; --r){
                h <<= ShiftBits;
                h += rank[decode(text, startp[i+1] - r)];
            }
            T[h/32] |= (1<<h%32);
        }
    }
    m2; if( m >= min ) m2 = m; else m2 = min;   // m2 = m' in our algorithm
    B2; if( B >= min ) B2 = B; else B2 = min;   // B2 = B'
    for (int q = m; q >= B; --q) {
        for (int j = 0; j < k; ++j) {  
            int len = startp[j+1] - startp[j];
            if( len >= m2 ){
                // computing hash for only case1 strings
                if( q == m ){    			
                    hash[j] = 0;
                    for (int r = B-1; r >= 0; --r) {
                        hash[j] <<= ShiftBits;
                        hash[j] += rank[decode(text, startp[j]+q-r-1)];
                    }
                }
                else{
                    hash[j] = hash[j] >> ShiftBits;
                    unsigned int temp = rank[decode(text, startp[j]+q-B)];
                    temp <<= ((B-1)*ShiftBits);
                    hash[j] += temp;
                }

                int h = hash[j];
                if( qrm[h] < q ){
                    qrm[h] = q;
                }
            }
            if( B2 <= len ){
                if( B2 <= len && len < m ){
                    if( q == len ){
                        hash[j] = 0;
                        for(int r = B-1; r >= 0; --r){
                            hash[j] <<= ShiftBits;
                            hash[j] += rank[decode(text, startp[j]+q-r-1)];
                        }
                    }
                    else if( q < len ){
                        hash[j] = hash[j] >> ShiftBits;
                        unsigned int temp = rank[decode(text, startp[j]+q-B)];
                        temp <<= ((B-1)*ShiftBits);
                        hash[j] += temp;
                    }
                    else continue;
                }
                int h = hash[j];
                if( B2 <= q && q < m && (T[h/32] & (1<<(h%32)))){
                    if( qList[h].empty() ){
                        qList[h].push_front( (unsigned int) q );
                    }
                    if( !qList[h].empty() && qList[h].front() != q ){
                        qList[h].push_front( (unsigned int) q );
                    }
                }
            }
        }
    }
    /*double aver_q = 0.0;
    for(int i = 0; i < Size; i++){
        aver_q += qrm[i];
    }
    cerr << "Average of qrm : " << aver_q / Size << endl;*/

    delete [] hash;
    delete [] T;
}

#define NUMLOCK 1024

void AuxData::Preprocessing_par(uchar *text, uint k, ulong startp[], int min, int threads, int _m) {
    m = _m;
    Alphabet = 4;
    rank['A'] = 0; rank['C'] = 1; rank['G'] = 2; rank['T'] = 3;
    B = log(2*m*k) / log(Alphabet);
    if(m < B) B = m;
    ShiftBits = (unsigned int)ceil(log((double)Alphabet) / log(2.0));

    cerr << "B : " << B << endl;
    cerr << "m : " << m << endl;
    BMask = 1;
    for (int i = 0; i < B*ShiftBits; i++) BMask <<= 1;
    BMask--;

    Size = (ulong) pow(pow(2.0, ShiftBits), B);

    qrm = (unsigned int*) memalign(64, Size*sizeof(unsigned int));

    unsigned int* T;
    if( min < m ){
        qList = new LinkedList [Size];
        T = (unsigned int*)( memalign(64, Size*sizeof(unsigned int)));
    }

    if( B > min ){
        Bprefixp = new int[B];
        Bprefixp[0] = 0; int p = 1;
        for(int i = 1; i < B; i++){
            p *= Alphabet;
            Bprefixp[i] = Bprefixp[i-1] + p;
        }   
        Bprefix = new struct tree_node *[Size];
        for (int j = 0; j < Size; j++){
            Bprefix[j] = NULL;
        }
    }
#pragma omp parallel for
    for (size_t i = 0; i < Size; ++i) {
        qrm[i] = B - 1;
    }

    omp_lock_t lock[NUMLOCK];

    for(int i = 0; i < NUMLOCK; i++){
        omp_init_lock(&lock[i]);
    }

    uint* hash = (unsigned int*) memalign(64, k*sizeof(unsigned int));
    if(min < m){
#pragma omp parallel for
        for(int i = 0; i < Size; i++){
            T[i] = 0;
        }

#pragma imp parallel for
        for(int i = 0; i < k; i++){
            int h = 0;
            for(int r = B; r >= 1; r--){
                h <<= ShiftBits;
                h += rank[decode(text, startp[i+1]-r)];
            }
            omp_set_lock(&lock[h%NUMLOCK]);
            T[h] = 1;
            omp_unset_lock(&lock[h%NUMLOCK]);
        }
    }
    m2; if(m >= min) m2 = m; else m2 = min;
    B2; if(B >= min) B2 = B; else B2 = min;
    for (int q = m; q >= B; --q) {
#pragma omp parallel for
        for (int j = 0; j < k; ++j) {  // loop through patterns
            int len = startp[j+1] - startp[j];
            if( len >= m2 ){
                if (q == m){
                    hash[j] = 0;
                    for (int r = B - 1; r >= 0; --r) {
                        hash[j] <<= ShiftBits;
                        hash[j] += rank[decode(text, startp[j] + q - r - 1)];
                    }
                }
                else{
                    hash[j] = hash[j] >> ShiftBits;
                    unsigned int temp = rank[decode(text, startp[j] + q - B)];
                    temp <<= ((B - 1)*ShiftBits);
                    hash[j] += temp;
                }
                int h = hash[j];

                omp_set_lock(&lock[h%NUMLOCK]);
                if (qrm[h] < q){
                    qrm[h] = q;
                }
                omp_unset_lock(&lock[h%NUMLOCK]);
            }
            if( B2 <= len ){

                if( B2 <= len && len < m ){
                    if( q == len ){
                        hash[j] = 0;
                        for(int r = B-1; r >= 0; r--){
                            hash[j] <<= ShiftBits;
                            hash[j] += rank[decode(text, startp[j]+q-r-1)];
                        }
                    }
                    else if( q < len ){
                        hash[j] = hash[j] >> ShiftBits;
                        unsigned int temp = rank[decode(text, startp[j]+q-B)];
                        temp <<= ((B-1)*ShiftBits);
                        hash[j] += temp;
                    }
                    else continue;
                }
                int h = hash[j];
                if( B2 <= q && q < m && T[h] == 1 ){
                    omp_set_lock(&lock[h%NUMLOCK]);
                    if ( qList[h].empty() ){
                        qList[h].push_front((unsigned int)q);
                    }
                    if ( !qList[h].empty() && qList[h].front() != q){
                        qList[h].push_front((unsigned int)q);
                    }
                    omp_unset_lock(&lock[h%NUMLOCK]);

                }

            }

        }
    }
    for(int i = 0; i < NUMLOCK; i++){
        omp_destroy_lock(&lock[i]);
    }
}
