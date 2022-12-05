#include <iostream>
#include <fstream>
#include <omp.h>
#include <cstdlib>
#include "Ctree.h"
#include <cstdio>
#include <cstring>
#include "AuxData.h"
#include <vector>

#define BUFFER 10


using namespace std;

void print_Sorted(uchar *a,uint *So);
void print(uchar *a);
bool isTerminator(uchar i);
bool isTerminatorfile(uchar i);
int choose_index(int a);
void display_choices(char *filename, int output, int threads,int min,int distribution_method);
void encode(uchar *final,ulong counter,int bitnum,char c);
ulong N;
uint K;
ulong *startS;
uint *Sorted;
struct stack_node **matched;




void Print_tree(struct tree_node *ptr, int depth){
    cout << "Depth : " << depth+ptr->pos << endl;    
    display_node(ptr);
    for(int i = 0; i < 5; i++){
        if( ptr->ptr[i] != NULL ){
            Print_tree(ptr->ptr[i], depth+1);
        }
    }
    cout << endl;
}

void APSP(char *filename, int output, int threads,int min,int distribution_method,int _c1, int _c2){
    //uint *Sorted= new uint[MAX_K];
    std::ifstream::pos_type posSize;
    std::ifstream file ((char *)filename, std::ios::in|std::ios::binary|std::ios::ate);
    ulong counter=0,bitnum=1,reminder=0,counterk=0,pos=0;
    //startS = new ulong[MAX_K];
    //struct stack_node **matched = new stack_node*[MAX_K];
    //startS[0]=0;
    uchar* text; 
    vector<unsigned int> startt;
    startt.push_back(0);

    int _m = 1000000000;
    ulong size;
    if (file.is_open())
    {
        posSize = file.tellg();
        size = posSize;
        char *memblock = new char [size/BUFFER + 1];
        text = new uchar[size/4];
        reminder = size%BUFFER;
        size = posSize;
        for(int i=0;i<BUFFER;i++){
            file.seekg (i*(size/BUFFER), std::ios::beg);
            file.read (memblock, size/BUFFER);
            for(ulong z=0;z<size/BUFFER;z++){
                if (memblock[z]!=SEPERATOR){
                    //cerr <<"encoding:"<<memblock[z]<<" pos:"<<pos<<endl;
                    encode(text,pos,bitnum,memblock[z]);
                    if (memblock[z]=='A' || memblock[z]=='C' || memblock[z]=='G' || memblock[z]=='T') {
                        counter++;
                        bitnum+=2;
                        //cerr<<"bitnum now "<<bitnum<<endl;
                        if (bitnum==9) {
                            bitnum=1;
                            pos++;
                        }
                    }
                }else {
                    startt.push_back(counter);
                    counterk++;
                }
            }
        }
        if (reminder>0){
            file.seekg (BUFFER*(size/BUFFER), std::ios::beg);
            file.read (memblock, reminder);
            for(ulong z=0;z<reminder;z++){
                if (memblock[z]!=SEPERATOR){
                    encode(text,pos,bitnum,memblock[z]);
                    if (memblock[z]=='A' || memblock[z]=='C' || memblock[z]=='G' || memblock[z]=='T') {
                        counter++;
                        bitnum+=2;
                        if (bitnum==9) {
                            bitnum=1;
                            pos++;
                        }
                    }
                }else {
                    startt.push_back(counter);
                    counterk++;
                }
            }
        }

        //startS[counterk]=counter;
        startt.push_back(counter);

        N = counter;
        file.close();
        cerr<<"K : "<< counterk<<endl;
        cerr<<"N : " << counter <<endl;
        K=counterk;


        Sorted= new uint[K+1];
        startS = new ulong[K+1];
        matched = new stack_node*[K+1];
        startS[0]=0; Sorted[0] = 0; matched[0] = NULL;
        for(int i = 1; i <= K; i++){
            Sorted[i] = i;
            startS[i] = startt[i];
            matched[i] = NULL;
            if( _m > startS[i] - startS[i-1] ) _m = startS[i] - startS[i-1];
        }
        startt.clear();
        delete [] memblock;

    }

    cerr << "min len : " << _m << endl;
    int c1 = _c1; int c2 = _c2;
    int m;
    //cerr << "N/K/c1 : " << N/K/c1 << "  N/K/c2 : " << N/K/c2 << endl;
    if( N/K/c1 <= _m && _m <= N/K/c2 ){
        m = _m;
    }
    else if( _m < N/K/c1 ){
        m = N/K/c1;
    }
    else if( N/K/c2 < _m ){
        m = N/K/c2;
    }
    double start, end;
    start = omp_get_wtime();
    omp_set_num_threads(threads);
    AuxData* ad = new AuxData();

    if( threads == 1 )
        ad->Preprocessing(text, K, startS, min, m);
    else
        ad->Preprocessing_par(text, K, startS, min, threads, m);

    struct tree_node **sptr = new tree_node*[K];
    for(int i = 0; i < K; i++) sptr[i] = NULL;

    end = omp_get_wtime();

    cerr<<"Running Time(sec) for AuxData : "<<end - start<<endl;
    struct tree_node *ptr;
    start = omp_get_wtime();
    ptr  =create_tree_modified(text,K,N,startS,matched);

    int counter1=0;
    traverse_tree_modified(ptr,Sorted,&counter1,matched,0, sptr, NULL, 0);
    if(ad->B > min)
        traverse_tree_prefix(ptr, text, startS, Sorted, 0, 0, ad);
    end = omp_get_wtime();
    cerr<<"Running Time(sec) for compact prefix tree, sorted, Bprefix : "<<end-start<<endl;
    find_all_pairs(ptr,text,N,K,startS,Sorted,output,threads,distribution_method,min,ad,sptr);

    delete [] startS;

}

bool strcompare(const char* a, const char* b, int l){
    for(int i = 0; i < l; i++){
        if( a[i] != b[i] ) return false;
    }
    return true;
}

int main(int argc,char *argv[])
{

    char *filename;
    int output;
    int method,processors,minlength,c1,c2;
    // Default values 
    c1 = 16; c2 = 8;
    method = 2;
    processors = omp_get_num_procs();
    minlength = 1;

    if (argc<2){
        return 0;
    }

    filename = argv[1];


    for(int i=2;i<argc;i++){
        if (argv[i][0]=='-'){
            if (argv[i][1]=='p')
                processors = atoi(argv[i+1]);
            else if ( strcompare(argv[i]+1, "om", 2) )
                minlength=atoi(argv[i+1]);
            else if ( strcompare(argv[i]+1, "output", 6) )
                output=atoi(argv[i+1]);
            else if ( strcompare(argv[i]+1, "constant", 8) ){
                c1=atoi(argv[i+1]);
                c2=atoi(argv[i+2]);
            }

        }
    }

    display_choices(filename,output,processors,minlength,method);

    APSP(filename,output,processors,minlength,method,c1,c2);
    return 0;
}



int choose_index(int a){
    if (a=='A') 
        return 1;
    else if (a=='C')
        return 2;
    else if (a=='G')
        return 3;
    else if (a=='T')
        return 4;
    else
        return 0;
}




bool isTerminatorfile(uchar i){

    return (i==SEPERATOR);
}



void display_choices(char *filename, int output, int threads,int min,int distribution_method){

    cerr << "== FastAPSP ================================" << endl;
    cerr << "File :" << filename << endl;
    cerr << "om option : " << min << endl;
    cerr << "output option : " << output << endl;
    cerr << "parallel option : " << threads << endl;
}


void encode(uchar *final,ulong counter,int bitnum, char c){
    if (c=='A'){
        final[counter] &= ~(1<<bitnum);
        final[counter] &= ~(1<<((bitnum+1)%8));
    }else if (c=='C'){
        final[counter] |= (1<<bitnum);
        final[counter] &= ~(1<<((bitnum+1)%8));
    }else if (c=='G'){
        final[counter] &= ~(1<<bitnum);
        final[counter] |= (1<<((bitnum+1)%8));
    }else if (c=='T'){
        final[counter] |= (1<<bitnum);
        final[counter] |= (1<<((bitnum+1)%8));
    }
}


