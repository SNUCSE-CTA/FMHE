#include <iostream>
#include "Ctree.h"
#include <omp.h>
#include <list>
using namespace std;


struct tree_node *create_tree_modified(uchar *text,uint k,ulong n, ulong startp[],struct stack_node **stacks){
    struct tree_node *root=create_node(0,0,0);
    for(uint i=0;i<k;i++){
        struct tree_node *curptr=root;
        struct tree_node *parent=root;
        ulong pos=1;long currpos=0;int choosen_child=1;
        ulong c= startp[i];
        //cout<<"k:"<<i<<endl;
        while (1){
            if (curptr->pos>=pos) { //still in the same node
                if (decode(text,c) /*text[c]*/==  decode(text,startp[curptr->from]+currpos) && c!=startp[i+1]/* added here recently*/   /*text[startp[curptr->from]+currpos]*/) { 
                    //putchar(text[c]); 
                    //cout<<"in a node : startp[sorted[curptr->from]]+currpos"<<startp[curptr->from]+currpos<<endl;
                    pos++;currpos++;c++;

                } else {  // Split 
                    //cout<<"enter splitting"<<endl;
                    struct tree_node *temp = create_node(pos-1,curptr->from,curptr->to); //new mother node
                    if (c==startp[i+1])
                        temp->ptr[0] = create_node(find_pos_modified(startp,i,k,n)-currpos,i,i);  //new node for the seperator if end of string is reached
                    else
                        temp->ptr[get_index(decode(text,c) /*text[c]*/)]= create_node(find_pos_modified(startp,i,k,n)-currpos,i,i);  //new node 
                    parent->ptr[choosen_child]=temp;
                    //cout<<"choosen_child:"<<choosen_child<<endl;
                    //cout<<"startp[sorted[curptr->from]]+currpos:"<<startp[sorted[curptr->from]]+currpos<<endl;
                    temp->ptr[get_index(/*text[startp[curptr->from]+currpos]*/ decode(text,startp[curptr->from]+currpos))]=curptr;
                    curptr->pos= curptr->pos-pos;
                    //cout<<"curpos->"<<curptr->pos<<endl;
                    remove_index(curptr,i);
                    break;
                }
            }else if (curptr->pos<pos){  // moving to another node either old or new node
                if (curptr->ptr[get_index(decode(text,c)/*text[c]*/)]!=NULL &&  c!=startp[i+1]){
                    choosen_child=get_index(decode(text,c)/*text[c]*/);
                    parent=curptr;
                    curptr= curptr->ptr[choosen_child];
                    //cout << "move toward a child node"<<endl;
                    if (c==startp[i+1]   /*text[c]==SEPERATOR*/){  /* useless now */
                        stacks[curptr->from]= create_stack_node(i,stacks[curptr->from]);
                        //cout <<"matched 2 " <<curptr->from << " "<<i<<endl;
                        break;
                    }
                    pos=1; currpos++;c++;

                }else {
                    if (/*text[c]==SEPERATOR*/ c==startp[i+1] && curptr->ptr[1]==NULL && curptr->ptr[2]==NULL && curptr->ptr[3]==NULL && curptr->ptr[4]==NULL){
                        stacks[curptr->from]= create_stack_node(i,stacks[curptr->from]);
                        //cout <<"create stack node" <<curptr->from << " "<<i<<endl;
                        break;
                    }else if (c!=startp[i+1]  && curptr->ptr[1]==NULL && curptr->ptr[2]==NULL && curptr->ptr[3]==NULL && curptr->ptr[4]==NULL && curptr!=root){
                        curptr->ptr[get_index(decode(text,c) /*text[c]*/)]= create_node(find_pos_modified(startp,i,k,n)-currpos,i,i);
                        curptr->ptr[0]= create_node(0,curptr->from,curptr->to);

                        //cout<<"create 2 new nodes"<<endl;
                        break;
                    }else {
                        if (c!=startp[i+1]   /*text[c]==SEPERATOR*/){
                            curptr->ptr[get_index(decode(text,c) /*text[c]*/)]= create_node(find_pos_modified(startp,i,k,n)-currpos,i,i);
                            //cout<<"create new nodes - not terminator"<<endl;
                        }else {
                            if (curptr->ptr[0]==NULL){
                                curptr->ptr[0]= create_node(0,i,i);
                                //cout<<"create new nodes - terminator "<<endl;
                            }else {
                                stacks[curptr->ptr[0]->from]= create_stack_node(i,stacks[curptr->ptr[0]->from]);
                                //cout <<"create stack node 2" <<curptr->from << " "<<i<<endl;
                            }
                        }
                        //display_node(curptr->ptr[get_index(text[c])]);

                        break;
                    }
                }
            }
        }

    }
    return root;
}

struct tree_node *create_node(int pos,int from,int to){
    struct tree_node*ptr = new struct tree_node();
    if (pos<0) pos=0;
    ptr->pos=pos;
    ptr->from=from;
    ptr->to=to;
    for(int i=0;i<5;i++)
        ptr->ptr[i]=NULL;
    return ptr;
}

void display_node(struct tree_node* p){
    //cout<<"Pos:"<<p->pos<<endl;
    cout<<"From:"<<p->from<<endl;
    cout<<"To:"<<p->to<<endl;
    for(int i = 1; i < 5; i++){
        cout << p->rem[i] << " ";  
    }
    cout << endl;
}

void display_tree(struct tree_node* p){

    display_node(p);
    for(int i=0;i<5;i++){
        if (p->ptr[i]!=NULL){
            cout<<"Branch ";
            putchar(get_char(i));
            cout<<endl;
            display_tree(p->ptr[i]);
        }
    }

}


int get_index(uchar c){
    if (c==SEPERATOR)
        return 0;
    else if (c=='A')
        return 1;
    else if (c=='C')
        return 2;
    else if (c=='G')
        return 3;
    else if (c=='T')
        return 4;
    else
        return -1;
}

uchar get_char(int c){
    if (c==0)
        return SEPERATOR;
    else if (c==1)
        return 'A';
    else if (c==2)
        return 'C';
    else if (c==3)
        return 'G';
    else if (c==4)
        return 'T';
    else
        return SEPERATOR;
}

ulong find_pos(ulong startpos[],uint sorted[],int i,ulong k,ulong n){
    if (sorted[i]<k-1)
        return startpos[sorted[i]+1]-startpos[sorted[i]]-1;
    else
        return n-startpos[sorted[i]]-1;
}

ulong find_pos_modified(ulong startpos[],int i,ulong k,ulong n){
    if ((uint)i<k-1)
        return startpos[i+1]-startpos[i]-1;
    else
        return n-startpos[i]-1;
}


void updateNodeRange(struct tree_node *ptr,int i){
    if (ptr->from > i ) ptr->from =i;
    if (ptr->to < i) ptr->to=i;
}

void remove_index(struct tree_node *ptr,int i){
    if (ptr->to==i) ptr->to=i-1;
}


void find_all_pairs(struct tree_node *ptr,uchar *T,ulong N,ulong k,ulong startpos[],uint sorted[],int output,int threads,int distribution_method,int min, AuxData* ad, struct tree_node** sptr){
    int **A;

    if (output==1){
        A= new int*[k];
        for (uint z=0;z<k;z++)
            A[z]= new int[k];

        for (uint z=0;z<k;z++){
            for (uint z1=0;z1<k;z1++)
                A[z][z1]=0;
        }
    }
    double starttime,endtime;
    starttime = omp_get_wtime();
    if (threads==1){
        matching_step(T,startpos,sorted,A,ptr,k,N,output,min, ad, sptr);
    }
    else{
        oList* ol = new oList[k];
        matching_step_par(threads, T, startpos, sorted, A, ptr, k, N, output, min, ad, ol, sptr);
        if( output == 2 || output == 3 ){
            for(int i = 0; i < k; i++){
                ol[i].output_step(i, sorted, output);
            }
        }
    }
    endtime = omp_get_wtime();

    cerr << "Running Time(sec) for matching and output step : " << endtime - starttime << endl;

    /*starttime = omp_get_wtime();
    for(int i = 0; i < k; i++){
        rl[i].Print_RangeList(i, sorted, output);
    }
    endtime = omp_get_wtime();

    cerr << "User Time for output step: " << endtime - starttime << endl;*/

    if (output==1){
        for (uint z=0;z<k;z++){
            for (uint z1=0;z1<k;z1++){
                if (z==z1)
                    printf("%5d ",0);
                else
                    printf("%5d ",A[z][z1]);
            }
            cout <<endl;
        }
    }

}


void print_depth(int d){
    for(int i = 0; i < d; i++){
        cout << "   ";
    }
}


void traverse_tree_prefix(struct tree_node *ptr, uchar *T, ulong startpos[], uint sorted[], int curpos, int hash, AuxData* ad){
    int rem = ad->B-1 - curpos;
    if( rem == 0 ) return;
    if( rem >= ptr->pos ){
        for(int i = 0; i < ptr->pos; i++){
            hash <<= ad->ShiftBits;
            hash += ad->rank[decode(T, startpos[sorted[ptr->from]] + curpos + i)];
            int Pp = ad->Bprefixp[curpos+i];
            ad->Bprefix[Pp+hash] = ptr;
        }
        if( rem == ptr->pos ) return;
        curpos += ptr->pos;
        for(int i = 1; i < 5; i++){
            if( ptr->ptr[i] != NULL ){
                int hash2 = hash;
                hash2 <<= ad->ShiftBits;
                hash2 += ad->rank[decode(T, startpos[sorted[ptr->ptr[i]->from]] + curpos)];
                int Pp = ad->Bprefixp[curpos];
                ad->Bprefix[Pp+hash2] = ptr->ptr[i];
                traverse_tree_prefix(ptr->ptr[i], T, startpos, sorted, curpos+1, hash2, ad);
            }
        }
    }
    else{
        for(int i = 0; i < rem; i++){
            hash <<= ad->ShiftBits;
            hash += ad->rank[decode(T, startpos[sorted[ptr->from]] + curpos + i)];
            int Pp = ad->Bprefixp[curpos+i];
            ad->Bprefix[Pp + hash] = ptr;
        }
        return;
    }
}

int traverse_tree_modified(struct tree_node *ptr,uint sorted [],int *counter,struct stack_node **stacks,int depth, struct tree_node **sptr, struct tree_node* parent, int sdepth){
    bool first=false;
    int ret = 0;
    ptr->sdepth = sdepth;
    for(int i=0;i<5;i++){
        if (ptr->ptr[i]!=NULL){
            if( i == 0 ){
                ptr->rem[i] = traverse_tree_modified(ptr->ptr[i], sorted, counter, stacks, depth+1, sptr, ptr, sdepth+ptr->pos+1) + ptr->ptr[i]->pos + 1;
            }
            else{
                ptr->rem[i] = traverse_tree_modified(ptr->ptr[i], sorted, counter, stacks, depth+1, sptr, NULL, sdepth+ptr->pos+1) + ptr->ptr[i]->pos + 1;
            }
            if( i == 0 ){ 
                ptr->rem[i]--;
            }
            if( ret < ptr->rem[i] ) ret = ptr->rem[i];
            if (!first){			
                ptr->from=ptr->ptr[i]->from;
            }
            else
                ptr->to=ptr->ptr[i]->to;
            first=true;
        }
    }

    if (!first){   
        sorted[*counter]=ptr->from;
        struct stack_node*temp= stacks[ptr->from]; 
        if( parent != NULL ) sptr[ptr->from] = parent;
        else sptr[ptr->from] = ptr;
        ptr->from=*counter;
        *counter=*counter+1;
        while (temp!=NULL){
            sorted[*counter]=temp->value;
            if( parent != NULL ) sptr[temp->value] = parent;
            else sptr[temp->value] = ptr;
            *counter=*counter+1;
            temp=temp->next;
        }

        ptr->to=(*counter)-1;

    }
    return ret;
}

struct stack_node* create_stack_node(uint i,struct stack_node *ptr){
    struct stack_node *temp=new struct stack_node();

    temp->value=i;
    if (ptr!=NULL)
        temp->next=ptr;
    else
        temp->next=NULL;

    return temp;
}


void traversal1(struct tree_node *ptr,int **A,struct tree_node *root,uchar *T,int curpos,int curpos2,ulong startpos[],uint sorted[],ulong N,ulong k,int output){

    //cout<<"Outer from:"<<ptr->from<<" to:"<<ptr->to<< endl;
    for(int i=0;i<=ptr->pos;i++){
        traversal2(root->ptr[get_index(T[startpos[sorted[ptr->from]]+curpos])],0,ptr,i,T,curpos,curpos2,k,N,startpos,sorted,A,output);
        curpos++;
    }

    for(int n=1;n<=4;n++){
        if (ptr->ptr[n]!=NULL) traversal1(ptr->ptr[n],A,root,T,curpos,0,startpos,sorted,N,k,output);
    }
}

void traversal2(struct tree_node* ptr2,int i2,struct tree_node *ptr1,int i1,uchar * T,int curpos,int curpos2,uint k,int N,ulong startpos[],uint sorted[],int **A,int output){

    if (ptr2==NULL) return;
    for (int n=i1;n<=ptr1->pos;n++){
        //cout <<"inner helper :ptr1 "<<ptr1->from <<"-"<<ptr1->to<<" ptr2:"<<ptr2->from <<"-"<<ptr2->to<<" i1="<<n<<" i2="<<i2<< " curpos1:"<<curpos<<" curpos2:"<<curpos2<<endl;	
        if (i2<=ptr2->pos){
            if (T[startpos[sorted[ptr1->from]]+curpos] != T[startpos[sorted[ptr2->from]]+curpos2])
                return;
        }else{
            ptr2=ptr2->ptr[get_index(T[startpos[sorted[ptr1->from]]+curpos])];			
            i2=0;
            if (ptr2==NULL) return;
        }
        curpos++;curpos2++;i2++;
    }


    for(int n=1;n<=4;n++){
        if (ptr1->ptr[n]!=NULL){
            traversal2(ptr2,i2,ptr1->ptr[n],0,T,curpos,curpos2,k,N,startpos,sorted,A,output);
        }
    }


    if (ptr1->ptr[1]==NULL && ptr1->ptr[2]==NULL && ptr1->ptr[3]==NULL && ptr1->ptr[4]==NULL){  //report a match.
        //cout<< "Report a match!"<<endl;
        for(int z=ptr1->from;z<=ptr1->to;z++){
            //cout << "curpos2:" <<curpos2<<endl;
            if (output==1){
                for(int z1=ptr2->from;z1<=ptr2->to;z1++){ 
                    if (A[sorted[z]][sorted[z1]]==0 && z1!=z){
                        //cout << "set value "<<z1<<" "<<z<<endl; 
                        A[sorted[z]][sorted[z1]]= curpos2;
                    }
                }
            }
        }
        return;
    }	

}


void matching_step_par(int threads, uchar *T, ulong startpos[], uint sorted[], int **A, struct  tree_node* ptr, ulong k, ulong N, int output, int min, AuxData* ad, oList* ol, struct tree_node** sptr)
{

#pragma omp parallel for schedule(dynamic) 
    for (int i = 0; i < k; i++){
        ulong j = startpos[i] + 1;
        ulong next = startpos[i + 1];
        ulong len = startpos[i+1] - startpos[i];
        //st->init();
        if( len >= min ){
            insert_output(i, sorted, len, sptr[i], output, &ol[i], A);
        }
        else continue;
        ulong ix = j + ad->m - 1;
        // remaining string length >= m
        // cerr << "Case 1 : " << endl;
        while (ix < next){
            if( ix-ad->m+1 + min > next ) break;
            unsigned int hash = 0;
            for (int r = ad->B - 1; r >= 0; r--){
                hash <<= ad->ShiftBits;
                hash += ad->rank[decode(T, ix - r)];
            }

            unsigned int q = ad->qrm[hash];
            if (q == ad->m){
                ulong v = ix - ad->m + 1; struct tree_node *curptr = ptr; int pos = 1; int curpos = 0; int sw = 0;
                while (1){
                    if (sw) break;
                    if (v == next && next - (ix - ad->m + 1) >= min){  // There is a match
                        insert_output(i, sorted, next - (ix - ad->m + 1), curptr, output, &ol[i], A);
                        break;
                    }

                    if (pos <= curptr->pos){
                        if (decode(T, v) == decode(T, startpos[sorted[curptr->from]] + curpos)){
                            pos++; curpos++; v++;
                        }
                        else
                            break;
                    }
                    else {
                        if (curptr->ptr[get_index(decode(T, v))] != NULL){
                            if (curptr->rem[get_index(decode(T, v))] < next - v){
                                sw = 1;
                            }
                            curptr = curptr->ptr[get_index(decode(T, v))];
                            pos = 1; curpos++; v++;
                        }
                        else
                            break;
                    }
                }
                ix ++;
            }
            ix += ad->m-q;
        }

        // B < remaining string length < m
        //cerr << "Case 2 : " << endl;
        unsigned int hash2 = 0;
        if( len >= ad->B && min < ad->m ){
            ix = next - ad->m + 1;
            for (int r = next - ad->B; r < next; r++){
                hash2 <<= ad->ShiftBits;
                hash2 += ad->rank[decode(T, r)];
            }

            Node* it = ad->qList[hash2].root;
            while (it != NULL){
                if( it->v + ad->B > len ){ it = it->next; continue; }
                ix = next - it->v;
                if( next-ix < min ){ it = it->next; continue;}
                ulong v = ix; struct tree_node *curptr = ptr; int pos = 1; int curpos = 0; int sw = 0;
                while (1){
                    if (sw) break;
                    if (v == next && next - ix >= min){  // There is a match
                        if( len != next-ix)
                            insert_output(i, sorted, next - ix, curptr, output, &ol[i], A);
                        break;
                    }
                    if (pos <= curptr->pos){
                        if (decode(T, v) == decode(T, startpos[sorted[curptr->from]] + curpos)){
                            pos++; curpos++; v++;
                        }
                        else
                            break;
                    }
                    else {
                        if (curptr->ptr[get_index(decode(T, v))] != NULL){
                            if (curptr->rem[get_index(decode(T, v))] < next - v){
                                sw = 1;
                            }
                            curptr = curptr->ptr[get_index(decode(T, v))];
                            pos = 1; curpos++; v++;
                        }
                        else
                            break;
                    }
                }
                it = it->next;
            }
        }
        // remaining string length < B
        //cerr << "Case 3 : " << endl;
        hash2 = 0;
        int tlen = len >= ad->B-1 ? ad->B-1 : len;
        ix = next - tlen;
        int BMask = 0;
        for(int r = ix; r < next; r++){
            hash2 <<= ad->ShiftBits;
            hash2 += ad->rank[decode(T, r)];
            BMask <<= ad->ShiftBits;
            BMask += (1<<ad->ShiftBits)-1;
        }
        for (int r = 0; r < tlen; r++){
            if( tlen-r < min ) break;
            int Pp = ad->Bprefixp[tlen-r-1];
            if (ad->Bprefix[Pp + hash2] != NULL && tlen - r >= min){
                if( len != tlen-r )
                    insert_output(i, sorted, tlen - r, ad->Bprefix[Pp + hash2], output, &ol[i], A);
            }
            BMask >>= ad->ShiftBits;
            hash2 &= BMask;
        }
    }
}

double case1_time, case2_time, case3_time;

void matching_step(uchar *T, ulong startpos[], uint sorted[], int **A, struct  tree_node* ptr, ulong k, ulong N, int output, int min, AuxData* ad, struct tree_node** sptr){
    int findcall = 0;
    double aver_find_len = 0.0;
    case1_time = case2_time = case3_time = 0.0;
    //double aver_shift = 0.0;
    for(int i = 0; i < k; i++){
        oList* range = new oList();
        ulong j = startpos[i] + 1;
        ulong next = startpos[i + 1];
        ulong len = startpos[i+1]-startpos[i];
        double start, end;
        start = omp_get_wtime();
        if( len >= min ){
            insert_output(i, sorted, len, sptr[i], output, range, A);
        }
        else continue;
        // ix means p+m-1 in our algorithm
        //cerr << "Case 1" << endl;
        ulong ix = j + ad->m - 1;
        while (ix < next){
            //cerr << "ix : " << ix << endl;
            if( ix-ad->m+1 + min > next ) break;
            unsigned int hash = 0;
            for (int r = ad->B - 1; r >= 0; r--){
                hash <<= ad->ShiftBits;
                hash += ad->rank[decode(T, ix - r)];
            }
            unsigned int q = ad->qrm[hash];
            if (q == ad->m){
                /*unsigned int hasht = 0; 
                for(ulong r = ix-ad->m+1; r < ix-ad->m+1+ad->B-1; r++){
                    hasht <<= ad->ShiftBits;
                    hasht += ad->rank[decode(T, r)];
                }*/
                //if( ad->Bprefix[ ad->Bprefixp[ad->B-2] + hasht] != NULL) {
                
                findcall++;       
                ulong v = ix - ad->m + 1; struct tree_node *curptr = ptr; int pos = 1; int curpos = 0; int sw = 0;
                //ulong v = ix - ad->m + 1 + ad->B-1; struct tree_node *curptr = ad->Bprefix[ ad->Bprefixp[ad->B-2]+hasht]; int pos = ad->B-1 - curptr->sdepth + 1; int curpos = ad->B-1; int sw = 0;
                while(1){
                    if (sw) break;
                    if (v == next && next-(ix-ad->m+1) >= min){  // There is a match
                        insert_output(i, sorted, next - (ix-ad->m+1), curptr, output, range, A);
                        break;
                    }

                    if (pos <= curptr->pos){
                        if (decode(T, v) == decode(T, startpos[sorted[curptr->from]] + curpos)){
                            pos++; curpos++; v++;
                        }
                        else
                            break;
                    }
                    else {
                        if (curptr->ptr[get_index(decode(T, v))] != NULL){
                            if (curptr->rem[get_index(decode(T, v))] < next - v){
                                sw = 1;
                            }
                            curptr = curptr->ptr[get_index(decode(T, v))];
                            pos = 1; curpos++; v++;
                        }
                        else
                            break;
                    }
                }
                aver_find_len += v - (ix-ad->m+1);
                //}
                ix++; //aver_shift += 1;
            }
            ix += ad->m-q; //aver_shift += (ad->m-q);
        }
        end = omp_get_wtime();
        case1_time += (end - start);
        //cerr << "Case 2" << endl;
        unsigned int hash2 = 0;
        if( len >= ad->B && min < ad->m ){
            // B < remaining string length < m
            ix = next - ad->m + 1;
            for (int r = next-ad->B; r < next; r++){
                hash2 <<= ad->ShiftBits;
                hash2 += ad->rank[decode(T, r)];
            }
            //cerr << "Case 2 : hash" << endl;
            Node* it = ad->qList[hash2].root;
            /*while( it != NULL ){
                cerr << it->v << endl;
                it = it->next;
            }
            it = ad->qList[hash2].root;*/
            start = omp_get_wtime();
            while( it != NULL ){
                if( it->v + ad->B > len ){it = it->next; continue; }
                ix = next - it->v;
                //cerr << "Case 2 : ix " << ix << endl;
                if( next-ix < min ){ it = it->next; continue; }
                findcall++;
                ulong v = ix; struct tree_node *curptr = ptr; int pos = 1; int curpos = 0; int sw = 0;
                while (1){
                    if (sw) break;
                    if (v == next && next-ix >= min){  // There is a match
                        if( len != next-ix ){
                            insert_output(i, sorted, next-ix, curptr, output, range, A);
                        }
                        break;
                    }
                    if (pos <= curptr->pos){
                        if (decode(T, v) == decode(T, startpos[sorted[curptr->from]] + curpos)){
                            pos++; curpos++; v++;
                        }
                        else
                            break;
                    }
                    else {
                        if (curptr->ptr[get_index(decode(T, v))] != NULL){
                            if (curptr->rem[get_index(decode(T, v))] < next - v){
                                sw = 1;
                            }
                            curptr = curptr->ptr[get_index(decode(T, v))];
                            pos = 1; curpos++; v++;
                        }
                        else
                            break;
                    }   
                } 
                aver_find_len += v - ix;
                it = it->next;
            }
            end = omp_get_wtime();
            case2_time += (end - start);
        }
        start = omp_get_wtime();
        // remaining string length < B
        // cerr << "Case 3" << endl;
        hash2 = 0;
        int tlen = len >= ad->B-1 ? ad->B-1 : len;
        ix = next - tlen;
        int BMask = 0;
        for(int r = ix; r < next; r++){
            hash2 <<= ad->ShiftBits;
            hash2 += ad->rank[decode(T, r)];
            BMask <<= ad->ShiftBits;
            BMask += (1<<ad->ShiftBits)-1;
        }
        for(int r = 0; r < tlen; r++){
            if( tlen-r < min ) break;
            int Pp = ad->Bprefixp[tlen-r-1];
            if( ad->Bprefix[Pp + hash2] != NULL && tlen-r >= min){
                if( len != tlen-r ){
                    insert_output(i, sorted, tlen-r, ad->Bprefix[Pp+hash2], output, range, A);
                }
            }
            BMask >>= ad->ShiftBits;
            hash2 &= BMask;
        }
        end = omp_get_wtime();
        case3_time += (end - start);
        range->output_step2(i, sorted, output);
        delete range;
    }
    //cerr << "# of findcall : " << findcall << endl;
    //cerr << "Aver_shift : " << aver_shift / k << endl;
    //cerr << "Aver_find_len : " << aver_find_len / findcall << endl;
    //cerr << "Case 1 : " << case1_time << endl;
    //cerr << "Case 2 : " << case2_time << endl;
    //cerr << "Case 3 : " << case3_time << endl;
}

void insert_output(int i, unsigned int* sorted, int value, struct tree_node* curptr, int output, oList* ol, int **A){
    if( output == 1 ){
        for(int r = curptr->from; r <= curptr->to; r++){
            if( A[i][sorted[r]] < value )
                A[i][sorted[r]] = value;
        }
    }
    if( output == 2 || output == 3 /*|| output == 4*/ ){
        ol->v.push_back(Range(curptr->from, curptr->to, value));
    }
    /*if( output == 5 ){
        return;
    }*/
}

char decode(uchar *final,ulong pos){
    ulong bbyte = pos/4;
    ulong bbit = (pos%4)*2+1;
    uchar c = final[bbyte];
    if (!(final[bbyte] & (1 << bbit)) && !(final[bbyte] & (1 << ((bbit+1)%8)))) return 'A' ;
    if ((final[bbyte] & (1 << bbit)) && !(final[bbyte] & (1 << ((bbit+1)%8)))) return 'C' ;
    if (!(final[bbyte] & (1 << bbit)) && (final[bbyte] & (1 << ((bbit+1)%8)))) return 'G' ;
    if ((final[bbyte] & (1 << bbit)) && (final[bbyte] & (1 << ((bbit+1)%8)))) return 'T' ;
    return 'y';
}


