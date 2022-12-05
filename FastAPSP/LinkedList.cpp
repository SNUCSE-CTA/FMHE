#include "LinkedList.h"


void LinkedList::push_front(unsigned int q){
    Node* t = new Node;
    t->v = q;
    Node* temp = root;
    root = t;
    t->next = temp;
}

unsigned int LinkedList::front(){
    if( root != NULL ){
        return root->v;
    }
    else return -1;
}

