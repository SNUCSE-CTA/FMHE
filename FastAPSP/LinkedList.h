#include <iostream>

using namespace std;

class Node{
public:
    unsigned int v;
    Node* next;
};

class LinkedList{
public:
    Node* root;
    LinkedList(){
        root = NULL;
    }
    void push_front(unsigned int q);

    unsigned int front();

    bool empty(){
        if( root == NULL ) return true;
        return false;
    }

};
