#ifndef _OLIST_H_
#define _OLIST_H_

#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include <cstdio>
#include <set>
#include <stack>

using namespace std;

class Range{
public:
    int a, b, c;        // [a, b] -> c, c in range[a, b]

    Range(int _a = -1, int _b = -1, int _c = -1){
        a = _a; b = _b; c = _c;
    }

};

class PosInfo{
public:
    int a, b, c;           // c at a position (b is an address before sorting)
    PosInfo(int _a = -1, int _b = -1, int _c = -1){
        a = _a; b = _b; c = _c;
    }
    bool operator < (const PosInfo& p) const {
        if( a == p.a ) return c > p.c;
        return a < p.a;
    }
};

struct Compare{
    bool operator() (const PosInfo& a, const PosInfo& b) const {
        if( a.a == b.a ) return a.c < b.c;
        return a.a < b.a;
    }
};

class oList{
public:
    vector<Range> v;
    int size; 
    vector<PosInfo> s, e;

    oList(){
    }
    ~oList(){
    }
    void output_step(int stridx, unsigned int*, int output);
    void output_step2(int, unsigned int*, int);
};

#endif
