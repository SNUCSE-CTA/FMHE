#include "oList.h"


void oList::output_step(int stridx, unsigned int* sorted, int output){
    int size = v.size();
    if( size == 0 ) return;
    /*if( output == 4 ){
        int cnt = 0;
        for(int i = 0; i < size; i++){
            for(int j = v[i].a; j <= v[i].b; j++){
                if( stridx != sorted[j] ) cnt++;
            }
        }
        printf("%d\n", cnt);
        return;   
    }*/
    if( output == 3 ){
        for(int i = 0; i < size; i++){
            for(int j = v[i].a; j <= v[i].b; j++){
                if( stridx != sorted[j] )
                    printf("%d %d %d\n", stridx, sorted[j], v[i].c);
            }
        }
        return;
    }

    for(int i = 0; i < size; i++){
        s.push_back(PosInfo(v[i].a, i, v[i].c));
        e.push_back(PosInfo(v[i].b, i, v[i].c));
    }
    s.push_back(PosInfo(1000000000, -1, -1));
    sort(s.begin(), s.end());
    sort(e.begin(), e.end());

    map<int, int> state;
    int now = s[0].a;
    state[v[s[0].b].c] = s[0].b;
    int sp, ep; sp = 1; ep = 0;
    int top = s[0].b;

    while( ep != size ){
        if( s[sp].a > e[ep].a ){
            for(int i = now; i <= e[ep].a; i++){
                if( stridx != sorted[i] )
                    printf("%d %d %d\n", stridx, sorted[i], v[top].c);
            }
            if( now < e[ep].a + 1) now = e[ep].a + 1;
            map<int, int>::iterator it = state.find(v[e[ep].b].c);
            state.erase(it);
            top = -1;
            if(!state.empty()){
                map<int, int>::reverse_iterator rit = state.rbegin();
                top = rit->second;
            }
            ep++;
        }
        else if(s[sp].a <= e[ep].a){
            if( top != -1 ){
                for(int i = now; i < s[sp].a; i++){
                    if( stridx != sorted[i] )
                        printf("%d %d %d\n", stridx, sorted[i], v[top].c);
                }
                if( now < s[sp].a ) now = s[sp].a;
                if( now == s[sp].a && v[top].c > v[s[sp].b].c){
                    if( stridx != sorted[now] )
                        printf("%d %d %d\n", stridx, sorted[now], v[top].c);
                    now++;
                }
            }
            else now = s[sp].a;

            state[v[s[sp].b].c] = s[sp].b;
            top = -1;
            if( !state.empty() ){
                map<int, int>::reverse_iterator rit = state.rbegin();
                top = rit->second;
            }
            sp++;
        }
    }
}


void oList::output_step2(int stridx, unsigned int* sorted, int output){
    int size = v.size();
    if( size == 0 ) return;
    /*if( output == 4 ){
        int cnt = 0;
        for(int i = 0; i < size; i++){
            for(int j = v[i].a; j <= v[i].b; j++){
                if( stridx != sorted[j] ) cnt++;
            }
        }   
        printf("%d\n", cnt);
        return;   
    }*/
    if( output == 3 ){
        for(int i = 0; i < size; i++){
            for(int j = v[i].a; j <= v[i].b; j++){
                if( stridx != sorted[j] )
                    printf("%d %d %d\n", stridx, sorted[j], v[i].c);
            }
        }
        return;
    }
    
    for(int i = 0; i < size; i++){
        s.push_back(PosInfo(v[i].a, i, v[i].c));
    }
    sort(s.begin(), s.end(), Compare() );
    
    //cerr << "size : " << size << endl;

    int now = s[0].a;
    stack<PosInfo> st; st.push(s[0]);
    for(int i = 1; i < size; i++){
        //cerr << "i : " << i << endl;
        if( !st.empty() ){
            PosInfo top = st.top();
            if( v[s[i].b].b <= v[top.b].b ){
                for(int j = now; j < v[s[i].b].a; j++){
                    if( stridx != sorted[j] )
                        printf("%d %d %d\n", stridx, sorted[j], v[top.b].c);
                }
                now = v[s[i].b].a;
                st.push(s[i]);
            }
            else{
                while(!st.empty()){
                    PosInfo top = st.top();
                    if( v[s[i].b].b <= v[top.b].b ){
                        for(int j = now; j < v[s[i].b].a; j++){
                            if( stridx != sorted[j] )
                                printf("%d %d %d\n", stridx, sorted[j], v[top.b].c);
                        }
                        break;
                    }
                    for(int j = now; j <= v[top.b].b; j++){
                        if( stridx != sorted[j] )
                            printf("%d %d %d\n", stridx, sorted[j], v[top.b].c);
                    }
                    now = v[top.b].b+1;
                    st.pop();
                }
                now = v[s[i].b].a;
                st.push(s[i]);
            }
        }
        else{
            st.push(s[i]);
            now = v[s[i].b].a;
        }
    }
    while(!st.empty()){
        PosInfo top = st.top();
        //cerr << "top : " << top.a << " " << top.b << " " << top.c << endl;
        for(int j = now; j <= v[top.b].b; j++){
            if( stridx != sorted[j] )
                printf("%d %d %d\n", stridx, sorted[j], v[top.b].c);
        }
        now = v[top.b].b+1;
        st.pop();
    }
}
