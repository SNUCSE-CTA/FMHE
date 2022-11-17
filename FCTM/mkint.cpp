#include <iostream>
#include <random>
#include <functional>
#include <climits>

using namespace std;

int main(){
	int n=10000000;
	default_random_engine gn;
	uniform_int_distribution<int> dst(INT_MIN, INT_MAX);
	auto ds = bind(dst, gn);
	for(int i=0;i<n;i++){
		cout<<ds()<<' ';
	}
	cout<<endl;
}
