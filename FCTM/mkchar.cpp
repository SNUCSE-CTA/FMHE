#include <iostream>
#include <random>
#include <functional>

using namespace std;

int main(){
	int n=10000000;
	default_random_engine gn;
	uniform_int_distribution<int> dst(-128, 127);
	auto ds = bind(dst, gn);
	for(int i=0;i<n;i++){
		cout<<ds()<<' ';
	}
	cout<<endl;
}
