#include <iostream>
#include <random>
#include <functional>

using namespace std;

int main(){
	int n=10000000;
	default_random_engine gn;
	uniform_real_distribution<double> dst(1, 10000);
	auto ds = bind(dst, gn);
	for(int i=0;i<n;i++){
		cout<<ds()<<endl;
	}
	cout<<endl;
}
