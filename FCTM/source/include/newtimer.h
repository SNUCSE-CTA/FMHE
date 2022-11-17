#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <chrono>

using namespace std;

class TIMER{
public:
	std::chrono::time_point<std::chrono::steady_clock> start;
	std::chrono::time_point<std::chrono::steady_clock> end;
//	std::chrono::time_point<std::chrono::high_resolution_clock> start;
//	std::chrono::time_point<std::chrono::high_resolution_clock> end;
	bool has_start = false;
};

static inline void timer_start(TIMER *t){
	t->has_start = true;
	t->start = std::chrono::steady_clock::now();
//	t->start = std::chrono::high_resolution_clock::now();
}

static inline void timer_stop(TIMER *t){
	t->end = std::chrono::steady_clock::now();
//	t->end = std::chrono::high_resolution_clock::now();
	if( t->has_start == false)
	{
		cerr << "WARNING: Timer.End() called without Timer.Start()" << endl;
		return;
	}
	t->has_start = false; //this line must be below WARNING check
	return;
}

static inline double timer_elapsed(TIMER *t)
{
	std::chrono::duration<double> elapsed = t->end - t->start;
	return elapsed.count();
}

#endif
