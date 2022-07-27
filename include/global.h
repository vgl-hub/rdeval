#ifndef GLOBAL_H
#define GLOBAL_H

#include <mutex>
#include <chrono>
#include <queue>
#include <thread>
#include <functional>

#include "log.h"
#include "threadpool.h"

//global time
extern std::chrono::high_resolution_clock::time_point start;

extern int hc_flag;
extern int hc_cutoff;
extern short int tabular_flag;
extern int outBubbles_flag;
extern int stats_flag;
extern int verbose_flag;
extern int maxThreads;

extern std::mutex mtx;
extern ThreadPool<std::function<bool()>> threadPool;
extern Log lg;

#endif /* GLOBAL_H */

