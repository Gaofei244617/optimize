#ifndef _GA_THREADSTATE_H_
#define _GA_THREADSTATE_H_

#include <thread>
#include <vector>
#include <utility>
#include <condition_variable>
#include <atomic>
#include "bool_array.h"

namespace opt
{
	// GA算法并行计算相关状态标志
	struct GA_ThreadState
	{
		std::mutex mtx;  // 互斥量
		std::condition_variable cv; // 条件变量, 用于线程同步
		bool crossReady;  // 交叉操作标志
		bool mutReady;    // 变异操作标志
		bool selectReady; // 选择操作标志
		int threadNum;    // 并行计算使用的线程数
		bool_array cross_flag; // 并行计算: 交叉线程状态标志
		bool_array mut_flag;   // 并行计算: 变异线程状态标志

		// 构造函数
		GA_ThreadState()
			:crossReady(false),
			mutReady(false),
			selectReady(false),
			threadNum(1),
			cross_flag(threadNum),
			mut_flag(threadNum)
		{}

		// 复制构造
		GA_ThreadState(const GA_ThreadState& other)
			:crossReady(other.crossReady),
			mutReady(other.mutReady),
			selectReady(other.selectReady),
			threadNum(other.threadNum),
			cross_flag(other.cross_flag),
			mut_flag(other.mut_flag)
		{}
	};
}

#endif
