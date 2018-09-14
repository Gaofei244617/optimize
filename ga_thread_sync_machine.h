#ifndef _GA_THREAD_SYNC_MACHINE_
#define _GA_THREAD_SYNC_MACHINE_

#include <thread>
#include <vector>
#include <utility>
#include <condition_variable>
#include "bool_array.h"
#include <iostream>

namespace opt
{
	template<class F>
	class GAGroup;

	// GA线程同步器
	template<class R, class... Args>
	class GAThreadSyncMachine
	{
	private:
		friend class GAGroup<R(Args...)>;

		std::mutex mtx;  // 互斥量
		std::condition_variable cv; // 条件变量, 用于线程同步
		bool crossReady;  // 交叉操作标志
		bool mutReady;    // 变异操作标志
		bool selectReady; // 选择操作标志
		int threadNum;    // 并行计算使用的线程数
		bool_array cross_flag; // 并行计算: 交叉线程状态标志
		bool_array mut_flag;   // 并行计算: 变异线程状态标志
		GAGroup<R(Args...)>* ga;

	public:
		// 常规构造
		GAThreadSyncMachine(GAGroup<R(Args...)>* ptr) 
			:crossReady(false),
			mutReady(false),
			selectReady(false),
			threadNum(1),
			cross_flag(threadNum),
			mut_flag(threadNum),
			ga(ptr)
		{}

		// 复制构造
		GAThreadSyncMachine(const GAThreadSyncMachine& other, GAGroup<R(Args...)>* ptr)
			:crossReady(other.crossReady),
			mutReady(other.mutReady),
			selectReady(other.selectReady),
			threadNum(other.threadNum),
			cross_flag(other.cross_flag),
			mut_flag(other.mut_flag),
			ga(ptr)
		{}

		// “交叉”线程同步
		void cross_sync(const int thread_seq)
		{
			cross_flag[thread_seq] = true;
			selectReady = false;
			if (cross_flag.is_all_true())
			{
				cross_flag.set_all(false);
				crossReady = false;
				mutReady = true;
			}
			else
			{
				mutReady = false;
			}
		}

		// “变异”线程同步
		void mutate_sync(const int thread_seq)
		{
			mut_flag[thread_seq] = true;
			crossReady = false;
			if (mut_flag.is_all_true())
			{
				mut_flag.set_all(false);
				mutReady = false;
				selectReady = true;
			}
			else
			{
				selectReady = false;
			}
		}

		// “选择”线程同步
		void select_sync()
		{
			selectReady = false;
			crossReady = true;
			mutReady = false;
		}
	};
}

#endif
