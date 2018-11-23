#ifndef _GA_THREAD_SYNC_
#define _GA_THREAD_SYNC_

#include <thread>
#include <vector>
#include <utility>
#include <condition_variable>
#include <tuple>
#include "bool_array.h"
#include "ga_info.h"

namespace opt
{
	template<class F> class GAGroup;

	// GA线程同步器
	template<class R, class... Args>
	class GAThreadSync
	{
		friend class GAGroup<R(Args...)>;

	private:
		std::mutex mtx;                                   // 互斥量
		std::condition_variable cv;                       // 条件变量, 用于线程同步
		bool crossReady;                                  // 交叉操作标志
		bool mutReady;                                    // 变异操作标志
		bool selectReady;                                 // 选择操作标志
		int threadNum;                                    // 并行计算使用的线程数
		bool_array cross_flag;                            // 并行计算: 交叉线程状态标志
		bool_array mut_flag;                              // 并行计算: 变异线程状态标志
		bool_array sel_flag;                              // 并行计算: 选择线程状态标志
		GAGroup<R(Args...)>* ga;                          // 同步器关联的GA种群

	public:
		// 常规构造
		GAThreadSync(GAGroup<R(Args...)>* ptr)
			:crossReady(false),
			mutReady(false),
			selectReady(false),
			threadNum(1),
			cross_flag(threadNum),
			mut_flag(threadNum),
			sel_flag(threadNum),
			ga(ptr)
		{}

		// 复制构造
		GAThreadSync(const GAThreadSync& other, GAGroup<R(Args...)>* ptr)
			:crossReady(other.crossReady),
			mutReady(other.mutReady),
			selectReady(other.selectReady),
			threadNum(other.threadNum),
			cross_flag(other.cross_flag),
			mut_flag(other.mut_flag),
			sel_flag(other.sel_flag),
			ga(ptr)
		{}

		// 移动构造
		GAThreadSync(GAThreadSync&& other) = delete;

		void setThreadNum(const int N);	            // 设置线程数量

		void cross_sync(const int thread_seq);		// "交叉"线程同步
		void mutate_sync(const int thread_seq);		// "变异"线程同步
		void select_sync(const int thread_seq);		// "选择"线程同步
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// 设置线程数量
	template<class R, class... Args>
	void GAThreadSync<R, Args...>::setThreadNum(const int N)
	{
		threadNum = N;
		cross_flag.set_length(N);
		mut_flag.set_length(N);
		sel_flag.set_length(N);
	}

	// “交叉”线程同步
	template<class R, class... Args>
	void GAThreadSync<R, Args...>::cross_sync(const int thread_seq)
	{
		cross_flag[thread_seq] = true;
		selectReady = false;

		// 若交叉线程组全部运行完
		if (cross_flag.is_all_true())
		{
			// 交换个体数组
			ga->switchIndivArray();

			cross_flag.set_all(false);
			crossReady = false;
			mutReady = true;
		}
		else
		{
			mutReady = false;
		}
	}

	// "变异"线程同步
	template<class R, class... Args>
	void GAThreadSync<R, Args...>::mutate_sync(const int thread_seq)
	{
		mut_flag[thread_seq] = true;
		crossReady = false;

		// 若变异线程组全部运行完
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

	// "选择"线程同步
	template<class R, class... Args>
	void GAThreadSync<R, Args...>::select_sync(const int thread_seq)
	{
		sel_flag[thread_seq] = true;
		mutReady = false;

		// 若选择线程组全部运行完
		if (sel_flag.is_all_true())
		{
			// 淘汰子代最差个体, 用父代最优个体取代
			ga->indivs[ga->group_state.worstIndex] = ga->bestIndivs.back();

			// 如果变异后出现更优个体
			if (ga->indivs[ga->group_state.bestIndex].fitness >= ga->bestIndivs.back().fitness)
			{
				ga->bestIndivs.push_back(ga->indivs[ga->group_state.bestIndex]);
			}
			else
			{
				ga->bestIndivs.push_back(ga->bestIndivs.back());
			}

			// 更新轮盘赌刻度线
			ga->updateRoulette(ga->groupSize);

			// 更新GA种群停止状态
			ga->updateStopState();

			ga->flushStopFlag();

			sel_flag.set_all(false);
			selectReady = false;
			crossReady = true;

			if (ga->group_state.sleep.signal == true)
			{
				ga->group_state.sleep.result = true;
			}

			// 调用外部监听函数
			if (ga->monitor)
			{
				GA_Info ga_info(ga->group_state.time, ga->group_state.nGene, ga->bestIndivs.back());
				ga->monitor(ga_info);
			}

			// 更新子代个体数量
			if (ga->resize)
			{
				std::size_t count = ga->resize(ga->group_state.nGene + 1);
				if (count < 0) { count = 0; }
				ga->groupSize = count + count % 2;
			}
		}
		else
		{
			crossReady = false;
		}
	}
}

#endif
