#ifndef _GA_THREAD_SYNC_
#define _GA_THREAD_SYNC_

#include <thread>
#include <vector>
#include <utility>
#include <condition_variable>
#include <tuple>
#include "bool_array.h"

namespace opt
{
	template<class F> class GAGroup;

	// GA线程同步器
	template<class R, class... Args>
	class GAThreadSync
	{
	private:
		friend class GAGroup<R(Args...)>;

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

		// 设置线程数量
		void setThreadNum(const int N);

		// “交叉”线程同步
		void cross_sync(const int thread_seq);
		// “变异”线程同步
		void mutate_sync(const int thread_seq);
		// “选择”线程同步
		void select_sync(const int thread_seq);
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
			// Do more things..........

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

	// “选择”线程同步
	template<class R, class... Args>
	void GAThreadSync<R, Args...>::select_sync(const int thread_seq)
	{
		sel_flag[thread_seq] = true;
		mutReady = false;
		
		// 若选择线程组全部运行完
		if (sel_flag.is_all_true())
		{
			// 淘汰子代最差个体, 用父代最优个体取代
			ga->indivs[ga->group_state.worstIndex] = ga->indivs[ga->groupSize];

			// 如果变异后出现更优个体
			if (ga->indivs[ga->group_state.bestIndex].fitness >= ga->indivs[ga->groupSize].fitness)
			{
				ga->indivs[ga->groupSize] = ga->indivs[ga->group_state.bestIndex];
			}

			// 记录最优解
			ga->bestIndivs.push_back(ga->indivs[ga->groupSize]);
			
			// 再次寻找子代最差和最优个体
            std::tie(ga->group_state.worstIndex, ga->group_state.bestIndex) = ga->selectPolarIndivs(0, 1);

			// 更新fitArrayCache数组
			ga->updateFitArrayCache();

			// 迭代次数加一
			(ga->group_state).nGene++;

			sel_flag.set_all(false);
			selectReady = false;
			crossReady = true;
		}
		else
		{
			crossReady = false;
		}
	}
}

#endif
