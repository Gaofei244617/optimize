#ifndef _PSO_THREAD_SYNC_
#define _PSO_THREAD_SYNC_

#include <thread>
#include <vector>
#include <utility>
#include <condition_variable>
#include <tuple>
#include "bool_array.h"
#include "pso_info.h"

namespace opt
{
	template<class F> class PSO;

	// GA线程同步器
	template<class R, class... Args>
	class PSOThreadSync
	{
		friend class PSO<R(Args...)>;

	private:
		std::mutex mtx;                                   // 互斥量
		std::condition_variable cv;                       // 条件变量, 用于线程同步
		bool iter_ready;                                  // 交叉操作标志
		std::size_t threadNum;                            // 并行计算使用的线程数
		bool_array iter_flag;                             // 并行计算: 线程状态标志
		PSO<R(Args...)>* pso;                             // 同步器关联的粒子群

	public:
		// 常规构造
		PSOThreadSync(PSO<R(Args...)>* ptr)
			:iter_ready(false),
			threadNum(1),
			iter_flag(threadNum),
			pso(ptr)
		{}

		// 复制构造
		PSOThreadSync(const PSOThreadSync& other, PSO<R(Args...)>* ptr)
			:iter_ready(other.iter_ready),
			threadNum(other.threadNum),
			iter_flag(other.iter_flag),
			pso(ptr)
		{}

		// 移动构造
		PSOThreadSync(PSOThreadSync&& other) = delete;

		// 设置线程数量
		void setThreadNum(const std::size_t N);
	};

	// 设置线程数量
	template<class R, class... Args>
	void PSOThreadSync<R, Args...>::setThreadNum(const std::size_t N)
	{
		threadNum = N;
		iter_flag.set_length(N);
	}
}

#endif
