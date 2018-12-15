#ifndef _PSO_H_
#define _PSO_H_

#include <chrono>
#include <vector>
#include "opt_time.h"
#include "pso_individual.h"
#include "pso_state.h"
#include "index_seq.h"
#include "pso_info.h"

namespace opt
{
	template<class F> class PSO;

	// PSO类，需提供适应度函数类型
	template<class R, class... Args>
	class PSO<R(Args...)>
	{
		using Bound = std::initializer_list<std::initializer_list<double>>;

	private:
		std::size_t groupSize;                                                // 粒子数量
		std::size_t nVars;                                                    // 适应度函数包含的变量个数
		PSO_Individual* indivs;                                               // 粒子个体
		R(*fitFunc)(Args...);                                                 // 适应度函数指针
		std::function<std::size_t(std::size_t)> reweight;                     // 种群数量动态调整
		std::function<void(const PSO_Info&)> monitor;                         // 外部监听器
		double(*bound)[2];                                                    // 每个变量的区间, 以数组指针表示
		std::vector<PSO_Individual> bestIndivs;                               // 记录每次迭代的最优个体
		PSO_State group_state;                                                // 粒子群状态

	public:
		PSO(R(*f)(Args...), const std::size_t size = 1000);                   // 构造函数，构造一个种群，需提供适应度函数及种群数量
		PSO(PSO<R(Args...)>& other);                                          // 拷贝构造
		PSO(PSO<R(Args...)>&& other);                                         // 移动构造
		PSO<R(Args...)>& operator=(const PSO<R(Args...)>& other) = delete;
		PSO<R(Args...)>& operator=(PSO<R(Args...)>&& other) = delete;
		~PSO();

		void setBoundary(double(*b)[2]);                                      // 设置变量区间
		void setBoundary(const Bound& b);                                     // 设置变量区间
		void setMaxGeneration(const unsigned int N);                          // 设置最大迭代次数
		void setMaxRuntime(const Second& time);                               // 设置最大运行时间（秒）
		void setStopTol(const long double t, const unsigned int N = 5);       // 设置最优解停止误差
		void setMonitor(const std::function<void(const PSO_Info&)>&);         // 设置外部监听函数
		void setReweight(const std::function<std::size_t(std::size_t)>&);     // 种群数量调整函数

		int getNVars()const;                                                  // 获取种群变量个数
		int getGeneration()const;                                             // 获得当前种群代数
		int getGroupSize()const;                                              // 获得当前种群个体数量
		std::vector<PSO_Individual> getBestIndivs();                          // 获取历次迭代的最优解
		int getStopCode();                                                    // 获取Stop Code

		void initGroup(const std::vector<PSO_Individual>& indivs = std::vector<PSO_Individual>());    // 初始化粒子群
		bool start();                                                         // 开始迭代进化
		void pause();                                                         // 暂停(为保证数据一致性，需在一次完整迭代后pause)
		void proceed();                                                       // 继续迭代
		void kill();                                                          // 结束迭代
		PSO<R(Args...)> clone();                                              // 克隆当前种群

	private:
		template<std::size_t... I>
		R callFitFunc(double* args, const opt::index_seq<I...>&);              // 调用适应度函数数
	};

	/******************************************* 构造与析构 ***********************************************************/
	// 构造函数，需提供供适应度函数及粒子群数量
	template<class R, class... Args>
	PSO<R(Args...)>::PSO(R(*f)(Args...), const std::size_t size)
		:groupSize(size),
		nVars(sizeof...(Args)),
		indivs(new PSO_Individual[groupSize]),
		fitFunc(f),
		bound(nullptr)
	{
		for (std::size_t i = 0; i < groupSize; i++)
		{
			indivs[i] = PSO_Individual(nVars);
		}

		// 变量区间
		bound = new double[nVars][2];
	}

	// 复制构造
	template<class R, class... Args>
	PSO<R(Args...)>::PSO(PSO<R(Args...)>& other)
		:groupSize(other.groupSize),
		nVars(other.nVars),
		indivs(new PSO_Individual[groupSize]),
		fitFunc(other.fitFunc),
		monitor(other.monitor),
		reweight(other.reweight),
		bound(nullptr),
		bestIndivs()
	{
		other.pause();

		for (std::size_t i = 0; i < groupSize; i++)
		{
			indivs[i] = (other.indivs)[i];
		}

		// 变量区间
		bound = new double[nVars][2];
		for (std::size_t i = 0; i < nVars; i++)
		{
			bound[i][0] = (other.bound)[i][0];
			bound[i][1] = (other.bound)[i][1];
		}

		this->bestIndivs = other.bestIndivs;

		other.proceed();
	}

	// 移动构造
	template<class R, class... Args>
	PSO<R(Args...)>::PSO(PSO<R(Args...)>&& other)
		:groupSize(other.groupSize),
		nVars(other.nVars),
		indivs(other.indivs),
		fitFunc(other.fitFunc),
		monitor(other.monitor),
		reweight(other.reweight),
		bound(other.bound),
		bestIndivs()
	{
		other.kill();        // 终止种群迭代

		other.indivs = nullptr;
		other.bound = nullptr;

		this->bestIndivs = std::move(other.bestIndivs);
	}

	// 析构函数
	template<class R, class... Args>
	PSO<R(Args...)>::~PSO()
	{
		delete[] indivs;
		delete[] bound;
	}
	
	/**************************************** Setter ************************************************************/
	// 设置所有变量区间,传入初始化列表
	template<class R, class... Args>
	void PSO<R(Args...)>::setBoundary(const Bound& b)
	{
		const std::size_t len = b.size();
		for (std::size_t i = 0; i < len; i++)
		{
			bound[i][0] = *((*(b.begin() + i)).begin());         // 下边界
			bound[i][1] = *((*(b.begin() + i)).begin() + 1);     // 上边界
		}
		group_state.setBoundFlag = true;
	}

	// 设置所有变量区间,传入数组指针
	template<class R, class... Args>
	void PSO<R(Args...)>::setBoundary(double(*b)[2])
	{
		for (int i = 0; i < nVars; i++)
		{
			bound[i][0] = b[i][0];       // 下边界
			bound[i][1] = b[i][1];       // 上边界
		}
		group_state.setBoundFlag = true;
	}

	// 设置最大迭代次数
	template<class R, class... Args>
	void PSO<R(Args...)>::setMaxGeneration(const unsigned int N)
	{
		group_state.setMaxGeneFlag = true;
		group_state.maxGene = N;
	}

	// 设置最大运行时间（秒）
	template<class R, class... Args>
	void PSO<R(Args...)>::setMaxRuntime(const Second& time)
	{
		group_state.setRuntimeFlag = true;
		group_state.maxRuntime = time;
	}

	// 设置最优解停止误差
	template<class R, class... Args>
	void PSO<R(Args...)>::setStopTol(const long double t, const unsigned int N)
	{
		group_state.converCount = N;
		group_state.setStopTolFlag = true;
		group_state.stopTol = t;
	}

	// 设置外部监听函数
	template<class R, class... Args>
	void PSO<R(Args...)>::setMonitor(const std::function<void(const PSO_Info&)>& func)
	{
		this->monitor = func;
	}

	// 种群数量调整函数
	template<class R, class... Args>
	void PSO<R(Args...)>::setReweight(const std::function<std::size_t(std::size_t)>& func)
	{
		this->reweight = func;
	}

	/**************************************Getter*********************************************************************/
	// 获取种群变量个数
	template<class R, class... Args>
	int PSO<R(Args...)>::getNVars()const
	{
		return nVars;
	}

	// 获取前种群代数
	template<class R, class... Args>
	int PSO<R(Args...)>::getGeneration()const
	{
		return group_state.nGene;
	}

	// 获取当前种群个体数量
	template<class R, class... Args>
	int PSO<R(Args...)>::getGroupSize()const
	{
		return groupSize;
	}

	// 获取历次迭代的最优解
	template<class R, class... Args>
	std::vector<PSO_Individual> PSO<R(Args...)>::getBestIndivs()
	{
		return bestIndivs;
	}

	// 获取停止代码
	template<class R, class... Args>
	int PSO<R(Args...)>::getStopCode()
	{
		return group_state.stopCode;
	}

	// 调用适应度函数
	template<class R, class... Args>
	template<std::size_t... I>
	inline R PSO<R(Args...)>::callFitFunc(double* args, const opt::index_seq<I...>&)
	{
		return fitFunc(args[I]...);
	}

}

#endif