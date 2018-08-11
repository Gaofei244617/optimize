#ifndef _Genetic_Algorithm_
#define _Genetic_Algorithm_

#include <string>
#include <utility>
#include <chrono>
#include <tuple>
#include <vector>
#include <cmath>
#include <thread>
#include "Individual.h"
#include "Range_Random.h"
#include "GA_Time.h"
#include "GA_State.h"
#include "cross_factor.h"
#include "mutate_factor.h"
#include "GA_ThreadPool.h"

/*************Test*******************/
#include <iostream>
#include "out_bool.h"
/************************************/

namespace opt
{
	template<class F>
	class GAGroup;

	// 种群类，需提供适应度函数类型
	template<class R, class... Args>
	class GAGroup<R(Args...)>
	{
		using GenBound = std::initializer_list<std::initializer_list<double>>;

	private:
		std::string name;                                                     // 种群名称
		int groupSize;                                                        // 初始种群个体数量，default = 1000；
		const int nVars;                                                      // 适应度函数包含的变量个数    
		Individual* indivs;                                                   // 种群个体, 指向groupSize + 1个个体数组(最后一个存放最优个体)
		Individual* tempIndivs;                                               // 子代个体缓存区
		R(*fitFunc)(Args...);                                                 // 适应度函数指针
		double(*bound)[2];                                                    // 每个变量(基因)的区间, 以数组指针表示
		double* fitArrayCache;                                                // fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		double mutateProb;                                                    // 个体变异概率,默认p = 0.1
		double crossProb;                                                     // 个体交叉概率, 默认p = 0.6
		
		GA_State state;                                                       // 种群状态
		GA_ThreadState thread_state;                                          // 

	public:
		std::vector<Individual> bestIndivs;                                   // 记录每次迭代的最优个体	

		GAGroup(R(*f)(Args...), const int size = 1000);                       // 构造函数，构造一个种群，需提供适应度函数及种群数量
		GAGroup(const GAGroup<R(Args...)>& other);                            // 拷贝构造
		GAGroup(GAGroup<R(Args...)>&& other);                                 // 移动构造
		~GAGroup();

		void setName(const std::string& str);                                 // 设置种群名称
		void setBoundary(double(*b)[2]);                                      // 设置变量区间
		void setBoundary(const GenBound& b);                                  // 设置变量区间
		void setMaxGeneration(int N);                                         // 设置最大迭代次数
		void setMaxRuntime(const Second& time);                               // 设置最大运行时间（秒）
		void setStopTol(long double t);                                       // 设置最优解停止误差
		void setMutateProb(double p);                                         // 设置基因变异概率
		void setCrossProb(double p);                                          // 设置交叉概率
		void setThreadNum(const int NUM);                                     // 设置并行计算的线程数，默认为1

		const std::string getName()const;                                          // 获取种群名称
		int getNVars()const;                                                  // 获取种群变量个数 
		int getGeneration()const;                                             // 获得当前种群代数
		int getGroupSize()const;                                              // 获得当前种群个体数量	
		const Individual& getIndivByIndex(int index)const;                    // 获取种群个体（通过下标）

		void test();                                                       ///////////////////////////////

		bool start();                                                         // 开始进化
		bool single_start();
		bool pause();                                                         // 停止进化
		bool proceed();                                                       // 继续迭代                                                    

	private:
		void initGroup();                                                     // 初始化种群个体
		bool crossover();                                                     // 交叉
		void cross_thread(const int seq);
		void mut_thread(const int seq);
		void select();
		bool mutate();                                                        // 变异
		void run();                                                           // 进化迭代
		std::pair<int, int> findPolarIndex();                                 // 寻找最差和最好的个体位置,返回<worst, best>
		int randomPickIndiv();                                                // 依据个体的fitness随机选取一个个体，返回个体位置
		template<std::size_t... I>
		R callFitFunc(double* args, const std::index_sequence<I...>&);        // 调用适应度函数
		bool flushStopFlag();                                                 // 判断是否结束迭代
	};

	/******************************************* Helper Functions*****************************************************/
	template<class R, class... Args>
	GAGroup<R(Args...)> createGAGroup(R(*func)(Args...), const int N = 1000)
	{
		return GAGroup<R(Args...)>(func, N);
	}
	/*****************************************************************************************************************/

	/******************************************* 构造与析构 ***********************************************************/
	// 构造函数，需提供供适应度函数及种群数量
	template<class R, class... Args>
	GAGroup<R(Args...)>::GAGroup(R(*f)(Args...), const int size)
		:name("None"),
		groupSize(size + size % 2),
		nVars(sizeof...(Args)),
		fitFunc(f),
		indivs(nullptr),
		tempIndivs(nullptr),
		bound(nullptr),
		fitArrayCache(nullptr),
		mutateProb(0.1),
		crossProb(0.6)
	{
		// 分配个体内存，但不构造个体(Indivadual无默认构造),最后一个位置用于存放最优个体
		this->indivs = static_cast<Individual*>(::operator new(sizeof(Individual) * (groupSize + 1)));
		this->tempIndivs = static_cast<Individual*>(::operator new(sizeof(Individual) * (groupSize + 1)));

		// 在已分配的内存上构造个体对象
		for (int i = 0; i < groupSize + 1; i++)
		{
			new(indivs + i) Individual(nVars);
			new(tempIndivs + i) Individual(nVars);
		}

		// 变量区间
		bound = new double[nVars][2];

		// fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		fitArrayCache = new double[groupSize + 1]();
	}

	// 复制构造
	template<class R, class... Args>
	GAGroup<R(Args...)>::GAGroup(const GAGroup<R(Args...)>& other)
		:name(other.name),
		groupSize(other.groupSize),
		nVars(other.nVars),
		fitFunc(other.fitFunc),
		mutateProb(other.mutateProb),
		crossProb(other.crossProb),
		indivs(nullptr),
		tempIndivs(nullptr),
		bound(nullptr),
		fitArrayCache(nullptr),
		state(other.state)
	{
		// 分配个体内存，但不构造个体(Indivadual无默认构造)
		this->indivs = static_cast<Individual*>(::operator new(sizeof(Individual) * (groupSize + 1)));
		this->tempIndivs = static_cast<Individual*>(::operator new(sizeof(Individual) * (groupSize + 1)));

		// 在已分配的内存上构造个体对象
		for (int i = 0; i < groupSize + 1; i++)
		{
			new(indivs + i) Individual(*(other.indivs + i));
			new(tempIndivs + i) Individual(*(other.tempIndivs + i));
		}

		// 变量区间
		bound = new double[nVars][2];
		for (int i = 0; i < nVars; i++)
		{
			bound[i][0] = (other.bound)[i][0];
			bound[i][1] = (other.bound)[i][1];
		}

		// fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		fitArrayCache = new double[groupSize + 1]();
		for (int i = 0; i < groupSize + 1; i++)
		{
			fitArrayCache[i] = (other.fitArrayCache)[i];
		}
	}

	// 移动构造
	template<class R, class... Args>
	GAGroup<R(Args...)>::GAGroup(GAGroup<R(Args...)>&& other)
		:name(other.name),
		groupSize(other.groupSize),
		nVars(other.nVars),
		indivs(other.indivs),
		tempIndivs(other.tempIndivs),
		fitFunc(other.fitFunc),
		bound(other.bound),
		fitArrayCache(other.fitArrayCache),
		mutateProb(other.mutateProb),
		crossProb(other.crossProb),
		state(other.state)
	{
		other.indivs = nullptr;
		other.tempIndivs = nullptr;
		other.bound = nullptr;
		other.fitArrayCache = nullptr;
	}

	// 析构函数
	template<class R, class... Args>
	GAGroup<R(Args...)>::~GAGroup()
	{
		// 析构indivs数组
		if (indivs != nullptr)
		{
			// 析构种群中的每一个个体（共有groupSize + 1个个体）
			for (int i = 0; i < groupSize + 1; i++)
			{
				(indivs + i) -> ~Individual();
			}
			// 释放个体对象指针
			::operator delete(indivs);
		}

		// 析构tempIndivs数组
		if (tempIndivs != nullptr)
		{
			// 析构种群中的每一个个体（共有groupSize + 1个个体）
			for (int i = 0; i < groupSize + 1; i++)
			{
				(tempIndivs + i) -> ~Individual();
			}
			// 释放个体对象指针
			::operator delete(tempIndivs);
		}

		// 释放Boundary数组指针
		delete[] bound;

		// 释放缓存数组
		delete[] fitArrayCache;
	}
	/*****************************************************************************************************************/

	/**************************************** Setter ************************************************************/
	// 设置种群名称
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setName(const std::string& str)
	{
		this->name = str;
	}

	// 设置所有变量区间,传入初始化列表
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setBoundary(const GenBound& b)
	{
		const size_t len = b.size();
		for (size_t i = 0; i < len; i++)
		{
			bound[i][0] = *((*(b.begin() + i)).begin());         // 下边界
			bound[i][1] = *((*(b.begin() + i)).begin() + 1);     // 上边界
		}
		state.setBoundFlag = true;
	}

	// 设置所有变量区间,传入数组指针
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setBoundary(double(*b)[2])
	{
		for (int i = 0; i < nVars; i++)
		{
			bound[i][0] = b[i][0];       // 下边界
			bound[i][1] = b[i][1];       // 上边界
		}
		state.setBoundFlag = true;
	}

	// 设置最大迭代次数
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setMaxGeneration(int N)
	{
		state.setMaxGeneFlag = true;
		state.maxGene = N;
	}

	// 设置最大运行时间（秒）
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setMaxRuntime(const Second& time)
	{
		state.setRuntimeFlag = true;
		state.maxRuntime = time;
	}

	// 设置最优解停止误差
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setStopTol(long double t)
	{
		state.setStopTolFlag = true;
		state.stopTol = t;
	}

	// 设置基因变异概率
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setMutateProb(double p)
	{
		mutateProb = p;
	}

	// 设置基因交叉概率
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setCrossProb(double p)
	{
		crossProb = p;
	}

	// 设置并行计算的线程数，默认为1
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setThreadNum(const int NUM)
	{
		if (NUM >= 1)
		{
			thread_state.threadNum = NUM;
			thread_state.cross_flag.set_length(NUM);
			thread_state.mut_flag.set_length(NUM);
		}
	}
	/*****************************************************************************************************************/

	/**************************************Getter*********************************************************************/
	// 获取种群名称
	template<class R, class... Args>
	const std::string GAGroup<R(Args...)>::getName()const
	{
		return this->name;
	}

	// 获取种群变量个数 
	template<class R, class... Args>
	int GAGroup<R(Args...)>::getNVars()const
	{
		return nVars;
	}

	// 获取前种群代数
	template<class R, class... Args>
	int GAGroup<R(Args...)>::getGeneration()const
	{
		return state.nGene;
	}

	// 获取当前种群个体数量
	template<class R, class... Args>
	int GAGroup<R(Args...)>::getGroupSize()const
	{
		return groupSize;
	}

	// 获取种群个体（通过下标）
	template<class R, class... Args>
	const Individual& GAGroup<R(Args...)>::getIndivByIndex(int index)const
	{
		return *(indivs + index);
	}

	/************************************************************************/
	/************************************************************************/
	/************************************************************************/
	// 获取最优个体
	template<class R, class... Args>
	void GAGroup<R(Args...)>::test()
	{
		using namespace std;

		// 优化结果的基因与fitness
		std::vector<double> vec;
		for (int i = 0; i < nVars; i++)
		{
			vec.push_back(indivs[groupSize].vars[i]);
		}
		vec.push_back(indivs[groupSize].fitness);

		// 输出优化结果
		std::cout << "Stop code: " << state.stopCode << std::endl;
		cout << vec[0] << endl;
		cout << vec[1] << endl;
		cout << "fitness: " << vec[2] << endl;
		cout << "**************************************" << endl;

		// 输出子代最优解进化过程
		cout << bestIndivs[0].fitness << endl;
		cout << endl;
		for (size_t i = 1; i < bestIndivs.size(); i++)
		{
			if (bestIndivs[i].fitness != bestIndivs[i - 1].fitness)
			{
				cout << bestIndivs[i].fitness << endl;
				cout << endl;
			}
		}
	}

	/******************************************************************************************************************/

	// 开始进化
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::start()
	{
		// 初始化种群
		initGroup();
		//std::cout << "aaaaa" << std::endl;

		std::vector<std::thread> vec_cross;
		std::vector<std::thread> vec_mut;

		if (state.runable())
		{
			//state.startTime = std::chrono::steady_clock::now();
			//// 单独的线程里进行种群迭代
			//std::thread t(&GAGroup<R(Args...)>::run, this);
			//t.join();
			//return true;

			// 构造种群交叉线程池
			for (int i = 0; i < thread_state.threadNum; i++)
			{
				vec_cross.emplace_back(&GAGroup<R(Args...)>::cross_thread, this, i);
			}

			// 构造种群变异线程池
			for (int i = 0; i < thread_state.threadNum; i++)
			{
				vec_mut.emplace_back(&GAGroup<R(Args...)>::mut_thread, this, i);
			}

			// 环境选择,  主线程中执行
			select();

			// 阻塞当前线程
			for (int i = 0; i < thread_state.threadNum; i++)
			{
				vec_cross[i].join();
				vec_mut[i].join();
			}

			return true;
		}
		return false;
	}

	// 单线程执行
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::single_start()
	{
		initGroup();

		if (state.runable())
		{
			state.startTime = std::chrono::steady_clock::now();
			std::thread t(&GAGroup<R(Args...)>::run, this);
			t.join();
			return true;
		}
		return false;
	}

	// 暂停进化
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::pause()
	{
		state.stopFlag = true;
		state.stopCode = 3;
		return true;
	}

	// 继续迭代
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::proceed()
	{
		state.stopFlag = false;
		state.stopCode = -1;
		std::thread t(&GAGroup<R(Args...)>::run, this);
		t.join();
		return true;
	}

	/******************************************* Private Functions *****************************************************/
	// 初始化种群
	template<class R, class... Args>
	void GAGroup<R(Args...)>::initGroup()
	{
		// 初始化前需要保证已设置基因变量区间
		if (state.setBoundFlag)
		{
			// 1.初始化每一个个体
			for (int i = 0; i < groupSize; i++)
			{
				// 设置个体初始基因
				for (int j = 0; j < nVars; j++)
				{
					indivs[i].vars[j] = random_real(bound[j][0], bound[j][1]);
				}
				// 个体适应度
				indivs[i].fitness = callFitFunc(indivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
			}

			// 2.寻找最差和最好适应度个体
			std::tie(state.worstIndex, state.bestIndex) = this->findPolarIndex();
			indivs[groupSize] = indivs[state.bestIndex];                 // indivs数组末尾位置用于存储最优解
			bestIndivs.push_back(indivs[groupSize]);                     // 记录最优解

			// 3.初始化fitArrayCache数组
			// fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
			for (int i = 1; i < groupSize; i++)
			{
				fitArrayCache[i] = fitArrayCache[i - 1] + (indivs[i - 1].fitness - indivs[state.worstIndex].fitness);
			}

			//4.设置初始化标志位
			state.initFlag = true;

			thread_state.crossReady = true;
			thread_state.mutReady = false;
			thread_state.selectReady = false;
		}
		// throw std::string("Fitness function or genes boundary is not set.");
	}

	template<class R, class... Args>
	bool GAGroup<R(Args...)>::crossover()
	{
		int Index_M = 0;             // 父个体
		int Index_F = 0;             // 母个体

		double rand_cross = 0;       // 随机数缓存    

		// 随机选取groupSize/2对个体交叉, fitness越大，被选中的概率越大
		for (int i = 0; i < groupSize; i += 2)
		{
			// 随机选取两个个体作为父代
			Index_M = randomPickIndiv();
			Index_F = randomPickIndiv();

			// 生成子代个体, 存放于缓存数组tempIndivs中
			for (int j = 0; j < nVars; j++)
			{
				// 依据交叉概率决定基因是否进行交叉
				rand_cross = random_real(0, 1);
				if (rand_cross <= crossProb)   // 交叉
				{
					// 基因交叉产生子代基因
					std::tie(tempIndivs[i].vars[j], tempIndivs[i + 1].vars[j]) = cross_SBX(indivs[Index_M].vars[j], indivs[Index_F].vars[j]);

					// 子代基因是否超过边界
					if (tempIndivs[i].vars[j] < bound[j][0] || tempIndivs[i].vars[j] > bound[j][1])
					{
						tempIndivs[i].vars[j] = random_real(bound[j][0], bound[j][1]);
					}
					if (tempIndivs[i + 1].vars[j] < bound[j][0] || tempIndivs[i + 1].vars[j] > bound[j][1])
					{
						tempIndivs[i + 1].vars[j] = random_real(bound[j][0], bound[j][1]);
					}
				}
				else    // 不交叉
				{
					// 继承父代基因，基因不发生交叉
					tempIndivs[i].vars[j] = indivs[Index_M].vars[j];
					tempIndivs[i + 1].vars[j] = indivs[Index_F].vars[j];
				}
			}

			// 计算子代个体适应度
			tempIndivs[i].fitness = callFitFunc(tempIndivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
			tempIndivs[i + 1].fitness = callFitFunc(tempIndivs[i + 1].vars, std::make_index_sequence<sizeof...(Args)>());
		}

		// 父代最优个体遗传到子代
		tempIndivs[groupSize] = indivs[groupSize];

		// 交换子代个体缓存区和父代种群个体存放区的指针
		Individual* temp = indivs;
		indivs = tempIndivs;
		tempIndivs = temp;

		// 寻找子代最差和最优个体
		std::tie(state.worstIndex, state.bestIndex) = this->findPolarIndex();

		// 如果子代出现更优个体
		if (indivs[state.bestIndex].fitness >= indivs[groupSize].fitness)
		{
			indivs[groupSize] = indivs[state.bestIndex];
		}

		return true;
	}

	// 变异
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::mutate()
	{
		double rand_num = 0; // 随机数

		// 每个个体的变异,最优解个体不发生变异
		for (int i = 0; i < groupSize; i++)
		{
			if (i != state.bestIndex)
			{
				// 每个基因的变异
				for (int j = 0; j < nVars; j++)
				{
					rand_num = random_real(0, 1);
					if (rand_num < mutateProb)
					{
						// 基因变异
						indivs[i].vars[j] = mutate_PM(indivs[i].vars[j], bound[j][0], bound[j][1]);
					}
				}
				// 个体适应度
				indivs[i].fitness = callFitFunc(indivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
			}
		}

		// 更新fitArrayCache数组
		// fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		for (int i = 1; i < groupSize; i++)
		{
			fitArrayCache[i] = fitArrayCache[i - 1] + (indivs[i - 1].fitness - indivs[state.worstIndex].fitness);
		}

		// 寻找子代最差和最优个体
		std::tie(state.worstIndex, state.bestIndex) = this->findPolarIndex();

		// 如果变异后出现更优个体
		if (indivs[state.bestIndex].fitness >= indivs[groupSize].fitness)
		{
			indivs[groupSize] = indivs[state.bestIndex];
		}

		return true;
	}

	// 进化迭代
	template<class R, class... Args>
	void GAGroup<R(Args...)>::run()
	{
		while (true)
		{
			crossover();
			mutate();

			// 记录最优解
			bestIndivs.push_back(indivs[groupSize]);        

			// 迭代次数加一
			state.nGene++;
			// 记录当前时间
			state.nowTime = std::chrono::steady_clock::now();

			// 判断最优解相较于上一次的波动值
			if (abs(bestIndivs[state.nGene - 1].fitness - bestIndivs[state.nGene].fitness) <= state.stopTol)
			{
				state.count++;
			}
			else
			{
				state.count = 0;
			}

			// 判断是否到达停止条件
			if (flushStopFlag())
			{
				return;
			}
		}
	}

	// 寻找最差和最好的个体位置,返回pair<worst, best>
	template<class R, class... Args>
	std::pair<int, int> GAGroup<R(Args...)>::findPolarIndex()
	{
		int worst = 0;
		int best = 0;

		// 遍历indivs的fitness，寻找fitness的极值，记录极值位置
		for (int i = 1; i < groupSize; i++)
		{
			if (indivs[i].fitness < indivs[worst].fitness)
			{
				worst = i;
			}

			if (indivs[i].fitness > indivs[best].fitness)
			{
				best = i;
			}
		}

		return std::make_pair(worst, best);
	}

	// 依据个体的fitness随机选取一个个体，返回个体位置
	template<class R, class... Args>
	int GAGroup<R(Args...)>::randomPickIndiv()
	{
		// fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		double temp = random_real(0, fitArrayCache[groupSize]);

		int lowIndex = 0;
		int upIndex = groupSize;
		int tempIndex;

		// 二分法查找随机个体位置,时间复杂度log_n
		// 通过逐步移动lowIndex和upIndex代表的下标位置，每次将搜索区间减小一半
		while (upIndex - lowIndex > 1)
		{
			tempIndex = (lowIndex + upIndex) / 2 + (lowIndex + upIndex) % 2;
			if (temp <= fitArrayCache[tempIndex])
			{
				upIndex = tempIndex;
			}
			else
			{
				lowIndex = tempIndex;
			}
		}
		return lowIndex;  // up = 1, low = 0 时随机个体是indivs[0]
	}

	template<class R, class... Args>
	template<std::size_t... I>
	inline R GAGroup<R(Args...)>::callFitFunc(double* args, const std::index_sequence<I...>&)
	{
		return fitFunc(args[I]...);
	}

	// 判断是否结束迭代
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::flushStopFlag()
	{
		// 是否达到最大迭代次数
		if (state.setMaxGeneFlag && state.nGene >= state.maxGene)
		{
			state.stopCode = 1;
		}

		// 是否达到最大迭代时间
		std::chrono::duration<double> evolTime = state.nowTime - state.startTime;
		if (state.setRuntimeFlag && evolTime.count() >= state.maxRuntime.value)
		{
			state.stopCode = 2;
		}

		// 是否连续5次最优解的波动小于停止误差
		if (state.setStopTolFlag && state.count == 5)
		{
			state.stopCode = 0;
		}

		// 根据stop code判断种群是否停止迭代
		if (state.stopCode != -1)
		{
			state.stopFlag = true;
		}

		return state.stopFlag;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	template<class R, class... Args>
	void GAGroup<R(Args...)>::cross_thread(const int seq)
	{
		//std::condition_variable cv_cross;
		int Index_M = 0;
		int Index_F = 0;
		double rand_cross = 0;

		while (true)
		{
			{
				std::unique_lock<std::mutex> lck(thread_state.mtx);
				//std::cout << "---* Cross before: "<< seq << std::endl;
				//std::cout << "CrossReady: " << out_bool(thread_state.crossReady) << endl;
				//std::cout << "MutateReady: " << out_bool(thread_state.mutReady) << endl;
				//std::cout << "SelecReady: " << out_bool(thread_state.selectReady) << endl;
				thread_state.cv.wait(lck, [this, &seq]() {return !(this->thread_state).cross_flag[seq] && (this->thread_state).crossReady || this->flushStopFlag(); });
			}// 离开作用域, thread_state.mtx的unlock()方法自动执行
			
			//std::this_thread::sleep_for(1s);

			 // 是否停止迭代
			if (flushStopFlag())
			{
				//cout << "Exit Cross..." << endl;
				return;
			}

			for (int i = seq; i < groupSize; i += 2 * thread_state.threadNum)
			{
				// 随机选取两个个体作为父代
				Index_M = randomPickIndiv();
				Index_F = randomPickIndiv();

				// 生成子代个体, 存放于缓存数组tempIndivs中
				for (int j = 0; j < nVars; j++)
				{
					// 依据交叉概率决定基因是否进行交叉
					rand_cross = random_real(0, 1);
					if (rand_cross <= crossProb)   // 交叉
					{
						// 基因交叉产生子代基因
						std::tie(tempIndivs[i].vars[j], tempIndivs[i + 1].vars[j]) = cross_SBX(indivs[Index_M].vars[j], indivs[Index_F].vars[j]);

						// 子代基因是否超过边界
						if (tempIndivs[i].vars[j] < bound[j][0] || tempIndivs[i].vars[j] > bound[j][1])
						{
							tempIndivs[i].vars[j] = random_real(bound[j][0], bound[j][1]);
						}
						if (tempIndivs[i + 1].vars[j] < bound[j][0] || tempIndivs[i + 1].vars[j] > bound[j][1])
						{
							tempIndivs[i + 1].vars[j] = random_real(bound[j][0], bound[j][1]);
						}
					}
					else    // 不交叉
					{
						// 继承父代基因，基因不发生交叉
						tempIndivs[i].vars[j] = indivs[Index_M].vars[j];
						tempIndivs[i + 1].vars[j] = indivs[Index_F].vars[j];
					}
				}

				// 计算子代个体适应度
				tempIndivs[i].fitness = callFitFunc(tempIndivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
				tempIndivs[i + 1].fitness = callFitFunc(tempIndivs[i + 1].vars, std::make_index_sequence<sizeof...(Args)>());
			}

			// 线程同步
			thread_state.mtx.lock();
			thread_state.cross_flag[seq] = true;
			thread_state.selectReady = false;
			//std::cout << "---* Cross after: " << seq << "\n" << std::endl;
			if (thread_state.cross_flag.is_all_true())
			{
				thread_state.cross_flag.set_all(false);
				thread_state.crossReady = false;
				thread_state.mutReady = true;
				//std::cout << "--------Cross finish once in thread: --------" << seq << std::endl;
			}
			else
			{
				thread_state.mutReady = false;
			}
			//std::cout << "\n++ Change in cross: " << seq << endl;
			//std::cout << "CrossReady: " << out_bool(thread_state.crossReady) << endl;
			//std::cout << "MutateReady: " << out_bool(thread_state.mutReady) << endl;
			//std::cout << "SelecReady: " << out_bool(thread_state.selectReady) << endl;
			thread_state.mtx.unlock();
			thread_state.cv.notify_all();
		}
	}

	////////////////////////////////////////////////////////
	template<class R, class... Args>
	void GAGroup<R(Args...)>::mut_thread(const int seq)
	{
		//std::condition_variable cv_mut;  // 条件变量
		double rand_num = 0; // 随机数

		while (true)
		{
			{
				// 临界区
				std::unique_lock<std::mutex> lck(thread_state.mtx);
				//std::cout << "---# Mutate before: " << seq << std::endl;
				//std::cout << "CrossReady: " << out_bool(thread_state.crossReady) << endl;
				//std::cout << "MutateReady: " << out_bool(thread_state.mutReady) << endl;
				//std::cout << "SelecReady: " << out_bool(thread_state.selectReady) << endl;
				thread_state.cv.wait(lck, [this, &seq]() {return !(this->thread_state).mut_flag[seq] && (this->thread_state).mutReady || this->flushStopFlag(); });
			}

			//std::this_thread::sleep_for(1s);

			// 是否停止迭代
			if (flushStopFlag())
			{
				//cout << "Exit Mutate..." << endl;
				return;
			}

			// 每个个体变异(最优解个体不发生变异)
			for (int i = seq; i < groupSize; i += thread_state.threadNum)
			{
				if (i != state.bestIndex)
				{
					// 每个基因的变异
					for (int j = 0; j < nVars; j++)
					{
						rand_num = random_real(0, 1);
						if (rand_num < mutateProb)
						{
							// 基因变异
							indivs[i].vars[j] = mutate_PM(indivs[i].vars[j], bound[j][0], bound[j][1]);
						}
					}
					// 个体适应度
					indivs[i].fitness = callFitFunc(indivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
				}
			}

			// 线程同步
			thread_state.mtx.lock();
			thread_state.mut_flag[seq] = true;
			thread_state.crossReady = false;
			//std::cout << "---# Mutate after: " << seq << "\n" << std::endl;
			if (thread_state.mut_flag.is_all_true())
			{
				thread_state.mut_flag.set_all(false);
				thread_state.mutReady = false;
				thread_state.selectReady = true;
				//std::cout << "--------Mutate finish once--------" << std::endl;
			}
			else
			{
				thread_state.selectReady = false;
			}
			//std::cout << "\n++ Change in mutate: " << seq << endl;
			//std::cout << "CrossReady: " << out_bool(thread_state.crossReady) << endl;
			//std::cout << "MutateReady: " << out_bool(thread_state.mutReady) << endl;
			//std::cout << "SelecReady: " << out_bool(thread_state.selectReady) << endl;
			thread_state.mtx.unlock();
			thread_state.cv.notify_all();
		}
	}

	// 选择
	template<class R, class ...Args>
	void GAGroup<R(Args...)>::select()
	{
		//std::condition_variable cv_select;

		while (true)
		{
			// 临界区
			{
				std::unique_lock<std::mutex> lck(thread_state.mtx);				
				//std::cout << "---^ Select before. " << std::endl;
				//std::cout << "CrossReady: " << out_bool(thread_state.crossReady) << endl;
				//std::cout << "MutateReady: " << out_bool(thread_state.mutReady) << endl;
				//std::cout << "SelecReady: " << out_bool(thread_state.selectReady) << endl;
				thread_state.cv.wait(lck, [this]() {return (this->thread_state).selectReady; });
			}
			
			//std::this_thread::sleep_for(1s);
			// 更新fitArrayCache数组
			// fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
			// fitArrayCache[0]始终为0
			for (int i = 1; i < groupSize; i++)
			{
				fitArrayCache[i] = fitArrayCache[i - 1] + (indivs[i - 1].fitness - indivs[state.worstIndex].fitness);
			}

			// 寻找子代最差和最优个体
			std::tie(state.worstIndex, state.bestIndex) = this->findPolarIndex();

			// 如果变异后出现更优个体
			if (indivs[state.bestIndex].fitness >= indivs[groupSize].fitness)
			{
				indivs[groupSize] = indivs[state.bestIndex];
			}

			// 线程同步
			thread_state.mtx.lock();
			state.nGene++;
			thread_state.selectReady = false;
			thread_state.crossReady = true;
			thread_state.mutReady = false;
			//std::cout << "---^ Select after." << "\n" << std::endl;
			//std::cout << "--------Select finish once--------" << std::endl;
			//std::cout << "第 " << state.nGene << " 代..." << endl;
			//std::cout << "\n++ Change in select: " << endl;
			//std::cout << "CrossReady: " << out_bool(thread_state.crossReady) << endl;
			//std::cout << "MutateReady: " << out_bool(thread_state.mutReady) << endl;
			//std::cout << "SelecReady: " << out_bool(thread_state.selectReady) << endl;
			thread_state.mtx.unlock();
			thread_state.cv.notify_all();
			// 是否停止迭代
			if (flushStopFlag())
			{
				//cout << "Exit Select..." << endl;
				return;
			}
		}
	}
}

#endif