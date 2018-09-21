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
#include "GA_GroupState.h"
#include "cross_factor.h"
#include "mutate_factor.h"
#include "ga_thread_sync.h"

/*************Test*******************/
#include <iostream>
#include "out_bool.h"
#include <iomanip>
//#include <crtdbg.h>
/************************************/

namespace opt
{
	template<class F> class GAGroup;

	// 种群类，需提供适应度函数类型
	template<class R, class... Args>
	class GAGroup<R(Args...)>
	{
		using GenBound = std::initializer_list<std::initializer_list<double>>;
		friend void GAThreadSync<R, Args...>::select_sync(const int thread_seq);

	private:
		std::string name;                                                     // 种群名称
		int groupSize;                                                        // 初始种群个体数量，default = 1000；
		const int nVars;                                                      // 适应度函数包含的变量个数    
		Individual* indivs;                                                   // 种群个体, 指向groupSize + 1个个体数组(最后一个存放最优个体)
		Individual* tempIndivs;                                               // 子代个体缓存区
		R(*fitFunc)(Args...);                                                 // 适应度函数指针
		double(*bound)[2];                                                    // 每个变量(基因)的区间, 以数组指针表示
		double* fitArrayCache;                                                // 轮盘赌刻度线: fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		double mutateProb;                                                    // 个体变异概率,默认p = 0.1
		double crossProb;                                                     // 个体交叉概率, 默认p = 0.6
		std::vector<Individual> bestIndivs;                                   // 记录每次迭代的最优个体	

		std::vector<std::thread> vec_cross;                                   // 交叉线程组
		std::vector<std::thread> vec_mut;                                     // 变异线程组
		std::vector<std::thread> vec_select;                                  // 环境选择线程组
		std::thread single_thread;                                            // 单线程
		
		GA_GroupState group_state;                                            // 种群状态
		GAThreadSync<R, Args...> thread_sync;                                 // 线程同步器  

	public:
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

		const std::string getName()const;                                     // 获取种群名称
		int getNVars()const;                                                  // 获取种群变量个数 
		int getGeneration()const;                                             // 获得当前种群代数
		int getGroupSize()const;                                              // 获得当前种群个体数量	
		const Individual& getIndivByIndex(const int index)const;              // 获取种群个体（通过下标）

		void test();                                                       ///////////////////////////////

		bool start();                                                         // 开始进化
		void wait_result();                                                   // 等待进化结果
		bool pause();                                                         // 停止进化
		bool proceed();                                                       // 继续迭代                                                    

	private:
		void initGroup();                                                     // 初始化种群个体

		void crossover();                                                     // 交叉(单线程模式)
		void mutate();                                                        // 变异(单线程模式)
		void select();                                                        // 选择(单线程模式)
		void run();                                                           // 进化迭代(单线程模式)
		
		void cross_thread(const int seq);                                     // 交叉(多线程模式)
		void mut_thread(const int seq);                                       // 变异(多线程模式)
		void sel_thread(const int seq);                                       // 选择(多线程模式)
		
		void updateFitArrayCache();                                           // 更新轮盘赌刻度线
		//void updateBestIndiv();                                               // 更新最优个体记录

		std::pair<int, int> selectPolarIndivs(const int seq, const int interval);   // 寻找最差和最好的个体位置,返回<worst, best>
		int randomPickIndiv();                                                // 依据个体的fitness随机选取一个个体，返回个体位置
		template<std::size_t... I>
		R callFitFunc(double* args, const std::index_sequence<I...>&);        // 调用适应度函数
		bool getStopFlag();                                                   // 判断是否结束迭代
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
		crossProb(0.6),
		thread_sync(this)
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

		//std::cout << std::hex << (int)indivs << std::endl;
		//std::cout << std::hex << (int)tempIndivs << std::endl;
		//std::cout << std::hex << (int)fitArrayCache << std::endl;
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
		group_state(other.group_state),
		thread_sync(other.thread_sync, this)
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
		group_state(other.group_state)
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
		group_state.setBoundFlag = true;
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
		group_state.setBoundFlag = true;
	}

	// 设置最大迭代次数
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setMaxGeneration(int N)
	{
		group_state.setMaxGeneFlag = true;
		group_state.maxGene = N;
	}

	// 设置最大运行时间（秒）
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setMaxRuntime(const Second& time)
	{
		group_state.setRuntimeFlag = true;
		group_state.maxRuntime = time;
	}

	// 设置最优解停止误差
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setStopTol(long double t)
	{
		group_state.setStopTolFlag = true;
		group_state.stopTol = t;
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
			thread_sync.setThreadNum(NUM);
		}
		else
		{
			throw std::string("Thread number error.");
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
		return group_state.nGene;
	}

	// 获取当前种群个体数量
	template<class R, class... Args>
	int GAGroup<R(Args...)>::getGroupSize()const
	{
		return groupSize;
	}

	// 获取种群个体（通过下标）
	template<class R, class... Args>
	const Individual& GAGroup<R(Args...)>::getIndivByIndex(const int index)const
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
		std::cout << "Thread number: " << thread_sync.threadNum << std::endl;
		std::cout << "Stop code: " << group_state.stopCode << std::endl;
		cout << vec[0] << endl;
		cout << vec[1] << endl;
		cout << "fitness: " << vec[2] << endl;

		// 输出子代最优解进化过程
		cout << "\n*********************************" << endl;
		cout << "子代最优解进化过程：" << endl;
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

		// 若种群满足进化条件
		if (group_state.runable())
		{
			group_state.startTime = std::chrono::steady_clock::now();
			// 多线程模式
			if (thread_sync.threadNum > 1)
			{
				// 构造种群交叉线程组
				for (int i = 0; i < thread_sync.threadNum; i++)
				{
					vec_cross.emplace_back(&GAGroup<R(Args...)>::cross_thread, this, i);
				}

				// 构造种群变异线程组
				for (int i = 0; i < thread_sync.threadNum; i++)
				{
					vec_mut.emplace_back(&GAGroup<R(Args...)>::mut_thread, this, i);
				}

				// 构造环境选择线程组
				for (int i = 0; i < thread_sync.threadNum; i++)
				{
					vec_select.emplace_back(&GAGroup<R(Args...)>::sel_thread, this, i);
				}

			}
			else // 单线程模式
			{
				single_thread = std::thread(&GAGroup<R(Args...)>::run, this);
			}
			return true;
		}
		return false;
	}

	// 阻塞当前线程, 等待计算结果
	template<class R, class... Args>
	void GAGroup<R(Args...)>::wait_result()
	{
		// 单线程模式
		if (thread_sync.threadNum == 1)
		{
			single_thread.join();
		}

		// 多线程模式
		if (thread_sync.threadNum > 1)
		{
			for (int i = 0; i < thread_sync.threadNum; i++)
			{
				vec_cross[i].join();
				vec_mut[i].join();
				vec_select[i].join();
			}
		}
	}

	//// 暂停进化
	//template<class R, class... Args>
	//bool GAGroup<R(Args...)>::pause()
	//{
	//	group_state.stopFlag = true;
	//	group_state.stopCode = 3;
	//	return true;
	//}

	//// 继续迭代
	//template<class R, class... Args>
	//bool GAGroup<R(Args...)>::proceed()
	//{
	//	group_state.stopFlag = false;
	//	group_state.stopCode = -1;
	//	std::thread t(&GAGroup<R(Args...)>::run, this);
	//	t.join();
	//	return true;
	//}

	/******************************************* Private Functions *****************************************************/
	// 初始化种群     
	template<class R, class... Args>
	void GAGroup<R(Args...)>::initGroup()
	{
		// 初始化前需要保证已设置基因变量区间
		if (group_state.setBoundFlag)
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
			std::tie(group_state.worstIndex, group_state.bestIndex) = this->selectPolarIndivs(0, 1);
			indivs[groupSize] = indivs[group_state.bestIndex];                 // indivs数组末尾位置用于存储最优解
			bestIndivs.push_back(indivs[groupSize]);                     // 记录最优解

			// 3.初始化fitArrayCache数组
			updateFitArrayCache();

			//4.设置初始化标志位
			group_state.initFlag = true;

			thread_sync.crossReady = true;
			thread_sync.mutReady = false;
			thread_sync.selectReady = false;
		}
		else
		{
			throw std::string("Fitness function or genes boundary is not set.");
		}
	}

	// 交叉(单线程模式)
	template<class R, class... Args>
	void GAGroup<R(Args...)>::crossover()
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
				else  // 不交叉
				{
					// 继承父代基因，基因不发生交叉
					tempIndivs[i].vars[j] = indivs[Index_M].vars[j];
					tempIndivs[i + 1].vars[j] = indivs[Index_F].vars[j];
				}
			}

		}

		// 父代最优个体遗传到子代
		tempIndivs[groupSize] = indivs[groupSize];

		// 交换子代个体缓存区和父代个体存放区的指针
		Individual* temp = indivs;
		indivs = tempIndivs;
		tempIndivs = temp;
	}

	// 变异(单线程模式)
	template<class R, class... Args>
	void GAGroup<R(Args...)>::mutate()
	{
		double rand_num = 0; // 随机数

		// 每个个体的变异
		for (int i = 0; i < groupSize; i++)
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
			
			// 计算每个个体适应度
			indivs[i].fitness = callFitFunc(indivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
		}
	}

	// 选择(单线程模式)
	template<class R, class... Args>
	void GAGroup<R(Args...)>::select()
	{
		// 寻找子代最差和最优个体
		std::tie(group_state.worstIndex, group_state.bestIndex) = this->selectPolarIndivs(0, 1);
		
		// 淘汰子代最差个体, 用父代最优个体取代
		indivs[group_state.worstIndex] = indivs[groupSize];

		// 如果子代出现更优个体
		if (indivs[group_state.bestIndex].fitness >= indivs[groupSize].fitness)
		{
			indivs[groupSize] = indivs[group_state.bestIndex];
		}

		// 再次寻找子代最差和最优个体
		std::tie(group_state.worstIndex, group_state.bestIndex) = this->selectPolarIndivs(0, 1);

		// 记录最优解
		bestIndivs.push_back(indivs[groupSize]);
	}

	// 进化迭代(单线程模式)
	template<class R, class... Args>
	void GAGroup<R(Args...)>::run()
	{
		while (true)
		{
			// 更新fitArrayCache数组
			updateFitArrayCache();

			crossover();
			mutate();
			select();

			// 迭代次数加一
			group_state.nGene++;
			
			// 记录当前时间
			group_state.nowTime = std::chrono::steady_clock::now();

			// 判断最优解相较于上一次的波动值
			if (abs(bestIndivs[group_state.nGene - 1].fitness - bestIndivs[group_state.nGene].fitness) <= group_state.stopTol)
			{
				group_state.count++;
			}
			else
			{
				group_state.count = 0;
			}

			// 判断是否到达停止条件
			if (getStopFlag())
			{
				return;
			}
		}
	}

	// 更新轮盘赌刻度线
	template<class R, class... Args>
	void GAGroup<R(Args...)>::updateFitArrayCache()
	{
		// 更新fitArrayCache数组
        // fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		for (int i = 1; i < groupSize; i++)
		{
			fitArrayCache[i] = fitArrayCache[i - 1] + (indivs[i - 1].fitness - indivs[group_state.worstIndex].fitness);
		}
	}

	// 寻找最差和最好的个体位置, 形参为线程ID(0 ~ N)及并行线程数量, 返回pair<worst, best>
	template<class R, class... Args>
	std::pair<int, int> GAGroup<R(Args...)>::selectPolarIndivs(const int seq, const int interval)
	{
		int worst = seq;
		int best = seq;

		// 遍历indivs的fitness，寻找fitness的极值，记录极值位置
		for (int i = seq + interval; i < groupSize; i += interval)
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
	bool GAGroup<R(Args...)>::getStopFlag()
	{
		// 是否达到最大迭代次数
		if (group_state.setMaxGeneFlag && group_state.nGene >= group_state.maxGene)
		{
			group_state.stopCode = 1;
		}

		// 是否达到最大迭代时间
		std::chrono::duration<double> evolTime = group_state.nowTime - group_state.startTime;
		if (group_state.setRuntimeFlag && evolTime.count() >= group_state.maxRuntime.value)
		{
			group_state.stopCode = 2;
		}

		// 是否连续5次最优解的波动小于停止误差
		if (group_state.setStopTolFlag && group_state.count == 5)
		{
			group_state.stopCode = 0;
		}

		// 根据stop code判断种群是否停止迭代
		if (group_state.stopCode != -1)
		{
			group_state.stopFlag = true;
		}

		return group_state.stopFlag;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	template<class R, class... Args>
	void GAGroup<R(Args...)>::cross_thread(const int seq)
	{
		int Index_M = 0;
		int Index_F = 0;
		double rand_cross = 0;

		while (true)
		{
			// 临界区
			{
				std::unique_lock<std::mutex> lck(thread_sync.mtx);
				thread_sync.cv.wait(lck, [this, &seq]() {return !(this->thread_sync).cross_flag[seq] && (this->thread_sync).crossReady || this->getStopFlag(); });
			}// 离开作用域, gm.thread_state.mtx的unlock()方法自动执行
			
			 // 是否停止迭代
			if (getStopFlag())
			{
				return;
			}

			// 父代个体交叉产生子代
			for (int i = seq; i < groupSize; i += 2 * thread_sync.threadNum)
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
			}

			// 线程同步
			thread_sync.mtx.lock();
			thread_sync.cross_sync(seq); // Change this function
			thread_sync.mtx.unlock();

			thread_sync.cv.notify_all();
		}
	}

	////////////////////////////////////////////////////////
	template<class R, class... Args>
	void GAGroup<R(Args...)>::mut_thread(const int seq)
	{
		double rand_num = 0; // 随机数

		while (true)
		{
			// 临界区
			{
				std::unique_lock<std::mutex> lck(thread_sync.mtx);
				thread_sync.cv.wait(lck, [this, &seq]() {return !(this->thread_sync).mut_flag[seq] && (this->thread_sync).mutReady || this->getStopFlag(); });
			}

			// 是否停止迭代
			if (getStopFlag())
			{
				return;
			}

			// 每个个体变异(最优解个体不发生变异)
			for (int i = seq; i < groupSize; i += thread_sync.threadNum)
			{
				if (i != group_state.bestIndex)
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
					// 计算个体适应度
					indivs[i].fitness = callFitFunc(indivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
				}
			}

			// 线程同步
			thread_sync.mtx.lock();
			thread_sync.mutate_sync(seq);
			thread_sync.mtx.unlock();

			thread_sync.cv.notify_all();
		}
	}

	// 选择
	template<class R, class ...Args>
	void GAGroup<R(Args...)>::sel_thread(const int seq)
	{
		int worst_temp = 0;
		int best_temp = 0;

		while (true)
		{
			// 临界区
			{
				std::unique_lock<std::mutex> lck(thread_sync.mtx);				
				thread_sync.cv.wait(lck, [this, &seq]() {return !(this->thread_sync).sel_flag[seq] && (this->thread_sync).selectReady || this->getStopFlag(); });
			}

			// 寻找子代最差和最优个体
			std::tie(worst_temp, best_temp) = this->selectPolarIndivs(seq, thread_sync.threadNum);

			// 线程同步
			thread_sync.mtx.lock();			
			//更新最优与最差个体位置;
			if (indivs[worst_temp].fitness < indivs[group_state.worstIndex].fitness)
			{
				group_state.worstIndex = worst_temp;
			}
			if (indivs[best_temp].fitness > indivs[group_state.bestIndex].fitness)
			{
				group_state.bestIndex = best_temp;
			}

			thread_sync.select_sync(seq);
			thread_sync.mtx.unlock();

			thread_sync.cv.notify_all();

			// 是否停止迭代
			if (getStopFlag())
			{
				return;
			}
		}
	}
}

#endif