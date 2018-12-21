#ifndef _Genetic_Algorithm_
#define _Genetic_Algorithm_

#include <string>
#include <chrono>
#include <tuple>
#include <vector>
#include <cmath>
#include <thread>
#include <functional>

#include "ga_individual.h"
#include "range_random.h"
#include "roulette.h"
#include "opt_time.h"
#include "ga_state.h"
#include "cross_factor.h"
#include "mutate_factor.h"
#include "ga_thread_sync.h"
#include "index_seq.h"
#include "ga_info.h"

////////
#include <iostream>

namespace opt
{
    template<class F> class GAGroup;

    // GA种群类，需提供适应度函数类型
    template<class R, class... Args>
    class GAGroup<R(Args...)>
    {
        using GenBound = std::initializer_list<std::initializer_list<double>>;
        friend class GAThreadSync<R, Args...>;

    private:
        std::string name;                                                     // 种群名称
        std::size_t groupSize;                                                // 初始种群个体数量，default = 1000；
        std::size_t groupCapacity;                                            //
        const std::size_t nVars;                                              // 适应度函数包含的变量个数
        GA_Individual* indivs;                                                // 种群个体, 指向groupSize个个体数组(最后一个存放最优个体)
        GA_Individual* tempIndivs;                                            // 子代个体缓存区
        R(*fitFunc)(Args...);                                                 // 适应度函数指针
        std::function<void(const GA_Info&)> monitor;                          // 外部监听器
        std::function<std::size_t(std::size_t)> resize;                       // 种群数量动态调整
        double(*bound)[2];                                                    // 每个变量(基因)的区间, 以数组指针表示
        Roulette<double> roulette;                                            // 轮盘赌对象
        double mutateProb;                                                    // 个体变异概率,默认p = 0.1
        double crossProb;                                                     // 个体交叉概率, 默认p = 0.6
        std::vector<GA_Individual> bestIndivs;                                // 记录每次迭代的最优个体

        std::vector<std::thread> vec_run;                                     // 交叉线程组
        std::thread single_thread;                                            // 单线程

        GA_State group_state;                                                 // 种群状态
        std::unique_ptr< GAThreadSync<R, Args...> > thread_sync;              // 线程同步器

    public:
        GAGroup(R(*f)(Args...), const std::size_t size = 1000);               // 构造函数，构造一个种群，需提供适应度函数及种群数量
        GAGroup(GAGroup<R(Args...)>& other);                                  // 拷贝构造
        GAGroup(GAGroup<R(Args...)>&& other);                                 // 移动构造
        GAGroup<R(Args...)>& operator=(const GAGroup<R(Args...)>& other) = delete;
        GAGroup<R(Args...)>& operator=(GAGroup<R(Args...)>&& other) = delete;
        ~GAGroup();

        void setName(const std::string& str);                                 // 设置种群名称
        void setBoundary(double(*b)[2]);                                      // 设置变量区间
        void setBoundary(const GenBound& b);                                  // 设置变量区间
        void setMaxGeneration(const unsigned int N);                          // 设置最大迭代次数
        void setMaxRuntime(const Second& time);                               // 设置最大运行时间（秒）
        void setStopTol(const long double t, const unsigned int N = 5);       // 设置最优解停止误差
        void setMutateProb(const double p);                                   // 设置基因变异概率
        void setCrossProb(const double p);                                    // 设置交叉概率
        void setThreadNum(const std::size_t NUM);                             // 设置并行计算的线程数，默认为1
        void setMonitor(const std::function<void(const GA_Info&)>&);          // 设置外部监听函数
        void setResize(const std::function<std::size_t(std::size_t)>&);       // 种群数量调整函数

        const std::string getName()const;                                     // 获取种群名称
        int getNVars()const;                                                  // 获取种群变量个数
        int getGeneration()const;                                             // 获得当前种群代数
        int getGroupSize()const;                                              // 获得当前种群个体数量
        std::vector<GA_Individual> getBestIndivs()const;                      // 获取历次迭代的最优解
        int getStopCode();                                                    // 获取Stop Code

        void initGroup(const std::vector<GA_Individual>& indivs = std::vector<GA_Individual>());    // 初始化种群个体
        bool start();                                                         // 开始迭代进化
        void pause();                                                         // 暂停(为保证数据一致性，需在一次完整迭代后pause)
        void proceed();                                                       // 继续迭代
        void kill();                                                          // 结束迭代
        void wait_result();                                                   // 阻塞当前线程,等待优化结果
        GAGroup<R(Args...)> clone();                                          // 克隆当前种群

    private:
        void crossover(const std::size_t seq);                                // 交叉(单线程模式)
        void mutate(const std::size_t seq);                                   // 变异(单线程模式)
        void select(const std::size_t seq);                                   // 选择(单线程模式)

        void run();                                                           // 单线程迭代
        void run_parallel(const std::size_t seq);                             // 并行迭代

        void updateRoulette(const std::size_t size);                          // 更新轮盘赌刻度线
        void updateStopState();                                               // 更新种群停止状态
        void switchIndivArray();

        std::pair<int, int> selectPolarIndivs(const int seq, const int interval);   // 寻找最差和最好的个体位置,返回<worst, best>
        template<std::size_t... I>
        R callFitFunc(double* args, const opt::index_seq<I...>&);              // 调用适应度函数
        bool flushStopFlag();                                                  // 判断是否结束迭代
    };

    /******************************************* 构造与析构 ***********************************************************/
    // 构造函数，需提供供适应度函数及种群数量
    template<class R, class... Args>
    GAGroup<R(Args...)>::GAGroup(R(*f)(Args...), const std::size_t size)
        :name("None"),
        groupSize(size + size % 2),
        groupCapacity(groupSize),
        nVars(sizeof...(Args)),
        indivs(new GA_Individual[groupSize]),
        tempIndivs(new GA_Individual[groupSize]),
        fitFunc(f),
        bound(nullptr),
        roulette(groupSize + 1),
        mutateProb(0.1),
        crossProb(0.6),
        group_state(),
        thread_sync(new GAThreadSync<R, Args...>(this))
    {
        for (std::size_t i = 0; i < groupSize; i++)
        {
            indivs[i] = GA_Individual(nVars);
            tempIndivs[i] = GA_Individual(nVars);
        }

        // 变量区间
        bound = new double[nVars][2];
    }

    // 复制构造
    template<class R, class... Args>
    GAGroup<R(Args...)>::GAGroup(GAGroup<R(Args...)>& other)
        :name(other.name),
        groupSize(other.groupSize),
        groupCapacity(groupSize),
        nVars(other.nVars),
        indivs(new GA_Individual[groupSize]),
        tempIndivs(new GA_Individual[groupSize]),
        fitFunc(other.fitFunc),
        monitor(other.monitor),
        resize(other.resize),
        bound(nullptr),
        roulette(),
        mutateProb(other.mutateProb),
        crossProb(other.crossProb),
        bestIndivs(),
        group_state(),
        thread_sync(nullptr)
    {
        other.pause();

        for (std::size_t i = 0; i < groupSize; i++)
        {
            indivs[i] = (other.indivs)[i];
            tempIndivs[i] = (other.indivs)[i];
        }

        // 变量区间
        bound = new double[nVars][2];
        for (std::size_t i = 0; i < nVars; i++)
        {
            bound[i][0] = (other.bound)[i][0];
            bound[i][1] = (other.bound)[i][1];
        }

        this->roulette = other.roulette;
        this->bestIndivs = other.bestIndivs;

        this->group_state = other.group_state;
        this->group_state.sleep.signal = false;
        this->group_state.sleep.result = false;
        this->group_state.stopCode = -1;   // 未开始迭代

        this->thread_sync.reset(new GAThreadSync<R, Args...>(*(other.thread_sync), this));

        other.proceed();
    }

    // 移动构造
    template<class R, class... Args>
    GAGroup<R(Args...)>::GAGroup(GAGroup<R(Args...)>&& other)
        :name(other.name),
        groupSize(other.groupSize),
        groupCapacity(other.groupCapacity),
        nVars(other.nVars),
        indivs(other.indivs),
        tempIndivs(other.tempIndivs),
        fitFunc(other.fitFunc),
        monitor(other.monitor),
        resize(other.resize),
        bound(other.bound),
        roulette(),
        mutateProb(other.mutateProb),
        crossProb(other.crossProb),
        bestIndivs(),
        group_state(other.group_state),
        thread_sync(nullptr)
    {
        other.kill();        // 终止种群迭代

        other.indivs = nullptr;
        other.tempIndivs = nullptr;
        other.bound = nullptr;

        this->roulette = std::move(other.roulette);
        this->bestIndivs = std::move(other.bestIndivs);

        this->group_state = other.group_state;
        this->group_state.sleep.signal = false;
        this->group_state.sleep.result = false;
        this->group_state.stopCode = -1;

        thread_sync.reset(new GAThreadSync<R, Args...>(*(other.thread_sync), this));
    }

    // 析构函数
    template<class R, class... Args>
    GAGroup<R(Args...)>::~GAGroup()
    {
        delete[] indivs;
        delete[] tempIndivs;
        delete[] bound;
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
    void GAGroup<R(Args...)>::setMaxGeneration(const unsigned int N)
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
    void GAGroup<R(Args...)>::setStopTol(const long double t, const unsigned int N)
    {
        group_state.converCount = N;
        group_state.setStopTolFlag = true;
        group_state.stopTol = t;
    }

    // 设置基因变异概率
    template<class R, class... Args>
    void GAGroup<R(Args...)>::setMutateProb(const double p)
    {
        mutateProb = p;
    }

    // 设置基因交叉概率
    template<class R, class... Args>
    void GAGroup<R(Args...)>::setCrossProb(const double p)
    {
        crossProb = p;
    }

    // 设置并行计算的线程数，默认为1
    template<class R, class... Args>
    void GAGroup<R(Args...)>::setThreadNum(const std::size_t NUM)
    {
        if (NUM >= 1)
        {
            thread_sync->setThreadNum(NUM);
        }
        if (NUM < 1)
        {
            throw std::string("Thread number error.");
        }
    }

    // 设置外部监听函数
    template<class R, class... Args>
    void GAGroup<R(Args...)>::setMonitor(const std::function<void(const GA_Info&)>& func)
    {
        this->monitor = func;
    }

    // 种群数量调整函数
    template<class R, class... Args>
    void GAGroup<R(Args...)>::setResize(const std::function<std::size_t(std::size_t)>& func)
    {
        this->resize = func;
    }

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

    // 获取历次迭代的最优解
    template<class R, class... Args>
    std::vector<GA_Individual> GAGroup<R(Args...)>::getBestIndivs()const
    {
        return bestIndivs;
    }

    // 获取停止代码
    template<class R, class... Args>
    int GAGroup<R(Args...)>::getStopCode()
    {
        return group_state.stopCode;
    }

    /******************************************************************************************************************/
    // 开始迭代
    template<class R, class... Args>
    bool GAGroup<R(Args...)>::start()
    {
        // 如果种群已被初始化则不再重复初始化
        if (group_state.initFlag == false)
        {
            initGroup();
        }

        // 若种群满足进化条件
        if (group_state.runable())
        {
            group_state.startTime = std::chrono::steady_clock::now();
            group_state.stopCode = 0; // 正在迭代

            if (thread_sync->threadNum > 1) // 多线程模式
            {
                // 构造种群交叉线程组
                for (std::size_t i = 0; i < thread_sync->threadNum; i++)
                {
                    vec_run.emplace_back(&GAGroup<R(Args...)>::run_parallel, this, i);
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

    // 进化迭代(单线程模式)
    template<class R, class... Args>
    void GAGroup<R(Args...)>::run()
    {
        while (true)
        {
            // 更新更新轮盘赌刻度线
            updateRoulette(groupSize);

            // 更新子代个体数量
            if (resize)
            {
                groupSize = resize(group_state.nGene + 1);
                if (groupSize < 0) { groupSize = 0; }
                groupSize = groupSize + groupSize % 2;
            }
            crossover(0);                     // 交叉
            switchIndivArray();               // 交换个体数组
            mutate(0);                        // 变异
            select(0);                        // 环境选择
            updateStopState();                // 更新停止状态

            // 调用外部监听函数
            if (monitor)
            {
                GA_Info ga_info(group_state.time, group_state.nGene, bestIndivs.back());
                monitor(ga_info);
            }

            // 判断是否到达停止条件
            if (flushStopFlag()) { return; }

            // 是否暂停迭代
            if (group_state.sleep.signal == true)
            {
                group_state.sleep.result = true;
                thread_sync->cv.notify_all();
            }

            // 是否暂停迭代(为保证数据一致性，需在一次完整迭代后pause)
            {
                std::unique_lock<std::mutex> lck(thread_sync->mtx);
                thread_sync->cv.wait(lck, [this]() { return !(this->group_state.sleep.signal); });
            }
        }
    }

    // 进化迭代(多线程模式)
    template<class R, class ...Args>
    void GAGroup<R(Args...)>::run_parallel(const std::size_t seq)
    {
        while (true)
        {
            /////////////////////////////// Crossover /////////////////////////////////////////////////////
            crossover(seq);

            // 线程同步
            thread_sync->mtx.lock();
            thread_sync->cross_sync(seq);
            thread_sync->mtx.unlock();

            thread_sync->cv.notify_all();

            {
                std::unique_lock<std::mutex> lck(thread_sync->mtx);
                thread_sync->cv.wait(lck, [this, &seq]() {
                    return !(this->thread_sync->mut_flag[seq]) && this->thread_sync->mutReady;
                    });
            }

            /////////////////////////////// Mutate /////////////////////////////////////////////////////
            mutate(seq);

            // 线程同步
            thread_sync->mtx.lock();
            thread_sync->mutate_sync(seq);
            thread_sync->mtx.unlock();
            thread_sync->cv.notify_all();

            {
                std::unique_lock<std::mutex> lck(thread_sync->mtx);
                thread_sync->cv.wait(lck, [this, &seq]() {
                    return !(this->thread_sync->sel_flag[seq]) && this->thread_sync->selectReady;
                    });
            }

            /////////////////////////////// Select /////////////////////////////////////////////////////
            int worst_temp = 0;
            int best_temp = 0;

            int thread_num = thread_sync->threadNum;

            // 寻找子代最差和最优个体
            std::tie(worst_temp, best_temp) = this->selectPolarIndivs(seq, thread_num);

            // 线程同步
            thread_sync->mtx.lock();

            // 更新最优与最差个体位置;
            if (indivs[worst_temp].fitness < indivs[group_state.worstIndex].fitness)
            {
                group_state.worstIndex = worst_temp;
            }
            if (indivs[best_temp].fitness > indivs[group_state.bestIndex].fitness)
            {
                group_state.bestIndex = best_temp;
            }

            thread_sync->select_sync(seq);
            thread_sync->mtx.unlock();
            thread_sync->cv.notify_all();

            // 为保证数据一致性，需在一次完整迭代后pause
            {
                std::unique_lock<std::mutex> lck(thread_sync->mtx);
                thread_sync->cv.wait(lck, [this, &seq]() {
                    return !(this->thread_sync->cross_flag[seq]) && this->thread_sync->crossReady && !(this->group_state.sleep.signal);
                    });
            }

            if (group_state.stopFlag == true) { return; }
        }
    }

    // 阻塞当前线程, 等待计算结果
    template<class R, class... Args>
    void GAGroup<R(Args...)>::wait_result()
    {
        if (thread_sync->threadNum == 1) // 单线程模式
        {
            if (single_thread.joinable())
            {
                single_thread.join();
            }
        }
        else  // 多线程模式
        {
            for (std::size_t i = 0; i < thread_sync->threadNum; i++)
            {
                if (vec_run[i].joinable())
                {
                    vec_run[i].join();
                }
            }
        }
    }

    // 克隆当前状态种群
    template<class R, class... Args>
    GAGroup<R(Args...)> GAGroup<R(Args...)>::clone()
    {
        return GAGroup<R(Args...)>(*this);
    }

    // 暂停(为保证数据一致性，需在一次完整迭代后pause)
    template<class R, class... Args>
    void GAGroup<R(Args...)>::pause()
    {
        // 如果种群未处于迭代运行状态
        if (group_state.stopCode != 0)
        {
            return;
        }

        this->group_state.sleep.signal = true;
        std::unique_lock<std::mutex> lck(thread_sync->mtx);
        thread_sync->cv.wait(lck, [this]() {return this->group_state.sleep.result; });
    }

    // 继续迭代
    template<class R, class... Args>
    void GAGroup<R(Args...)>::proceed()
    {
        // 如果种群未处于迭代运行状态
        if (group_state.stopCode != 0)
        {
            return;
        }

        this->group_state.sleep.signal = false;
        this->group_state.sleep.result = false;
        this->thread_sync->cv.notify_all();
    }

    // 停止迭代
    template<class R, class... Args>
    void GAGroup<R(Args...)>::kill()
    {
        // 如果种群未处于迭代运行状态
        if (group_state.stopCode != 0)
        {
            return;
        }

        this->group_state.stopCode = 4;
        this->wait_result();
    }

    /******************************************* Private Functions *****************************************************/
    // 初始化种群
    template<class R, class... Args>
    void GAGroup<R(Args...)>::initGroup(const std::vector<GA_Individual>& init_indivs)
    {
        // 初始化前需要保证已设置基因变量区间
        if (group_state.setBoundFlag)
        {
            // 1.初始化每一个个体
            for (std::size_t i = 0; i < init_indivs.size(); i++)
            {
                indivs[i] = init_indivs[i];
                indivs[i].fitness = callFitFunc(indivs[i].vars, opt::make_index_seq<sizeof...(Args)>());
            }
            for (std::size_t i = init_indivs.size(); i < groupSize; i++)
            {
                // 设置个体初始基因
                for (std::size_t j = 0; j < nVars; j++)
                {
                    indivs[i].vars[j] = random_real(bound[j][0], bound[j][1]);
                }
                // 个体适应度
                indivs[i].fitness = callFitFunc(indivs[i].vars, opt::make_index_seq<sizeof...(Args)>());
            }

            // 2.寻找最差和最好适应度个体
            std::tie(group_state.worstIndex, group_state.bestIndex) = this->selectPolarIndivs(0, 1);
            bestIndivs.push_back(indivs[group_state.bestIndex]);                     // 记录最优解

            // 3.更新轮盘赌刻度线
            updateRoulette(groupSize);

            // 4.设置初始化标志位
            group_state.initFlag = true;

            // 5.初始化线程同步器
            thread_sync->crossReady = true;
            thread_sync->mutReady = false;
            thread_sync->selectReady = false;
        }
        else
        {
            throw std::string("Fitness function or genes boundary is not set.");
        }
    }

    // 交叉(单线程模式)
    template<class R, class... Args>
    void GAGroup<R(Args...)>::crossover(const std::size_t seq)
    {
        int Index_M = 0;             // 父个体
        int Index_F = 0;             // 母个体
        double rand_cross = 0;       // 随机数缓存

        int thread_num = thread_sync->threadNum;          // 线程数

        // 随机选取个体交叉, fitness越大，被选中的概率越大
        for (std::size_t i = seq * 2; i < groupSize; i += 2 * thread_num)
        {
            // 随机选取两个个体作为父代
            Index_M = roulette.roll();
            Index_F = roulette.roll();

            // 生成子代个体, 存放于缓存数组tempIndivs中
            for (std::size_t j = 0; j < nVars; j++)
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
    }

    // 变异(单线程模式)
    template<class R, class... Args>
    void GAGroup<R(Args...)>::mutate(const std::size_t seq)
    {
        double rand_num = 0; // 随机数

        int thread_num = thread_sync->threadNum;  // 线程数

        // 每个个体的变异
        for (std::size_t i = seq; i < groupSize; i += thread_num)
        {
            // 每个基因的变异
            for (std::size_t j = 0; j < nVars; j++)
            {
                rand_num = random_real(0, 1);
                if (rand_num < mutateProb)
                {
                    // 单个基因变异
                    indivs[i].vars[j] = mutate_PM(indivs[i].vars[j], bound[j][0], bound[j][1]);
                }
            }

            // 计算每个个体适应度
            indivs[i].fitness = callFitFunc(indivs[i].vars, opt::make_index_seq<sizeof...(Args)>());
        }
    }

    // 选择(单线程模式)
    template<class R, class... Args>
    void GAGroup<R(Args...)>::select(const std::size_t seq)
    {
        int thread_num = thread_sync->threadNum; // 线程数

        // 寻找子代最差和最优个体
        std::tie(group_state.worstIndex, group_state.bestIndex) = this->selectPolarIndivs(seq, thread_num);

        // 淘汰子代最差个体, 用父代最优个体取代
        indivs[group_state.worstIndex] = bestIndivs.back();

        // 如果子代出现更优个体
        if (indivs[group_state.bestIndex].fitness >= bestIndivs.back().fitness)
        {
            bestIndivs.push_back(indivs[group_state.bestIndex]);
        }
        else
        {
            bestIndivs.push_back(bestIndivs.back());
        }
    }

    // 更新轮盘赌刻度线
    template<class R, class... Args>
    void GAGroup<R(Args...)>::updateRoulette(const std::size_t size)
    {
        int worst = 0;

        // 遍历indivs的fitness，寻找fitness的最小值,记录最小值位置
        for (std::size_t i = 1; i < size; i++)
        {
            if (indivs[i].fitness < indivs[worst].fitness)
            {
                worst = i;
            }
        }

        // 重新划分刻度线
        roulette.reset(size + 1);

        // 更新轮盘赌刻度线
        for (std::size_t i = 1; i < size + 1; i++)
        {
            roulette[i] = roulette[i - 1] + (indivs[i - 1].fitness - indivs[worst].fitness);
        }
    }

    // 更新种群停止状态
    template<class R, class... Args>
    void GAGroup<R(Args...)>::updateStopState()
    {
        // 迭代次数加一
        group_state.nGene++;

        // 记录当前时间
        group_state.nowTime = std::chrono::steady_clock::now();

        // 是否达到最大迭代时间
        std::chrono::duration<double> evolTime = group_state.nowTime - group_state.startTime;
        group_state.time = evolTime.count(); // 秒

        // 判断最优个体fitness值较上一代的波动情况
        if (abs(bestIndivs[group_state.nGene - 1].fitness - bestIndivs[group_state.nGene].fitness) <= group_state.stopTol)
        {
            group_state.count++;
        }
        else
        {
            group_state.count = 0;
        }
    }

    // 更新个体数组指针
    template<class R, class... Args>
    void GAGroup<R(Args...)>::switchIndivArray()
    {
        // 交换子代个体缓存区和父代个体存放区的指针
        GA_Individual* temp = indivs;
        indivs = tempIndivs;
        tempIndivs = temp;
    }

    // 寻找最差和最好的个体位置, 形参为线程ID(0 ~ N)及并行线程数量, 返回pair<worst, best>
    template<class R, class... Args>
    std::pair<int, int> GAGroup<R(Args...)>::selectPolarIndivs(const int seq, const int interval)
    {
        int worst = seq;
        int best = seq;

        // 遍历indivs的fitness，寻找fitness的极值，记录极值位置
        for (std::size_t i = seq + interval; i < groupSize; i += interval)
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

    // 调用适应度函数
    template<class R, class... Args>
    template<std::size_t... I>
    inline R GAGroup<R(Args...)>::callFitFunc(double* args, const opt::index_seq<I...>&)
    {
        return fitFunc(args[I]...);
    }

    // 判断是否结束迭代
    template<class R, class... Args>
    bool GAGroup<R(Args...)>::flushStopFlag()
    {
        // Stop Code : -1-未开始迭代; 0-正在迭代; 1-达到最大迭代次数; 2-达到最大迭代时间; 3-最优解收敛于稳定值; 4-人为停止迭代

        // 是否达到最大迭代次数
        if (group_state.setMaxGeneFlag && group_state.nGene >= group_state.maxGene)
        {
            group_state.stopCode = 1;
        }

        if (group_state.setRuntimeFlag && group_state.time >= group_state.maxRuntime.value)
        {
            group_state.stopCode = 2;
        }

        // 是否连续5次最优解的波动小于停止误差
        if (group_state.setStopTolFlag && group_state.count == group_state.converCount)
        {
            group_state.stopCode = 3;
        }

        // 根据stop code判断种群是否停止迭代
        if (group_state.stopCode != -1 && group_state.stopCode != 0)
        {
            group_state.stopFlag = true;
        }

        return group_state.stopFlag;
    }
}

#endif
