#ifndef _PSO_H_
#define _PSO_H_

#include <chrono>
#include <vector>
#include "opt_time.h"
#include "pso_individual.h"
#include "pso_state.h"
#include "index_seq.h"
#include "pso_info.h"
#include "range_random.h"
#include "pso_thread_sync.h"

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
        PSO_Individual* indivs;                                               // 粒子群数组
        R(*fitFunc)(Args...);                                                 // 适应度函数指针
        std::function<double(std::size_t)> reweight;                          // 种群数量动态调整
        std::function<void(const PSO_Info&)> monitor;                         // 外部监听器
        double(*bound)[2];                                                    // 每个变量的区间, 以数组指针表示
        std::vector<PSO_Individual> bestIndivs;                               // 记录每次迭代的最优个体
        PSO_Individual best_indiv;                                            // 最优粒子
        PSO_State group_state;                                                // 粒子群状态

        std::vector<std::thread> vec_run;                                     // 交叉线程组
        std::thread single_thread;                                            // 单线程

        std::unique_ptr< PSOThreadSync<R, Args...> > thread_sync;             // 线程同步器

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
        void setThreadNum(const std::size_t NUM);                             // 设置并行计算的线程数，默认为1
        void setMonitor(const std::function<void(const PSO_Info&)>&);         // 设置外部监听函数
        void setReweight(const std::function<std::size_t(std::size_t)>&);     // 种群数量调整函数

        int getNVars()const;                                                  // 获取种群变量个数
        int getGeneration()const;                                             // 获得当前种群代数
        int getGroupSize()const;                                              // 获得当前种群个体数量
        std::vector<PSO_Individual> getBestIndivs()const;                     // 获取历次迭代的最优解
        int getStopCode();                                                    // 获取Stop Code

        void initGroup(const std::vector<PSO_Individual>& indivs = std::vector<PSO_Individual>());    // 初始化粒子群
        bool start();                                                         // 开始迭代优化
        void pause();                                                         // 暂停(为保证数据一致性，需在一次完整迭代后pause)
        void proceed();                                                       // 继续迭代
        void kill();                                                          // 结束迭代
        void wait_result();                                                   // 阻塞当前线程,等待优化结果
        PSO<R(Args...)> clone();                                              // 克隆当前粒子群

    private:
        template<std::size_t... I>
        R callFitFunc(double* args, const opt::index_seq<I...>&);              // 调用适应度函数数

        std::size_t findBest();                                                // 寻找最优个体
        void iter();                                                           // 进行一次迭代
        void updateStopState();                                                // 更新粒子群停止状态
        void run();
        bool flushStopFlag();                                                  // 判断是否结束迭代
    };

    /******************************************* 构造与析构 ***********************************************************/
    // 构造函数，需提供供适应度函数及粒子群数量
    template<class R, class... Args>
    PSO<R(Args...)>::PSO(R(*f)(Args...), const std::size_t size)
        :groupSize(size),
        nVars(sizeof...(Args)),
        indivs(new PSO_Individual[groupSize]),
        fitFunc(f),
        reweight([](std::size_t n) {return 0.5; }),
        bound(nullptr),
        thread_sync(new PSOThreadSync<R, Args...>(this))
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
        bestIndivs(),
        thread_sync(nullptr)
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
        thread_sync.reset(new PSOThreadSync<R, Args...>(*(other.thread_sync), this));
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
        bestIndivs(),
        thread_sync(nullptr)
    {
        other.kill();        // 终止种群迭代

        other.indivs = nullptr;
        other.bound = nullptr;

        this->bestIndivs = std::move(other.bestIndivs);

        thread_sync.reset(new PSOThreadSync<R, Args...>(*(other.thread_sync), this));
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

    // 设置线程数量
    template<class R, class ...Args>
    void PSO<R(Args...)>::setThreadNum(const std::size_t NUM)
    {
        if (NUM > 1)
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
    std::vector<PSO_Individual> PSO<R(Args...)>::getBestIndivs()const
    {
        return bestIndivs;
    }

    // 获取停止代码
    template<class R, class... Args>
    int PSO<R(Args...)>::getStopCode()
    {
        return group_state.stopCode;
    }

    // 初始化粒子群
    template<class R, class... Args>
    void PSO<R(Args...)>::initGroup(const std::vector<PSO_Individual>& init_indivs)
    {
        // 初始化前需要保证已设置基因变量区间
        if (group_state.setBoundFlag)
        {
            // 1.初始化每一个粒子
            for (std::size_t i = 0; i < init_indivs.size(); i++)
            {
                indivs[i] = init_indivs[i];
                indivs[i].fitness = callFitFunc(indivs[i].xs, opt::make_index_seq<sizeof...(Args)>());
            }
            for (std::size_t i = init_indivs.size(); i < groupSize; i++)
            {
                // 设置个体位置
                for (std::size_t j = 0; j < nVars; j++)
                {
                    indivs[i].xs[j] = random_real(bound[j][0], bound[j][1]);
                    indivs[i].best_xs[j] = indivs[i].xs[j];
                }
                // 个体适应度
                indivs[i].fitness = callFitFunc(indivs[i].xs, opt::make_index_seq<sizeof...(Args)>());
            }

            // 2.寻找最优个体
            std::size_t bestIndex = this->findBest();
            bestIndivs.push_back(indivs[bestIndex]);                     // 记录最优解
            best_indiv = indivs[bestIndex];                              // 最优个体

            // 3.设置粒子飞行速度
            double c2 = 2;   // 学习因子

            // 计算各个粒子飞行速度，初始化粒子速度没有惯性部分和自身部分，只有社会部分
            for (std::size_t i = 0; i < groupSize; i++)
            {
                // 随机数
                double r2 = opt::random_real(0, 1);
                for (std::size_t j = 0; j < nVars; j++)
                {
                    indivs[i].vs[j] = c2 * r2 * (indivs[bestIndex].xs[j] - indivs[i].xs[j]);

                    // 最大飞行速度
                    double V_max = 0.15 * (bound[j][1] - bound[j][0]);

                    // 检查速度是否过大
                    if (indivs[i].vs[j] > V_max)
                    {
                        indivs[i].vs[j] = V_max;
                    }
                    if (indivs[i].vs[j] < -V_max)
                    {
                        indivs[i].vs[j] = -V_max;
                    }
                }
            }

            // 4.设置初始化标志位
            group_state.initFlag = true;
        }
        else
        {
            throw std::string("Fitness function or variables boundary is not set.");
        }
    }

    // 开始迭代优化
    template<class R, class ...Args>
    bool PSO<R(Args...)>::start()
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
                    //vec_run.emplace_back(&GAGroup<R(Args...)>::run_parallel, this, i);
                }
            }
            else // 单线程模式
            {
                single_thread = std::thread(&PSO<R(Args...)>::run, this);
            }
            return true;
        }
        return false;
    }

    // 暂停(为保证数据一致性，需在一次完整迭代后pause)
    template<class R, class ...Args>
    void PSO<R(Args...)>::pause()
    {
        // 如果粒子群未处于迭代运行状态
        if (group_state.stopCode != 0)
        {
            return;
        }

        this->group_state.sleep.signal = true;
        std::unique_lock<std::mutex> lck(thread_sync->mtx);
        thread_sync->cv.wait(lck, [this]() {return this->group_state.sleep.result; });
    }

    // 继续迭代
    template<class R, class ...Args>
    void PSO<R(Args...)>::proceed()
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

    // 结束迭代
    template<class R, class ...Args>
    void PSO<R(Args...)>::kill()
    {
        // 如果种群未处于迭代运行状态
        if (group_state.stopCode != 0)
        {
            return;
        }

        this->group_state.stopCode = 4;
        this->wait_result();
    }

    // 阻塞当前线程,等待优化结果
    template<class R, class ...Args>
    void PSO<R(Args...)>::wait_result()
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

    // 克隆当前粒子群
    template<class R, class ...Args>
    PSO<R(Args...)> PSO<R(Args...)>::clone()
    {
        return PSO<R(Args...)>(*this);
    }

    // 寻找最优个体
    template<class R, class ...Args>
    std::size_t PSO<R(Args...)>::findBest()
    {
        std::size_t index = 0;
        for (std::size_t i = 1; i < groupSize; i++)
        {
            if (indivs[i].fitness > indivs[index].fitness)
            {
                index = i;
            }
        }
        return index;
    }

    // 进行一次迭代
    template<class R, class ...Args>
    void PSO<R(Args...)>::iter()
    {
        double weight = reweight(group_state.nGene + 1);         // 惯性系数

        double c1 = 1.49618;                                     // 学习因子
        double c2 = 1.49618;

        double r1 = 0.5;
        double r2 = 0.5;

        double V_max = 0;                                        // 最大飞行速度

        double temp = 0;
        double temp_1 = 0;
        double temp_2 = 0;
        double temp_3 = 0;

        // 更新速度
        for (std::size_t j = 0; j < nVars; j++)
        {
            V_max = 0.15 * (bound[j][1] - bound[j][0]);

            for (std::size_t i = 0; i < groupSize; i++)
            {
                r1 = opt::random_real(0, 1);
                r2 = opt::random_real(0, 1);

                temp_1 = weight * indivs[i].vs[j];                                // 惯性部分
                temp_2 = c1 * r1 * (indivs[i].best_xs[j] - indivs[i].xs[j]);      // 自我部分
                temp_3 = c2 * r2 * (best_indiv.xs[j] - indivs[i].xs[j]);          // 社会部分

                temp = temp_1 + temp_2 + temp_3;                                  // 飞行速度

                if (temp > V_max)
                {
                    temp = V_max;
                }
                if (temp < -V_max)
                {
                    temp = -V_max;
                }

                // 更新位置
                indivs[i].xs[j] = indivs[i].xs[j] + indivs[i].vs[j];
                if (indivs[i].xs[j] < bound[j][0])
                {
                    indivs[i].xs[j] = bound[j][0];
                }
                if (indivs[i].xs[j] > bound[j][1])
                {
                    indivs[i].xs[j] = bound[j][1];
                }
                // 更新速度
                indivs[i].vs[j] = temp;
            }
        }

        // 更新粒子适应度
        for (std::size_t i = 0; i < groupSize; i++)
        {
            indivs[i].fitness = callFitFunc(indivs[i].xs, opt::make_index_seq<sizeof...(Args)>());
            double best_fit = callFitFunc(indivs[i].best_xs, opt::make_index_seq<sizeof...(Args)>());
            if (indivs[i].fitness > best_fit)
            {
                for (std::size_t j = 0; j < nVars; j++)
                {
                    indivs[i].best_xs[j] = indivs[i].xs[j];
                }
            }
        }

        std::size_t bestIndex = this->findBest();                    // 寻找当前粒子群最优个体
        bestIndivs.push_back(indivs[bestIndex]);                     // 记录最优解

        if (indivs[bestIndex].fitness > best_indiv.fitness)          // 最优个体体
        {
            best_indiv = indivs[bestIndex];
        }
    }

    // 更新粒子群停止状态
    template<class R, class ...Args>
    void PSO<R(Args...)>::updateStopState()
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

    // 进化迭代(单线程模式)
    template<class R, class ...Args>
    void PSO<R(Args...)>::run()
    {
        while (true)
        {
            iter();                               // 进行一次迭代
            updateStopState();                    // 更新停止状态

            // 调用外部监听函数
            if (monitor)
            {
                PSO_Info pso_info(group_state.time, group_state.nGene, best_indiv);
                monitor(pso_info);
            }

            // 判断是否到达停止条件
            if (flushStopFlag())
            {
                return;
            }

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

    // 判断是否结束迭代
    template<class R, class ...Args>
    bool PSO<R(Args...)>::flushStopFlag()
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

    // 调用适应度函数
    template<class R, class... Args>
    template<std::size_t... I>
    inline R PSO<R(Args...)>::callFitFunc(double* args, const opt::index_seq<I...>&)
    {
        return fitFunc(args[I]...);
    }
}

#endif
