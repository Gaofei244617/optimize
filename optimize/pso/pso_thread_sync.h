#ifndef _PSO_THREAD_SYNC_
#define _PSO_THREAD_SYNC_

#include <thread>
#include <vector>
#include <utility>
#include <condition_variable>
#include <tuple>
#include "..\bool_array.h"
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
        std::size_t threadNum;                            // 并行计算使用的线程数
        bool_array iter_flag;                             // 并行计算: 线程状态标志
        PSO<R(Args...)>* pso;                             // 同步器关联的粒子群

    public:
        // 常规构造
        PSOThreadSync(PSO<R(Args...)>* ptr)
            :threadNum(1),
            iter_flag(threadNum),
            pso(ptr)
        {}

        // 复制构造
        PSOThreadSync(const PSOThreadSync& other, PSO<R(Args...)>* ptr)
            :threadNum(other.threadNum),
            iter_flag(other.iter_flag),
            pso(ptr)
        {}

        // 移动构造
        PSOThreadSync(PSOThreadSync&& other) = delete;

        // 设置线程数量
        void setThreadNum(const std::size_t N);

        // "迭代"线程同步
        void iter_sync(const std::size_t thread_id);
    };

    // 设置线程数量
    template<class R, class... Args>
    void PSOThreadSync<R, Args...>::setThreadNum(const std::size_t N)
    {
        threadNum = N;
        iter_flag.set_length(N);
    }

    // 迭代线程同步
    template<class R, class ...Args>
    void PSOThreadSync<R, Args...>::iter_sync(const std::size_t thread_id)
    {
        iter_flag[thread_id] = true;

        // 若迭代线程组全部运行完
        if (iter_flag.is_all_true())
        {
            //std::cout << "A" << thread_id << std::endl;
            // 更新全局最优解
            pso->updateGlobalBest();

            // 更新粒子群状态
            pso->updateStopState();

            // 更新stop code 与 stop flag
            pso->flushStopFlag();

            iter_flag.set_all(false);
            //iter_ready = true;

            if (pso->group_state.sleep.signal == true)
            {
                pso->group_state.sleep.result = true;
            }

            // 调用外部监听函数
            if (pso->monitor)
            {
                PSO_Info pso_info(pso->group_state.time, pso->group_state.nGene, pso->best_indiv);
                pso->monitor(pso_info);
            }
        }
        else
        {
            //iter_ready = false;
        }
    }
}

#endif
