#ifndef _PSO_STATE_H_
#define _PSO_STATE_H_

#include <chrono>
#include <atomic>
#include "..\opt_time.h"

namespace opt
{
    // 粒子状态标（随着迭代过程可能发生变化）
    struct PSO_State
    {
        long double stopTol;                                                  // 停止迭代误差
        unsigned int maxGene;                                                 // 最大迭代次数
        unsigned int converCount;                                             // 停止迭代的收敛次数
        Second maxRuntime;                                                    // 最大迭代时间(秒)

        bool setBoundFlag;                                                    // 是否设置变量区间
        bool setRuntimeFlag;                                                  // 是否设置最大运行时间
        bool setStopTolFlag;                                                  // 是否设置停止误差
        bool setMaxGeneFlag;                                                  // 是否设置最大迭代次数

        bool initFlag;                                                        // 是否初始化粒子群
        bool stopFlag;                                                        // 迭代停止标志, true:达到停止条件, false:未达到停止条件
        SleepFlag sleep;                                                      //
        std::atomic<short> stopCode;                                          // 迭代停止原因,-1-未开始迭代; 0-正在迭代; 1-达到最大迭代次数; 2-达到最大迭代时间; 3-最优解收敛于稳定值; 4-人为停止迭代
        unsigned int count;                                                   // 最优解波动连续小于停止误差的次数, 波动值连续5次小于停止误差即停止迭代

        unsigned int nGene;                                                   // 当前粒子群代数
        std::chrono::time_point<std::chrono::steady_clock> startTime;         // 开始迭代的时间点
        std::chrono::time_point<std::chrono::steady_clock> nowTime;           // 当前迭代的时间点
        double time;                                                          // 粒子群迭代时间：秒
        //int worstIndex;                                                       // 适应度最差的个体位置
        //int bestIndex;                                                        // 适应度最好的个体位置

        /* 构造函数*/
        PSO_State()
            :stopTol(0),
            maxGene(0),
            converCount(5),
            maxRuntime(0),
            setBoundFlag(false),
            setRuntimeFlag(false),
            setStopTolFlag(false),
            setMaxGeneFlag(false),
            initFlag(false),
            stopFlag(false),
            sleep(),
            stopCode(-1),
            count(0),
            nGene(0),
            time(0)
            //worstIndex(0),
            //bestIndex(0)
        {}

        // 复制构造
        PSO_State(const PSO_State& other)
            :stopTol(other.stopTol),
            maxGene(other.maxGene),
            converCount(other.converCount),
            maxRuntime(other.maxRuntime),
            setBoundFlag(other.setBoundFlag),
            setRuntimeFlag(other.setRuntimeFlag),
            setStopTolFlag(other.setStopTolFlag),
            setMaxGeneFlag(other.setMaxGeneFlag),
            initFlag(other.initFlag),
            stopFlag(other.stopFlag),
            sleep(other.sleep),
            stopCode(other.stopCode.load()),
            count(other.count),
            nGene(other.nGene),
            time(other.time)
            //worstIndex(other.worstIndex),
            //bestIndex(other.bestIndex)
        {}

        // 赋值函数
        PSO_State& operator=(const PSO_State& other)
        {
            if (this != &other)
            {
                stopTol = other.stopTol;
                maxGene = other.maxGene;
                converCount = other.converCount;
                maxRuntime = other.maxRuntime;
                setBoundFlag = other.setBoundFlag;
                setRuntimeFlag = other.setRuntimeFlag;
                setStopTolFlag = other.setStopTolFlag;
                setMaxGeneFlag = other.setMaxGeneFlag;
                initFlag = other.initFlag;
                stopFlag = other.stopFlag;
                sleep = other.sleep;
                stopCode = other.stopCode.load();
                count = other.count;
                nGene = other.nGene;
                time = other.time;
                //worstIndex = other.worstIndex;
                //bestIndex = other.bestIndex;
            }

            return *this;
        }

        // 判断粒子群是否处于可迭代状态
        // 可迭代条件: (1)设置了变量区间(2)且完成粒子群初始化(3)且至少设置了最大运行时间、停止误差、最大迭代次数之一
        bool runable()
        {
            return setBoundFlag && initFlag && (setRuntimeFlag || setStopTolFlag || setMaxGeneFlag);
        }
    };
}

#endif
