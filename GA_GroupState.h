#ifndef _GA_STATE_H_
#define _GA_STATE_H_

#include <chrono>
#include <atomic>
#include "opt_time.h"

namespace opt
{
	// 种群状态标（随着迭代过程可能发生变化）
	struct GA_GroupState
	{
		long double stopTol;                                                  // 停止迭代误差
		unsigned int maxGene;                                                 // 最大迭代次数
		unsigned int converCount;                                             // 停止迭代的收敛次数
		Second maxRuntime;                                                    // 最大迭代时间

		bool setBoundFlag;                                                    // 是否设置变量区间
		bool setRuntimeFlag;                                                  // 是否设置最大运行时间
		bool setStopTolFlag;                                                  // 是否设置停止误差
		bool setMaxGeneFlag;                                                  // 是否设置最大迭代次数
		bool initFlag;                                                        // 是否初始化种群
		std::atomic<bool> stopFlag;                                           // 迭代停止标志, true:达到停止条件, false:未达到停止条件
		std::atomic<short> stopCode;                                          // 迭代停止原因,-1:未停止; 0:最优解收敛于稳定值; 1:达到最大迭代次数; 2:达到最大迭代时间; 3.人为停止迭代
		unsigned int count;                                                   // 最优解波动连续小于停止误差的次数, 波动值连续5次小于停止误差即停止迭代

		unsigned int nGene;                                                   // 当前种群代数
		std::chrono::time_point<std::chrono::steady_clock> startTime;         // 开始迭代的时间点
		std::chrono::time_point<std::chrono::steady_clock> nowTime;           // 当前迭代的时间点
		int worstIndex;                                                       // 适应度最差的个体位置
		int bestIndex;                                                        // 适应度最好的个体位置


		/* 构造函数*/
		GA_GroupState()
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
			stopCode(-1),
			count(0),
			nGene(0),
			worstIndex(0),
			bestIndex(0)
		{}

		/* 复制构造, 因为存在atomic成员, 默认复制构造被删除*/
		GA_GroupState(const GA_GroupState& other)
			:stopTol(other.stopTol),
			maxGene(other.maxGene),
			converCount(other.converCount),
			maxRuntime(other.maxRuntime),
			setBoundFlag(other.setBoundFlag),
			setRuntimeFlag(other.setRuntimeFlag),
			setStopTolFlag(other.setStopTolFlag),
			setMaxGeneFlag(other.setMaxGeneFlag),
			initFlag(other.initFlag),
			stopFlag(other.stopFlag.load()),
			stopCode(other.stopCode.load()),
			count(other.count),
			nGene(other.nGene),
			worstIndex(other.worstIndex),
			bestIndex(other.bestIndex)
		{}
		
		// 判断种群是否处于可迭代状态
		// 可迭代条件: (1)设置了变量区间(2)且完成种群初始化(3)且至少设置了最大运行时间、停止误差、最大迭代次数之一
		bool runable()
		{
			return setBoundFlag && initFlag && (setRuntimeFlag || setStopTolFlag || setMaxGeneFlag);
		}
	};
}

#endif

