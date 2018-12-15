#ifndef _PSO_STATE_H_
#define _PSO_STATE_H_

#include <chrono>
#include <atomic>
#include "opt_time.h"

namespace opt
{
	//struct SleepFlag
	//{
	//	bool signal;
	//	bool result;
	//	SleepFlag() :signal(false), result(false) {}
	//};

	// ����״̬�꣨���ŵ������̿��ܷ����仯��
	struct PSO_State
	{
		long double stopTol;                                                  // ֹͣ�������
		unsigned int maxGene;                                                 // ����������
		unsigned int converCount;                                             // ֹͣ��������������
		Second maxRuntime;                                                    // ������ʱ��(��)

		bool setBoundFlag;                                                    // �Ƿ����ñ�������
		bool setRuntimeFlag;                                                  // �Ƿ������������ʱ��
		bool setStopTolFlag;                                                  // �Ƿ�����ֹͣ���
		bool setMaxGeneFlag;                                                  // �Ƿ���������������

		bool initFlag;                                                        // �Ƿ��ʼ������Ⱥ
		bool stopFlag;                                                        // ����ֹͣ��־, true:�ﵽֹͣ����, false:δ�ﵽֹͣ����
		//SleepFlag sleep;                                                      //
		std::atomic<short> stopCode;                                          // ����ֹͣԭ��,-1-δ��ʼ����; 0-���ڵ���; 1-�ﵽ����������; 2-�ﵽ������ʱ��; 3-���Ž��������ȶ�ֵ; 4-��Ϊֹͣ����
		unsigned int count;                                                   // ���ŽⲨ������С��ֹͣ���Ĵ���, ����ֵ����5��С��ֹͣ��ֹͣ����

		unsigned int nGene;                                                   // ��ǰ����Ⱥ����
		std::chrono::time_point<std::chrono::steady_clock> startTime;         // ��ʼ������ʱ���
		std::chrono::time_point<std::chrono::steady_clock> nowTime;           // ��ǰ������ʱ���
		double time;                                                          // ����Ⱥ����ʱ�䣺��
		//int worstIndex;                                                       // ��Ӧ�����ĸ���λ��
		//int bestIndex;                                                        // ��Ӧ����õĸ���λ��

		/* ���캯��*/
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
			//sleep(),
			stopCode(-1),
			count(0),
			nGene(0),
			time(0)
			//worstIndex(0),
			//bestIndex(0)
		{}

		// ���ƹ���
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
			//sleep(other.sleep),
			stopCode(other.stopCode.load()),
			count(other.count),
			nGene(other.nGene),
			time(other.time)
			//worstIndex(other.worstIndex),
			//bestIndex(other.bestIndex)
		{}

		// ��ֵ����
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
				//sleep = other.sleep;
				stopCode = other.stopCode.load();
				count = other.count;
				nGene = other.nGene;
				time = other.time;
				//worstIndex = other.worstIndex;
				//bestIndex = other.bestIndex;
			}

			return *this;
		}

		// �ж�����Ⱥ�Ƿ��ڿɵ���״̬
		// �ɵ�������: (1)�����˱�������(2)���������Ⱥ��ʼ��(3)�������������������ʱ�䡢ֹͣ������������֮һ
		bool runable()
		{
			return setBoundFlag && initFlag && (setRuntimeFlag || setStopTolFlag || setMaxGeneFlag);
		}
	};
}

#endif
