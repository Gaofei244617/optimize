#ifndef _GA_STATE_H_
#define _GA_STATE_H_

#include <chrono>
#include <atomic>
#include "GA_Time.h"

namespace opt
{
	// ��Ⱥ״̬�꣨���ŵ������̿��ܷ����仯��
	struct GA_GroupState
	{
		long double stopTol;                                                  // ֹͣ�������
		unsigned int maxGene;                                                 // ����������
		unsigned int converCount;                                             // ֹͣ��������������
		Second maxRuntime;                                                    // ������ʱ��

		bool setBoundFlag;                                                    // �Ƿ����ñ�������
		bool setRuntimeFlag;                                                  // �Ƿ������������ʱ��
		bool setStopTolFlag;                                                  // �Ƿ�����ֹͣ���
		bool setMaxGeneFlag;                                                  // �Ƿ���������������
		bool initFlag;                                                        // �Ƿ��ʼ����Ⱥ
		std::atomic<bool> stopFlag;                                           // ����ֹͣ��־, true:�ﵽֹͣ����, false:δ�ﵽֹͣ����
		std::atomic<short> stopCode;                                          // ����ֹͣԭ��,-1:δֹͣ; 0:���Ž��������ȶ�ֵ; 1:�ﵽ����������; 2:�ﵽ������ʱ��; 3.��Ϊֹͣ����
		unsigned int count;                                                   // ���ŽⲨ������С��ֹͣ���Ĵ���, ����ֵ����5��С��ֹͣ��ֹͣ����

		unsigned int nGene;                                                   // ��ǰ��Ⱥ����
		std::chrono::time_point<std::chrono::steady_clock> startTime;         // ��ʼ������ʱ���
		std::chrono::time_point<std::chrono::steady_clock> nowTime;           // ��ǰ������ʱ���
		int worstIndex;                                                       // ��Ӧ�����ĸ���λ��
		int bestIndex;                                                        // ��Ӧ����õĸ���λ��


		/* ���캯��*/
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

		/* ���ƹ���, ��Ϊ����atomic��Ա, Ĭ�ϸ��ƹ��챻ɾ��*/
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
		
		// �ж���Ⱥ�Ƿ��ڿɵ���״̬
		// �ɵ�������: (1)�����˱�������(2)�������Ⱥ��ʼ��(3)�������������������ʱ�䡢ֹͣ������������֮һ
		bool runable()
		{
			return setBoundFlag && initFlag && (setRuntimeFlag || setStopTolFlag || setMaxGeneFlag);
		}
	};
}

#endif // 
