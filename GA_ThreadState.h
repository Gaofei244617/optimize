#ifndef _GA_THREADSTATE_H_
#define _GA_THREADSTATE_H_

#include <thread>
#include <vector>
#include <utility>
#include <condition_variable>
#include <atomic>
#include "bool_array.h"

namespace opt
{
	// GA�㷨���м������״̬��־
	struct GA_ThreadState
	{
		std::mutex mtx;  // ������
		std::condition_variable cv; // ��������, �����߳�ͬ��
		bool crossReady;  // ���������־
		bool mutReady;    // ���������־
		bool selectReady; // ѡ�������־
		int threadNum;    // ���м���ʹ�õ��߳���
		bool_array cross_flag; // ���м���: �����߳�״̬��־
		bool_array mut_flag;   // ���м���: �����߳�״̬��־

		// ���캯��
		GA_ThreadState()
			:crossReady(false),
			mutReady(false),
			selectReady(false),
			threadNum(1),
			cross_flag(threadNum),
			mut_flag(threadNum)
		{}

		// ���ƹ���
		GA_ThreadState(const GA_ThreadState& other)
			:crossReady(other.crossReady),
			mutReady(other.mutReady),
			selectReady(other.selectReady),
			threadNum(other.threadNum),
			cross_flag(other.cross_flag),
			mut_flag(other.mut_flag)
		{}
	};
}

#endif
