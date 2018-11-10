#ifndef _ROULETTE_H_
#define _ROULETTE_H_

#include "range_random.h"

namespace opt
{
	// ���̶ķ�װ��
	template<class T>
	class Roulette
	{
	private:
		int mark_num;                                           // �̶�������
		T* tick_mark;                                           // ���̶̿���(�豣֤��ֵ��������)

	public:
		Roulette();                                             // ���캯�� 
		Roulette(const int N);                                  // ���캯�� 
		Roulette(const Roulette<T>& other);                     // ��������
		Roulette(Roulette<T>&& other);                          // �ƶ�����
		~Roulette();                                            // ��������

		int roll();                                             // �������ת��
		//void setMark(const int INDEX, const T& NUM);          // ����̶��ߵ�ֵ
		T& operator[](const int& i);
	};

	// ���캯��
	template<class T>
	Roulette<T>::Roulette() 
		:mark_num(0),
		tick_mark(nullptr)
	{}
	
	// ���캯��
	template<class T>
	Roulette<T>::Roulette(const int N)
		: mark_num(N),
		tick_mark(new T[N]())
	{}

	// ���ƹ���
	template<class T>
	Roulette<T>::Roulette(const Roulette<T>& other)
		:mark_num(other.mark_num),
		tick_mark(new T[mark_num]())
	{
		for (int i = 0; i < mark_num; i++)
		{
			tick_mark[i] = other.tick_mark[i];
		}
	}

	// �ƶ�����
	template<class T>
	Roulette<T>::Roulette(Roulette<T>&& other)
		:mark_num(other.mark_num),
		tick_mark(other.tick_mark)
	{
		other.tick_mark = nullptr;
	}

	// ��������
	template<class T>
	Roulette<T>::~Roulette()
	{
		delete[] tick_mark;
	}

	// �������ת��ָ��һ������
	template<class T>
	int Roulette<T>::roll()
	{
		if (mark_num > 0)
		{
			double temp = random_real(tick_mark[0], tick_mark[mark_num - 1]);

			int lowIndex = 0;
			int upIndex = mark_num - 1;
			int tempIndex;

			// ���ַ�����,ʱ�临�Ӷ�log_n
			// ͨ�����ƶ�lowIndex��upIndex������±�λ�ã�ÿ�ν����������Сһ��
			while (upIndex - lowIndex > 1)
			{
				tempIndex = (lowIndex + upIndex) / 2 + (lowIndex + upIndex) % 2;
				if (temp <= tick_mark[tempIndex])
				{
					upIndex = tempIndex;
				}
				else
				{
					lowIndex = tempIndex;
				}
			}
			return lowIndex;  // up = 1, low = 0 ʱ���������indivs[0]
		}
		else
		{
			return 0;
		}
	}

	//// ����̶��ߵ�ֵ
	//template<class T>
	//inline void Roulette<T>::setMark(const int INDEX, const T& NUM)
	//{
	//	tick_mark[INDEX] = NUM;
	//}

	// ��ȡ���̶Ŀ̶���
	template<class T>
	T& Roulette<T>::operator[](const int& i)
	{
		return tick_mark[i];
	}
}

#endif