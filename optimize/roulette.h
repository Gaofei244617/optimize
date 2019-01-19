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
		std::size_t mark_num;                                   // �̶�������
		T* tick_mark;                                           // ���̶̿���(�豣֤��ֵ��������)

	public:
		Roulette();                                             // ���캯��
		Roulette(const std::size_t N);                          // ���캯��
		Roulette(const Roulette<T>& other);                     // ��������
		Roulette(Roulette<T>&& other);                          // �ƶ�����
		Roulette<T>& operator=(const Roulette<T>& other);       // ��ֵ����
		Roulette<T>& operator=(Roulette<T>&& other);            // �ƶ���ֵ
		~Roulette();                                            // ��������

		int roll();                                             // �������ת��
		T& operator[](const std::size_t& i);
		void reset(const std::size_t N);                        // ���»��̶ֿ���
	};

	// ���캯��
	template<class T>
	Roulette<T>::Roulette()
		:mark_num(0),
		tick_mark(nullptr)
	{}

	// ���캯��
	template<class T>
	Roulette<T>::Roulette(const std::size_t N)
		: mark_num(N),
		tick_mark(new T[N]())
	{}

	// ���ƹ���
	template<class T>
	Roulette<T>::Roulette(const Roulette<T>& other)
		:mark_num(other.mark_num),
		tick_mark(new T[mark_num]())
	{
		for (std::size_t i = 0; i < mark_num; i++)
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

	// ��ֵ����
	template<class T>
	Roulette<T>& Roulette<T>::operator=(const Roulette<T>& other)
	{
		if (this != &other)
		{
			mark_num = other.mark_num;
			tick_mark = new T[mark_num]();
			for (std::size_t i = 0; i < mark_num; i++)
			{
				tick_mark[i] = other.tick_mark[i];
			}
		}
		return *this;
	}

	// �ƶ���ֵ
	template<class T>
	Roulette<T>& Roulette<T>::operator=(Roulette<T>&& other)
	{
		if (this != &other)
		{
			mark_num = other.mark_num;
			tick_mark = other.tick_mark;
			other.tick_mark = nullptr;
		}
		return *this;
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

	// ��ȡ���̶Ŀ̶���
	template<class T>
	T& Roulette<T>::operator[](const std::size_t& i)
	{
		return tick_mark[i];
	}

	//
	template<class T>
	void Roulette<T>::reset(const std::size_t N)
	{
		if (N != mark_num)
		{
			delete[] tick_mark;
			tick_mark = new T[N]();
			mark_num = N;
		}
	}
}

#endif