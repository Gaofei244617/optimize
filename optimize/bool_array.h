#ifndef _BOOL_ARRAY_
#define _BOOL_ARRAY_

namespace opt
{
	// bool����������
	class bool_array
	{
	private:
		bool* array;
		std::size_t Length;
		std::size_t Capacity;

	public:
		// ���캯��������Ԫ�س�ʼ��Ϊfalse
		bool_array() :array(nullptr), Length(0), Capacity(0) {}

		bool_array(const std::size_t N)
			:array(new bool[N]),
			Length(N),
			Capacity(N)
		{
			for (std::size_t i = 0; i < N; i++)
			{
				array[i] = false;
			}
		}

		// ���ƹ���
		bool_array(const bool_array& other)
			:array(new bool[other.Length]),
			Length(other.Length),
			Capacity(Length)
		{
			for (std::size_t i = 0; i < Length; i++)
			{
				array[i] = (other.array)[i];
			}
		}

		// �ж��Ƿ�����Ԫ��Ϊtrue
		bool is_all_true()const
		{
			if (Length > 0)
			{
				bool temp = true;
				for (std::size_t i = 0; i < Length; i++)
				{
					temp = temp && array[i];
				}
				return temp;
			}
			return false;
		}

		// Ϊ����Ԫ��������ֵͬ
		void set_all(const bool b)
		{
			for (std::size_t i = 0; i < Length; i++)
			{
				array[i] = b;
			}
		}

		// ����bool��������
		void set_length(const std::size_t N)
		{
			/*** �����������(1) N > Capacity; (2.1) Length < N <= Capacity; (2.2) N <= Length ***/
			// (1) N > Capacity�����,��Ҫ��Ҫ���·����ڴ�
			if (N > Capacity)
			{
				bool* temp = new bool[N];
				for (std::size_t i = 0; i < Capacity; i++)
				{
					temp[i] = array[i];
				}
				for (std::size_t i = Capacity; i < N; i++)
				{
					temp[i] = false;
				}

				Capacity = N;
				Length = N;

				bool* temp2 = temp;
				temp = array;
				array = temp2;
				delete[] temp;
				temp = temp2 = nullptr;
			}
			else
			{
				// (2.1) Length < N <= Capacity; (2.2) N <= Length, ����Ҫ���·����ڴ�
				if (N > Length)
				{
					for (std::size_t i = Length; i < N; i++)
					{
						array[i] = false;
					}
				}
				Length = N;
			}
		}

		// �±����������
		bool& operator[](const std::size_t n)
		{
			return array[n];
		}

		// ��������
		~bool_array()
		{
			delete[] array;
		}
	};
}

#endif