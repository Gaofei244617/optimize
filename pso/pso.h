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

namespace opt
{
	template<class F> class PSO;

	// PSO�࣬���ṩ��Ӧ�Ⱥ�������
	template<class R, class... Args>
	class PSO<R(Args...)>
	{
		using Bound = std::initializer_list<std::initializer_list<double>>;

	private:
		std::size_t groupSize;                                                // ��������
		std::size_t nVars;                                                    // ��Ӧ�Ⱥ��������ı�������
		PSO_Individual* indivs;                                               // ����Ⱥ����
		R(*fitFunc)(Args...);                                                 // ��Ӧ�Ⱥ���ָ��
		std::function<double(std::size_t)> reweight;                          // ��Ⱥ������̬����
		std::function<void(const PSO_Info&)> monitor;                         // �ⲿ������
		double(*bound)[2];                                                    // ÿ������������, ������ָ���ʾ
		std::vector<PSO_Individual> bestIndivs;                               // ��¼ÿ�ε��������Ÿ���
		PSO_Individual best_indiv;                                            // ��������
		PSO_State group_state;                                                // ����Ⱥ״̬

	public:
		PSO(R(*f)(Args...), const std::size_t size = 1000);                   // ���캯��������һ����Ⱥ�����ṩ��Ӧ�Ⱥ�������Ⱥ����
		PSO(PSO<R(Args...)>& other);                                          // ��������
		PSO(PSO<R(Args...)>&& other);                                         // �ƶ�����
		PSO<R(Args...)>& operator=(const PSO<R(Args...)>& other) = delete;
		PSO<R(Args...)>& operator=(PSO<R(Args...)>&& other) = delete;
		~PSO();

		void setBoundary(double(*b)[2]);                                      // ���ñ�������
		void setBoundary(const Bound& b);                                     // ���ñ�������
		void setMaxGeneration(const unsigned int N);                          // ��������������
		void setMaxRuntime(const Second& time);                               // �����������ʱ�䣨�룩
		void setStopTol(const long double t, const unsigned int N = 5);       // �������Ž�ֹͣ���
		void setMonitor(const std::function<void(const PSO_Info&)>&);         // �����ⲿ��������
		void setReweight(const std::function<std::size_t(std::size_t)>&);     // ��Ⱥ������������

		int getNVars()const;                                                  // ��ȡ��Ⱥ��������
		int getGeneration()const;                                             // ��õ�ǰ��Ⱥ����
		int getGroupSize()const;                                              // ��õ�ǰ��Ⱥ��������
		std::vector<PSO_Individual> getBestIndivs();                          // ��ȡ���ε��������Ž�
		int getStopCode();                                                    // ��ȡStop Code

		void initGroup(const std::vector<PSO_Individual>& indivs = std::vector<PSO_Individual>());    // ��ʼ������Ⱥ
		bool start();                                                         // ��ʼ��������
		void pause();                                                         // ��ͣ(Ϊ��֤����һ���ԣ�����һ������������pause)
		void proceed();                                                       // ��������
		void kill();                                                          // ��������
		PSO<R(Args...)> clone();                                              // ��¡��ǰ��Ⱥ

	private:
		template<std::size_t... I>
		R callFitFunc(double* args, const opt::index_seq<I...>&);              // ������Ӧ�Ⱥ�����

		std::size_t findBest();                                                // Ѱ�����Ÿ���
		void iter();                                                           // ����һ�ε���
		void updateStopState();                                                // ��������Ⱥֹͣ״̬
		void run();
	};

	/******************************************* ���������� ***********************************************************/
	// ���캯�������ṩ����Ӧ�Ⱥ���������Ⱥ����
	template<class R, class... Args>
	PSO<R(Args...)>::PSO(R(*f)(Args...), const std::size_t size)
		:groupSize(size),
		nVars(sizeof...(Args)),
		indivs(new PSO_Individual[groupSize]),
		fitFunc(f),
		reweight([](std::size_t n) {return 0.5; }),
		bound(nullptr)
	{
		for (std::size_t i = 0; i < groupSize; i++)
		{
			indivs[i] = PSO_Individual(nVars);
		}

		// ��������
		bound = new double[nVars][2];
	}

	// ���ƹ���
	template<class R, class... Args>
	PSO<R(Args...)>::PSO(PSO<R(Args...)>& other)
		:groupSize(other.groupSize),
		nVars(other.nVars),
		indivs(new PSO_Individual[groupSize]),
		fitFunc(other.fitFunc),
		monitor(other.monitor),
		reweight(other.reweight),
		bound(nullptr),
		bestIndivs()
	{
		other.pause();

		for (std::size_t i = 0; i < groupSize; i++)
		{
			indivs[i] = (other.indivs)[i];
		}

		// ��������
		bound = new double[nVars][2];
		for (std::size_t i = 0; i < nVars; i++)
		{
			bound[i][0] = (other.bound)[i][0];
			bound[i][1] = (other.bound)[i][1];
		}

		this->bestIndivs = other.bestIndivs;

		other.proceed();
	}

	// �ƶ�����
	template<class R, class... Args>
	PSO<R(Args...)>::PSO(PSO<R(Args...)>&& other)
		:groupSize(other.groupSize),
		nVars(other.nVars),
		indivs(other.indivs),
		fitFunc(other.fitFunc),
		monitor(other.monitor),
		reweight(other.reweight),
		bound(other.bound),
		bestIndivs()
	{
		other.kill();        // ��ֹ��Ⱥ����

		other.indivs = nullptr;
		other.bound = nullptr;

		this->bestIndivs = std::move(other.bestIndivs);
	}

	// ��������
	template<class R, class... Args>
	PSO<R(Args...)>::~PSO()
	{
		delete[] indivs;
		delete[] bound;
	}

	/**************************************** Setter ************************************************************/
	// �������б�������,�����ʼ���б�
	template<class R, class... Args>
	void PSO<R(Args...)>::setBoundary(const Bound& b)
	{
		const std::size_t len = b.size();
		for (std::size_t i = 0; i < len; i++)
		{
			bound[i][0] = *((*(b.begin() + i)).begin());         // �±߽�
			bound[i][1] = *((*(b.begin() + i)).begin() + 1);     // �ϱ߽�
		}
		group_state.setBoundFlag = true;
	}

	// �������б�������,��������ָ��
	template<class R, class... Args>
	void PSO<R(Args...)>::setBoundary(double(*b)[2])
	{
		for (int i = 0; i < nVars; i++)
		{
			bound[i][0] = b[i][0];       // �±߽�
			bound[i][1] = b[i][1];       // �ϱ߽�
		}
		group_state.setBoundFlag = true;
	}

	// ��������������
	template<class R, class... Args>
	void PSO<R(Args...)>::setMaxGeneration(const unsigned int N)
	{
		group_state.setMaxGeneFlag = true;
		group_state.maxGene = N;
	}

	// �����������ʱ�䣨�룩
	template<class R, class... Args>
	void PSO<R(Args...)>::setMaxRuntime(const Second& time)
	{
		group_state.setRuntimeFlag = true;
		group_state.maxRuntime = time;
	}

	// �������Ž�ֹͣ���
	template<class R, class... Args>
	void PSO<R(Args...)>::setStopTol(const long double t, const unsigned int N)
	{
		group_state.converCount = N;
		group_state.setStopTolFlag = true;
		group_state.stopTol = t;
	}

	// �����ⲿ��������
	template<class R, class... Args>
	void PSO<R(Args...)>::setMonitor(const std::function<void(const PSO_Info&)>& func)
	{
		this->monitor = func;
	}

	// ��Ⱥ������������
	template<class R, class... Args>
	void PSO<R(Args...)>::setReweight(const std::function<std::size_t(std::size_t)>& func)
	{
		this->reweight = func;
	}

	/**************************************Getter*********************************************************************/
	// ��ȡ��Ⱥ��������
	template<class R, class... Args>
	int PSO<R(Args...)>::getNVars()const
	{
		return nVars;
	}

	// ��ȡǰ��Ⱥ����
	template<class R, class... Args>
	int PSO<R(Args...)>::getGeneration()const
	{
		return group_state.nGene;
	}

	// ��ȡ��ǰ��Ⱥ��������
	template<class R, class... Args>
	int PSO<R(Args...)>::getGroupSize()const
	{
		return groupSize;
	}

	// ��ȡ���ε��������Ž�
	template<class R, class... Args>
	std::vector<PSO_Individual> PSO<R(Args...)>::getBestIndivs()
	{
		return bestIndivs;
	}

	// ��ȡֹͣ����
	template<class R, class... Args>
	int PSO<R(Args...)>::getStopCode()
	{
		return group_state.stopCode;
	}

	// ��ʼ������Ⱥ
	template<class R, class... Args>
	void PSO<R(Args...)>::initGroup(const std::vector<PSO_Individual>& init_indivs)
	{
		// ��ʼ��ǰ��Ҫ��֤�����û����������
		if (group_state.setBoundFlag)
		{
			// 1.��ʼ��ÿһ������
			for (std::size_t i = 0; i < init_indivs.size(); i++)
			{
				indivs[i] = init_indivs[i];
				indivs[i].fitness = callFitFunc(indivs[i].xs, opt::make_index_seq<sizeof...(Args)>());
			}
			for (std::size_t i = init_indivs.size(); i < groupSize; i++)
			{
				// ���ø���λ��
				for (std::size_t j = 0; j < nVars; j++)
				{
					indivs[i].xs[j] = random_real(bound[j][0], bound[j][1]);
					indivs[i].best_xs[j] = indivs[i].xs[j];
				}
				// ������Ӧ��
				indivs[i].fitness = callFitFunc(indivs[i].xs, opt::make_index_seq<sizeof...(Args)>());
			}

			// 2.Ѱ�����Ÿ���
			std::size_t bestIndex = this->findBest();
			bestIndivs.push_back(indivs[bestIndex]);                     // ��¼���Ž�
			best_indiv = indivs[bestIndex];                              // ���Ÿ���

			// 3.�������ӷ����ٶ�			
			double c2 = 2;   // ѧϰ����

			// ����������ӷ����ٶȣ���ʼ�������ٶ�û�й��Բ��ֺ������֣�ֻ����Ჿ��
			for (int i = 0; i < groupSize; i++)
			{
				// �����
				double r2 = opt::random_real(0, 1);
				for (int j = 0; j < nVars; j++)
				{
					indivs[i].vs[j] = c2 * r2 * (indivs[bestIndex].xs[j] - indivs[i].xs[j]);

					// �������ٶ�
					double V_max = 0.15 * (bound[j][1] - bound[j][0]);

					// ����ٶ��Ƿ����
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

			// 4.���ó�ʼ����־λ
			group_state.initFlag = true;
		}
		else
		{
			throw std::string("Fitness function or variables boundary is not set.");
		}
	}

	// Ѱ�����Ÿ���
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

	// ����һ�ε���
	template<class R, class ...Args>
	void PSO<R(Args...)>::iter()
	{
		double weight = reweight(group_state.nGene + 1);         // ����ϵ��

		double c1 = 1.49618;                                     // ѧϰ����
		double c2 = 1.49618;

		double r1 = 0.5;
		double r2 = 0.5;

		double V_max = 0;                                        // �������ٶ�

		double temp = 0;
		double temp_1 = 0;
		double temp_2 = 0;
		double temp_3 = 0;

		// �����ٶ�
		for (int j = 0; j < nVars; j++)
		{
			V_max = 0.15 * (bound[j][1] - bound[j][0]);

			for (int i = 0; i < groupSize; i++)
			{
				r1 = opt::random_real(0, 1);
				r2 = opt::random_real(0, 1);

				temp_1 = weight * indivs[i].vs[j];                                // ���Բ���
				temp_2 = c1 * r1 * (indivs[i].best_xs[j] - indivs[i].xs[j]);      // ���Ҳ���
				temp_3 = c2 * r2 * (best_indiv.xs[j] - indivs[i].xs[j]);          // ��Ჿ��

				temp = temp_1 + temp_2 + temp_3;                                  // �����ٶ�

				if (temp > V_max)
				{
					temp = V_max;
				}
				if (temp < -V_max)
				{
					temp = -V_max;
				}

				// ����λ��
				indivs[i].xs[j] = indivs[i].xs[j] + indivs[i].vs[j];
				if (indivs[i].xs[j] < bound[j][0])
				{
					indivs[i].xs[j] = bound[j][0];
				}
				if (indivs[i].xs[j] > bound[j][1])
				{
					indivs[i].xs[j] = bound[j][1];
				}
				// �����ٶ�
				indivs[i].vs[j] = temp;
			}
		}

		// ����������Ӧ��
		for (int i = 0; i < groupSize; i++)
		{
			indivs[i].fitness = callFitFunc(indivs[i].xs, opt::make_index_seq<sizeof...(Args)>());
			double best_fit = callFitFunc(indivs[i].best_xs, opt::make_index_seq<sizeof...(Args)>());
			if (indivs[i].fitness > best_fit)
			{
				for (int j = 0; j < nVars; j++)
				{
					indivs[i].best_xs[j] = indivs[i].xs[j];
				}
			}
		}

		std::size_t bestIndex = this->findBest();                    // Ѱ�ҵ�ǰ����Ⱥ���Ÿ���
		bestIndivs.push_back(indivs[bestIndex]);                     // ��¼���Ž�

		if (indivs[bestIndex].fitness > best_indiv.fitness)          // ���Ÿ�����
		{
			best_indiv = indivs[bestIndex];
		}
	}

	// ��������Ⱥֹͣ״̬
	template<class R, class ...Args>
	void PSO<R(Args...)>::updateStopState()
	{
		// ����������һ
		group_state.nGene++;

		// ��¼��ǰʱ��
		group_state.nowTime = std::chrono::steady_clock::now();

		// �Ƿ�ﵽ������ʱ��
		std::chrono::duration<double> evolTime = group_state.nowTime - group_state.startTime;
		group_state.time = evolTime.count(); // ��

		// �ж����Ÿ���fitnessֵ����һ���Ĳ������
		if (abs(bestIndivs[group_state.nGene - 1].fitness - bestIndivs[group_state.nGene].fitness) <= group_state.stopTol)
		{
			group_state.count++;
		}
		else
		{
			group_state.count = 0;
		}
	}

	// ������Ӧ�Ⱥ���
	template<class R, class... Args>
	template<std::size_t... I>
	inline R PSO<R(Args...)>::callFitFunc(double* args, const opt::index_seq<I...>&)
	{
		return fitFunc(args[I]...);
	}
}

#endif