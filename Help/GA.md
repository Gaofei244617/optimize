# GA�ࣺGAGroup

## �ඨ��

$$a^2=b^2+c^2$$
$$\int_{0}^{\frac{\pi}{2}}(x^2+3)dx$$

```cpp
namespace opt
{
    template<class F> class GAGroup; // �޾���ʵ��
    template<class R, class... Args> class GAGroup<R(Args...)>;  // ʹ����Ӧ�Ⱥ���ƫ�ػ�GAGroup
}
```

## ���캯��

```cpp
// ���캯����f:��Ӧ�Ⱥ���,size:��������
GAGroup(R(*f)(Args...), const int size = 1000);
// ��������
GAGroup(const GAGroup<R(Args...)>& other);  
// �ƶ�����
GAGroup(GAGroup<R(Args...)>&& other);
```

## ��������

```cpp
~GAGroup();
```

## ��Ա����

### (1) ������Ⱥ����

```cpp
void setName(const std::string& str);
```

������Ⱥ���ơ�

```cpp
void setBoundary(double(*b)[2]);
```

���ñ������䡣

 ```cpp
void setBoundary(const GenBound& b);
```

���ñ������䡣

 ```cpp
void setMaxGeneration(const unsigned int N);
```

����������������

```cpp
void setMaxRuntime(const Second& time);
```

�����������ʱ�䣨�룩��

```cpp
void setStopTol(long double t, unsigned int N = 5);
```

�������Ž�ֹͣ��

```cpp
void setMutateProb(double p);
```

���û��������ʡ�

```cpp
void setCrossProb(double p);
```

���ý�����ʡ�

```cpp
void setThreadNum(const int NUM);
```

���ò��м�����߳�����default=1��

### (2) ��ȡ��Ⱥ����

```cpp
const std::string getName()const;
```

��ȡ��Ⱥ���ơ�

```cpp
int getNVars()const;
```

��ȡ��Ⱥ����������

```cpp
int getGeneration()const;
```

��õ�ǰ��Ⱥ������

```cpp
int getGroupSize()const;
```

��õ�ǰ��Ⱥ����������

```cpp
std::vector<Individual> getBestIndivs();
```

��ȡ���ε��������Ž⡣

```cpp
int getStopCode();
```

��ȡStop Code��

### (3) �û�����

```cpp
bool start();
```

���������������̡�

```cpp
void wait_result();
```

������ǰ�̣߳��ȴ��Ż������