# ��Ⱥ����

```cpp
namespace opt
{
    struct Individual
    {
        const int nVars;                                   // ����(����)����
        double* vars;                                      // �����������(����)
        double fitness;                                    // ����Ի�������Ӧ��

        Individual(int n);                                 // ���캯��
        Individual(const Individual& other);               // ��������
        Individual& operator=(const Individual& other);    // ��ֵ����
        ~Individual();                                     // ����
    };
}
```