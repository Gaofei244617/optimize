#ifndef _INDEX_SEQ_
#define _INDEX_SEQ_

namespace opt
{
	// �����ڽ�Nչ��Ϊ(0,1,...,N-1)

	template<std::size_t... I>
	struct index_seq {};

	template<std::size_t N, std::size_t... I>
	struct make_index_seq :public make_index_seq<N - 1, N - 1, I...> {};

	template<std::size_t... I>
	struct make_index_seq<0, I...> :public index_seq<I...> {};
}

#endif