#ifndef _INDEX_SEQ_
#define _INDEX_SEQ_

namespace opt
{
	template<size_t... I>
	struct index_seq {};

	template<size_t N, size_t... I>
	struct make_index_seq :public make_index_seq<N - 1, N - 1, I...> {};

	template<size_t... I>
	struct make_index_seq<0, I...> :public index_seq<I...> {};
}

#endif