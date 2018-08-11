#pragma once
#include <string>

std::string out_bool(const bool b)
{
	if (b)
	{
		return std::string("true");
	}
	else
	{
		return std::string("false");
	}
}
