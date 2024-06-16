#pragma once

#include "cmath"

#define M_PI	3.14159265358979323846

namespace MATH_utils
{
	inline double cube_root(double x)
	{
		if (x > 0)
		{
			return std::pow(x, 1.0 / 3.0);
		}
		else {
			return -std::pow(-x, 1.0 / 3.0);
		}
	}

	inline double sign(double x)
	{
		if (x < 0)
			return -1.;

		if (x > 0)
			return 1.;

		return 0;
	}
}
