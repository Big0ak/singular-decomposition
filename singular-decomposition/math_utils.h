#pragma once

namespace MATH_utils
{
	double cube_root(double x)
	{
		if (x > 0)
		{
			return std::pow(x, 1.0 / 3.0);
		}
		else {
			return -std::pow(-x, 1.0 / 3.0);
		}
	}

	double sign(double x)
	{
		if (x < 0)
			return -1.;

		if (x > 0)
			return 1.;

		return 0;
	}
}