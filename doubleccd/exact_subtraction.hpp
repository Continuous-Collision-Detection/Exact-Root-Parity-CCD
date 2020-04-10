// Require IEEE-754 compliance
// To compile on MSVC: use /fp:strict
// On GNU GCC: use -frounding-math
// In any case: compile for 64bit target to avoid X87 ugly extended precision

#include <utility>
#include <fenv.h>

#ifdef _MSC_VER

#pragma fenv_access (on)

inline void _setFPUModeToRoundUP() { _controlfp(_RC_UP, _MCW_RC); }
inline void _setFPUModeToRoundNEAR() { _controlfp(_RC_NEAR, _MCW_RC); }

#else

#pragma STDC FENV_ACCESS ON

inline void _setFPUModeToRoundUP() { fesetround(FE_UPWARD); }
inline void _setFPUModeToRoundNEAR() { fesetround(FE_TONEAREST); }
#endif

static double sterbenzDisplacement(double x, double y)
{
	// Assume we are rounding towards +infinity
	if (x == 0 || y == 0 || x == y) return 0.0;
	if (x > 0 && y > 0 && 0.5 * y <= x && x <= 2 * y) return 0.0;
	if (x < 0 && y < 0 && 0.5 * y >= x && x >= 2 * y) return 0.0;
	return std::fmax((y - 2 * x), (x - 2 * y));
}

static double maxCommonDisplacement(const std::vector<std::pair<double, double>>& subtractions)
{
	double ad, cd = 0.0;
	for (std::pair<double, double> p : subtractions)
		if ((ad = sterbenzDisplacement(p.first, p.second)) > cd) cd = ad;
	return cd;
}

// This function takes a vector of pairs <x_i, y_i>, calculates a value z,
// and increases all the values in the vector by z.
// In the modified vector, we are guaranteed that x_i - y_i is error free for 
// all the pairs.
static void displaceSubtractions(std::vector<std::pair<double, double>>& subtractions)
{
	_setFPUModeToRoundUP();

	double z = maxCommonDisplacement(subtractions);
	for (std::pair<double, double>& p : subtractions)
	{
		p.first += z;
		p.second += z;
	}

	_setFPUModeToRoundNEAR();
}

// This function takes a vector of pairs <x_i, y_i>, calculates a value z,
// and adds z to all the values in the vector. Then, subtracts z from all
// the values of the vector. This induces a loss of precision in the initial
// values guaranteeing that x_i - y_i is error free for all the pairs.
static void perturbSubtractions(std::vector<std::pair<double, double>>& subtractions)
{
	_setFPUModeToRoundUP();

	double z = maxCommonDisplacement(subtractions);
	for (std::pair<double, double>& p : subtractions)
	{
		p.first += z;
		p.second += z;
	}

	for (std::pair<double, double>& p : subtractions)
	{
		p.first -= z;
		p.second -= z;
	}

	_setFPUModeToRoundNEAR();
}
static double displaceSubtractions_double(std::vector<std::pair<double, double>>& subtractions)
{
    _setFPUModeToRoundUP();

    double z = maxCommonDisplacement(subtractions);
    for (std::pair<double, double>& p : subtractions) {
        p.first += z;
        p.second += z;
    }

    _setFPUModeToRoundNEAR();
    return z;
}