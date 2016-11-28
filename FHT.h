// Split-Radix Fast Hartley Transform
// Inplace with power of 2 data length
// Author: 2015-2016 Eugeniy Sokol, Magnitogorsk

#pragma once

#include <cmath>
#include <complex>

#ifndef consts_pi
#define consts_pi
const double pi =    3.14159265358979323846264338327950;
const double pi4 =   0.78539816339744830961566084581988;
const double sqrt2 = 1.41421356237309504880168872420970;
const double cos45 = 0.70710678118654752440084436210485;
const double cos22 = 0.92387953251128675612818318939679;
const double sin22 = 0.38268343236508977172845998403040;
#endif

class FHT
{
public:
	FHT::FHT(int size);
	FHT::~FHT();

	void FHT::transform(double* data, bool scaled=true);
	void FHT::back_transform(double* fht_data);
	void FHT::convolute(double* fht_data, double* fht_fir, bool scaled=false);
	std::complex<double> FHT::get_frequency(double* fht_data, int number);

private:
	static int fht_instance_count[];
	static double* fht_trig_tables[];
	static int* fht_revbin_tables[];
	static int fht_revbin_counts[];

	int size;
	int ldn;
	//--------------------------------------------
	void FHT::revbin_permute(double* data);
	void FHT::step_addsub(double* data, int ldn);
	void FHT::step_rotate(double* data, int ldn);
	//--------------------------------------------
	__inline void FHT::step_addsub4(double* data);
	__inline void FHT::step_addsub8(double* data);
	__inline void FHT::step_rotate8(double* data);
	__inline void FHT::addsub(double& u, double& v);
	__inline void FHT::addsub(double& u, double& v, const double scale);
	__inline void FHT::rotate(double& u, double& v, const double c, const double s);
};