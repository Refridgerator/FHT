#pragma once

// Fast Hartley Transform
// Inplace with power of 2 data length

#include <math.h>

#ifndef consts_pi
#define consts_pi
const double pi =    3.14159265358979323846264338327950;
const double pi4 =   0.78539816339744830961566084581988;
const double sqrt2 = 1.41421356237309504880168872420970;
const double cos45 = 0.70710678118654752440084436210485;
const double cos22 = 0.92387953251128675612818318939679;
const double sin22 = 0.38268343236508977172845998403040;
#endif

//------------------------------------------------------------
// Global Variables
//
int fht_instance_count[]=  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double* fht_trig_tables[]= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int* fht_revbin_tables[]=  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int fht_revbin_counts[]=   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//------------------------------------------------------------

class FHT
{
private:

public:
	int size;
	int ldn;

	FHT::FHT(int size)
	{
		this->size = size;
		ldn =(int)(log((double)size)/log(2.));
		fht_instance_count[ldn]++;

		int i, j, index, k;
		
		// trigonometry table
		k=2;
		for(int cnt=1;cnt<=size/4;cnt*=2)
		{
			if(fht_trig_tables[k]==0)
			{
				fht_trig_tables[k] = static_cast<double*>(_aligned_malloc(cnt*4*sizeof(double),64));
				double* trig_table = fht_trig_tables[k];

				index = 0;
				double nn = cnt*4;
				for (i = 1; i <= cnt; i++)
				{
					trig_table[index + 0] = sin(i * pi / nn);
					trig_table[index + 1] = cos(i * pi / nn);

					trig_table[index + 2] = sin(3.*i * pi / nn);
					trig_table[index + 3] = cos(3.*i * pi / nn);

					index+=4;
				}
			}
			k++;
		}

		// permute table
		if(fht_revbin_tables[ldn]==0)
		{
			// count of non-symmetric numbers
			int k=0;
			if(ldn%2==0)
				k=(int)(pow(2,ldn-1)-pow(2,ldn/2-1));
			else
				k=(int)(pow(2,ldn-1)-pow(2,(ldn+1)/2-1));

			fht_revbin_tables[ldn] = static_cast<int*>(_aligned_malloc(k*2*sizeof(int),64));
			fht_revbin_counts[ldn] = k;
			index = 0;
			for (int x = 1; x < size; x++)
			{
				int t=x, r=0;
				for(int i=0;i<ldn;i++)
				{
					r<<=1;
					r+=t&1;
					t>>=1; 
				}
				if(r>x)
				{
					fht_revbin_tables[ldn][index]=x;
					fht_revbin_tables[ldn][index+1]=r;
					index+=2;
				}
			}
		}
	}

	FHT::~FHT()
	{
		fht_instance_count[ldn]--;
		if(fht_instance_count[ldn]==0)
		{
			_aligned_free(fht_revbin_tables[ldn]);
			fht_revbin_tables[ldn]=0;
		}
		int k0 = 31;
		while(k0>0 && fht_instance_count[k0]==0) k0--;
		k0+=1;
		for(int k=k0;k<32;k++)
		{
			if(fht_trig_tables[k]!=0)
			{
				_aligned_free(fht_trig_tables[k]);
				fht_trig_tables[k] = 0;
			}
		}
	}

	//------------------------------------------------------------
	void FHT::transform(double* data, bool scaled=true)
	{
		r4step_addsub(data, ldn);
		if (scaled)
			for (int i = 0; i < size; i++) data[i] *= 1.0/(double)size;
		revbin_permute(data);
	}

	void FHT::back_transform(double* data)
	{
		r4step_addsub(data, ldn);
		revbin_permute(data);
	}

	void FHT::convolute(double* fht_data, double* fht_fir, bool scaled=false)
	{
		double mulf = 0.5*1.0;
		double mul0 = 1.0;
		if(scaled)
		{
			mulf*=size;
			mul0*=size;
		}
		int i = 1, j = size-1;
		do
		{
			double yp = (fht_data[i] + fht_data[j]) * mulf;
			double ym = (fht_data[i] - fht_data[j]) * mulf;
			fht_data[i] = (fht_fir[i] * yp + fht_fir[j] * ym);
			fht_data[j] = (fht_fir[j] * yp - fht_fir[i] * ym);
			i++; j--;
		}
		while (i<j);

		fht_data[0] *= fht_fir[0] * mul0;
		fht_data[size/2] *= fht_fir[size/2] * mul0;
	}
	//------------------------------------------------------------

private:
	void FHT::revbin_permute(double* data)
	{
		int rev_count = fht_revbin_counts[ldn];
		int* rev_table = fht_revbin_tables[ldn];
		for (int i = 0; i <rev_count; i++)
		{
			int x = rev_table[i*2];
			int r = rev_table[i*2+1];
			double vx = data[x];
			double vr = data[r];
			data[r] = vx;
			data[x] = vr;
		}
	}

	void FHT::r4step_addsub(double* data, int ldn) // real a[0..n-1] input,result
	{
		if(ldn==4)
		{
			for(int i=0;i<4;i++) 
			{
				double a = data[i];
				double b = data[i+4];
				double c = data[i+8];
				double d = data[i+12];

				double ac=a+c;
				double amc=a-c;
				double bd=b+d;
				double bmd=b-d;

				data[i]    =ac+bd;
				data[i+4]  =ac-bd;
				data[i+8]  =amc;
				data[i+12] =bmd;		
			}

			addsub(data[5],data[7],cos45);

			N21(data);
			N21(data+4);

			r4step_rotate8(data+8);
			N21(data+8);
			N21(data+8+4);
			return;
		}

		if(ldn==3)
		{
			r4step_addsub8(data);

			N21(data);
			N21(data+4);
			return;
		}

		if(ldn==2)
		{
			N21(data);
			return;
		}

		if(ldn==1)
		{
			addsub(data[0],data[1]);
			return;
		}

		//---------------------------
		int n = pow(2,ldn-1);
		int nh = n/2;

		for(int i=0;i<nh;i++) 
		{
			double a = data[i];
			double b = data[i+nh];
			double c = data[i+n];
			double d = data[i+n+nh];

			double ac=a+c;
			double amc=a-c;
			double bd=b+d;
			double bmd=b-d;

			data[i]=ac+bd;
			data[i+nh]=ac-bd;
			data[i+n]=amc;
			data[i+n+nh]=bmd;		
		}

		// recursion
		r4step_addsub(data, ldn-2);
		r4step_rotate(data+nh, ldn-2);
		r4step_rotate(data+n, ldn-1);
	}

	void FHT::r4step_rotate(double* data, int ldn) // real a[0..sz-1] input,result
	{
		int n = pow(2,ldn-1);
		int nh = n/2;

		addsub(data[0], data[n]);
		data[nh]*= sqrt2;
		data[n+nh]*= sqrt2;

		double* trig_table = fht_trig_tables[ldn];
		
		int index = 1;
		int step = n-2;
		while(step>0)
		{
			double a = data[index];
			double b = data[index+step];
			double c = data[index+n];
			double d = data[index+n+step];

			double ab = a+b;
			double amb = a-b;
			double cd = c+d;
			double cmd = c-d;

			double ss = trig_table[0];
			double cc = trig_table[1];
			data[index] =      ab * cc - cmd * ss;
			data[index+step] = ab * ss + cmd * cc;

			double ss3 = trig_table[2];
			double cc3 = trig_table[3];
			data[index+n] =      amb * cc3 + cd * ss3;
			data[index+n+step] = amb * ss3 - cd * cc3;

			trig_table+=4;
			index++;
			step-=2;
		}

		// recursion
		r4step_addsub(data, ldn-1);
		r4step_addsub(data+n, ldn-1);
	}

	//------------------------------------------------------------
	__inline void FHT::r4step_rotate8(double* data) // real a[0..sz-1] input,result
	{
		addsub(data[0], data[4]);
		data[2]*= sqrt2;
		data[6]*= sqrt2;

		double a = data[1];
		double b = data[3];
		double c = data[5];
		double d = data[7];

		double ab = a+b;
		double amb = a-b;
		double cd = c+d;
		double cmd = c-d;

		data[1] = ab * cos22 - cmd * sin22;
		data[3] = ab * sin22 + cmd * cos22;
		data[5] = cd * cos22 + amb * sin22;
		data[7] = amb * cos22 - cd * sin22;
	}

	__inline void FHT::r4step_addsub8(double* data)
	{
		addsub(data[0],data[4]);
		addsub(data[1],data[5]);
		addsub(data[2],data[6]);
		addsub(data[3],data[7]);
		addsub(data[5],data[7],cos45);
	}

	__inline void addsub(double& u, double& v)
	{
		double tempu = u;
		u = u + v;
		v = tempu - v;
	}

	__inline void addsub(double& u, double& v, const double scale)
	{
		double tempu = u;
		u = (u + v)*scale;
		v = (tempu - v)*scale;
	}

	__inline void rotate(double& u, double& v, const double c, const double s)
	{
		double tempu = u;
		u = u*c + v*s;
		v = tempu*s - v*c;
	}

	__inline void N21(double* data)
	{
		double v02=data[0]+data[2];
		double v13=data[1]+data[3];
		double v02m=data[0]-data[2];
		double v13m=data[1]-data[3];
		data[0]=v02+v13;
		data[1]=v02-v13;
		data[2]=v02m+v13m;
		data[3]=v02m-v13m;
	}

};