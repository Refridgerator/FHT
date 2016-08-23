#include <tchar.h>
#include <iostream>

#include "FHT.h"

void show(double* data, int size)
{
	for(int i=0;i<size;i++)
		std::cout << floor(data[i]*1000.+.5)/1000. << "\t";
	std::cout << std::endl << std::endl;
}

int _tmain(int argc, _TCHAR* argv[])
{
	double fir[]= {11,22,33,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0};
	double data[]={1,0,1,0,0,0,3,0,    0,0,0,0,0,0,0,0, 10,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0};

	int size = 32;
	FHT fht(size);
	
	std::cout << "Convolute test using FHT." << std::endl << std::endl;
	std::cout << "fir:" << std::endl;
	show(fir, size);
	std::cout << "source:" << std::endl;
	show(data, size);

	//	
	fht.transform(fir, false);
	fht.transform(data);		
	fht.convolute(data, fir);
	fht.back_transform(data);
	//
	std::cout << "result:" << std::endl;
	show(data, size);
	
	std::cin.get();
	return 0;
}
