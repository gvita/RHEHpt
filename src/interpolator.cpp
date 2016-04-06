#include "interpolator.h"
#include <string>
#include <fstream>
#include <iostream>

namespace RHEHpt{

Interpolator::Interpolator(unsigned choice):_choice(choice)
{
//	std::cout << "creating interpolator..." << std::endl;
	std::vector<std::string> filenames{"src/top.dat","src/bottom.dat","src/interference.dat"};
	std::size_t n_choices = filenames.size();
	
	_gsl_acc = gsl_interp_accel_alloc();
	
	for (unsigned int i = 0; i < n_choices; ++i){
		std::vector<double> xp;
		std::vector<double> nlo;
		std::vector<double> nnlo;
		double tmp;
		std::vector< gsl_spline* > interpolator;
		
		/* read data from file */
		std::ifstream ing(filenames[i].c_str());
//		std::cout << "start reading " << std::endl;
		ing >> tmp;
		while(ing){
			xp.push_back(tmp);

			ing >> tmp;
			nlo.push_back(tmp);
		
			ing >> tmp;
			nnlo.push_back(tmp);
			
			ing >> tmp;
		}
//		std::cout << "end reading. " << xp.size() << std::endl;		
		/*** create nlo interpolator ***/
		interpolator.push_back( gsl_spline_alloc( gsl_interp_akima, xp.size() ) );
		gsl_spline_init( interpolator[0], &xp[0], &nlo[0], xp.size() );
//		std::cout << "first create. " << xp.size() << std::endl;
		/*** create nnlo interpolator ***/
		interpolator.push_back( gsl_spline_alloc( gsl_interp_akima, xp.size() ) );
		gsl_spline_init( interpolator[1], &xp[0], &nnlo[0], xp.size() );

//		std::cout << "created interpolator " << i << std::endl;	

		/*** assign it to class member ***/
		_interpolators.push_back(interpolator);
	}
}

Interpolator::~Interpolator(){
	for (auto choice : _interpolators) // foreach choice
		for ( auto order : choice) // foreach order
			gsl_spline_free(order);
	gsl_interp_accel_free(_gsl_acc);
}

double Interpolator::operator()(unsigned order,double xp) const
{
	/* There is a mismatch of 2 between the interpolator order and order since order=0 and order=1 are special cases not handled by interpolation. Analgously there is a mismatch of 1 between the choices since _choice = 0 is a special case.*/
//	std::cout << "Interpolator(" << _choice << "," << order << "," << xp << ")" << std::endl;
	if(_choice == 0){ // tot is the sum of the three
		return gsl_spline_eval( _interpolators[0][order-2], xp, _gsl_acc ) + 
		       gsl_spline_eval( _interpolators[1][order-2], xp, _gsl_acc ) + 
		       gsl_spline_eval( _interpolators[2][order-2], xp, _gsl_acc );
	}
	else
		return gsl_spline_eval( _interpolators[_choice-1][order-2], xp, _gsl_acc );
}

}
