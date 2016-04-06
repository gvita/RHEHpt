#ifndef __INTERPOLATOR_H__
#define __INTERPOLATOR_H__

#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

namespace RHEHpt{

class Interpolator
{
	public:
		Interpolator (unsigned choice = 0);
		virtual ~Interpolator ();
		
		double operator()(unsigned order, double xp) const;

		void set_choice(unsigned choice){ _choice = choice; }
		
	private:
		unsigned _choice; 
		gsl_interp_accel* _gsl_acc;
		
		/* matrix of gsl interpolators one row per _choice and one column per order */
		std::vector < std::vector < gsl_spline *> > _interpolators; 
	
};


} /* close RHEHpt namespace* */

#endif /* __INTERPOLATOR_H__ */

