#ifndef __RHEHpt_H__
#define __RHEHpt_H__

#include <cmath>
#include <string>
#include <memory>
#include "LOBaur.h"
#include "Luminosity.h"
#include "NLO_PL.h"
#include "LHAPDF/LHAPDF.h"

namespace RHEHpt
{
    constexpr long double Gf = 0.00001166364; // Fermi constant in GeV^-2
    constexpr long double gev2_to_pb = 389379304.; // GeV^-2 to pb conversion factor == (hc)^2 
   	
class RHEHpt
{
    public:
	   	typedef const std::vector< double >& vector_ref;
	   	typedef const std::vector< std::complex< double > >& cvector_ref;
		typedef LHAPDF::PDF* PDF_ptr;
		typedef std::function< std::complex<long double>(std::complex<long double>)> c_function;

        RHEHpt(double CME, const std::string& PDFname = "NNPDF30_nnlo_as_0118", double MH = 125.09,
        	   double MT = 173.3, double MB = 4.18, double MUR=125.09, double MUF=125.09, unsigned choice = 1, bool verbose = false);
		double core_hpt_infinite(double N) const;
        double hpt_infinite(double pt,double N) const;
        double hpt_finite(double pt,double N) const;
	
	
        double C(unsigned i,unsigned j,double pt) const; // coefficient of the series expansion of the resummed result

		std::vector< double > Integral_coeffs(unsigned order,double xp,bool heavy_quark = false) const;
		std::vector< double > Xp_prefactor_coeffs(unsigned order,double xp) const;
		
		long double M(long double x,unsigned i) const;

        template<class T> T _an_dim(T as_N) const{ // LO expansion of the LLx anomalous dimension (gamma_s)
			const T CA = 3.;
			const T lambda = CA * as_N ;
			return lambda/M_PI;
		}

		/* 	Given a vector produced by pt_distr_series_terms, pt_distr_series returns the series up to a specific order. 
			In this way if you compute the NNNLO distr and at some point you need the NNLO distr you have to recompute nothing */
		std::vector< double > pt_distr_series_terms(unsigned order, double xp, double N,
													bool heavy_quark = false, bool Nspace = true) const;

		template<class T> std::vector< T > pt_distr_series_terms(T N, vector_ref integral_coeff_list, vector_ref xp_coeff_list,
																 unsigned order, bool Nspace = true) const {
		/* 	This version returns a vector of pt_distr_series_terms, given a vector produced by Integral_coeffs and another
			one produced by Xp_prefactor_coeffs. In this way the N dipendence is factorized from the xp	depend terms,
			which is crucial when computing the inverse mellin transform of the pt distribution	without recomputing each
			coefficient for each N.

			IMPORTANT NOTE: This is formally correct up to order 5 then you need to expand also _an_dim

			Note:
			- xp_coeff_list[i] = 1/xp * ((2log(xp))^i)/i! (the coeff expansion of xp^(2M - 1) = sum_i M^i xp_coeff_list[i])
			- Integral_coeffs[i] = sum_{j,k} C(j,k) s.t. j + k = i and j != k 

			In this way:
				integral_coeff_list[i] * xp_coeff_list[j] --> M^(i+j)
		*/

			std::vector<T> terms;
			T M_N, x;
			if ( Nspace ){ 
				M_N = _an_dim(_as/N);
				terms.push_back(M_N*integral_coeff_list[1]*xp_coeff_list[0]); // integral_coeff_list[i] * xp_coeff_list[j] --> M^(i+j)
			}
			else{
				x = N;
				terms.push_back( M(x,1)*integral_coeff_list[1]*xp_coeff_list[0]); // M(x,n) = inverse mellin of M^n at x
			}
			if(order > 1){
				for(unsigned i = 2 ; i < order + 1 ; ++i){
		
					// get couples of int (x,y) s.t. x + y = i
					std::vector< std::array<unsigned,2> > indices = symmetric_partition_of(i);
					std::vector< T > tmp( indices.size() );
		
					// compute all coefficients of the i-th order
					std::transform(indices.cbegin(),indices.cend(),tmp.begin(),[&](const std::array<unsigned,2>& index_pair) -> T{
						//N.B. by definition of partition_of(), index_pair[0] >= index_pair[1]
						if ( Nspace )
							return std::pow(M_N,i)*(  integral_coeff_list[index_pair[0]]*xp_coeff_list[ index_pair[1] ]) ;
						else
							return M(x,i)*(  integral_coeff_list[index_pair[0]]*xp_coeff_list[ index_pair[1] ]) ;
					});
					terms.insert(terms.end(),tmp.begin(),tmp.end());
				}
			}
			if (_verbose){
				std::cout << "terms " << std::endl;
				for (auto it : terms){
					std::cout << it << " " ;
				}
				std::cout << std::endl;
			}
			return terms;

		}

		
		template<class T> T pt_distr_series(const std::vector< T >& terms,unsigned order) const {
			unsigned useful_number_of_terms = ( order*(order + 3) ) / 2 -1; // sum_{n=1}^order (n + 1)
			if(_verbose){
				std::cout << "useful_number_of_terms = " << useful_number_of_terms << std::endl;
				for (unsigned int i = 0 ; i < useful_number_of_terms ; ++i){
					std::cout << terms[i] << " " ;
				}
				std::cout << std::endl;	
			}
			T zero = static_cast<T> (0.);
			return SIGMA0() * std::accumulate(terms.cbegin(),terms.cbegin() + useful_number_of_terms,zero);
		}


		double pt_distr_series(unsigned order, double xp, double N, bool heavy_quark = false,bool Nspace = true) const;

		long double pt_distr_hadro(long double pt, unsigned int order = 1, bool heavyquark = false);
		std::vector<double> pt_distr_hadro(const std::vector< double >& ptgrid, unsigned int order = 1, bool heavyquark = false);	
		inline double get_yt() const{ return (_mT*_mT)/(_mH*_mH); }
		inline double get_yb() const{ return (_mB*_mB)/(_mH*_mH); }
		inline double get_alphas() const{ return _as; }
		inline PDF_ptr get_PDF() const { return _PDF; }
		inline double get_scale() const{ return _CME; }
		inline double get_mH() const{ return _mH; }
		inline double get_xp(double pt) const { return std::pow(pt/_mH,2.);}
		inline double get_x() const{ return std::pow(_mH/_CME,2.); }
		inline double get_tau() const {return std::pow(_mH/_CME,2.); }

		void set_mh(double mh){ _mH = mh;}		
		void set_mt(double mt){ _mT = mt;}
		void set_mb(double mb){ _mB = mb;}
		void set_mur(double mur){ 
		  _muR=mur;
		  _as=_PDF->alphasQ(mur);
		  Exact_FO_fullmass.SetAS(_as);
		  Exact_FO_PL.SetAS(_as);
		}
	
		void set_muf(double muf){
		  _muF=muf;
		  _Lum.Cheb_Lum(muf);
		}

        void set_CME_collider(double Q){
            _CME = Q;
            _s = Q*Q;
        }
		
		void set_choice(unsigned int CHOICE){
			if(CHOICE > 3 )
		    	std::cout << "Error invalid choice; set a number from 0 to 3" << std::endl;  
			else{
				_choice = CHOICE;
				Exact_FO_fullmass.setchoice(_choice);
			}
		}

		// dimensionless LO cross section of the pointlike inclusive production. It plays the role of a normalization factor.
		inline long double SIGMA0() const { return (_as*_as) * std::sqrt(2)*Gf/(576.*M_PIl); } 

        virtual ~RHEHpt();
                
        std::complex<long double> Lum_gg_N(std::complex<long double> N){
        	return _Lum.CLum_gg_N(N);
        }
        
        long double Lum_gg_x(long double x){
        	return _Lum.CLum_gg_x(x);
        }
        
        /* utility class _parameters. This has to be public because it has to be seen by the Cuba integrand. */
   		#include "parameters.h"
        
		//Fixed Order part
		long double sigma_part(long double CME_part, long double pt, unsigned int order=1, bool heavyquark=true);
		long double sigma_hadro_FO_fullmass(long double pt);
		std::vector<double> sigma_hadro_FO_pointlike (std::vector<double>& ptgrid, unsigned int order = 1, int channel = 1);
    	LOBaur Exact_FO_fullmass;
		NLOPL  Exact_FO_PL;
        
        
    private:
		/* utility function */
    	std::vector< std::array<unsigned,2> > partition_of(unsigned n) const;
    	std::vector< std::array<unsigned,2> > symmetric_partition_of(unsigned n) const;

        double R(double as_N) const;

        double _CME; // center of mass energy of the collider
        double _s;  // hadronic mandelstam variable
        double _mH; // Higgs mass
        double _mT; // Top mass
        double _mB; // Bottom mass
        double _as; // Current value of alphas, default is alphas(MH)
        double _muR;
		double _muF;
        bool _verbose;
		unsigned int _choice;
		LHAPDF::PDF* _PDF;
		std::function < long double(long double, long double, long double, long double , long double, long double)> F_0f;
		std::function < long double(long double, long double, long double, long double , long double, long double)> D_0f;
		Luminosity _Lum;
};

extern "C" {
	void hqt_(double* sroot, double* amh, double* resummscale,  double* mur, double* muf, int* flag1,
			  const char* name, int* namelen, int* mem, double* ptgrid, double* result,int* gridsize, int* channel2);
}


} //end NAMESPACE RHEHpt

#endif /* __RHEHpt_H__ */
