#ifndef __RHEHpt_H__
#define __RHEHpt_H__
#include <cmath>
#include <string>
#include <memory>
#include "LOBaur.h"
#include "LHAPDF/LHAPDF.h"

namespace RHEHpt
{
    constexpr long double Gf = 0.00001166364; // Fermi constant in GeV^-2
    constexpr long double gev2_to_pb = 389379304.; // GeV^-2 to pb conversion factor == (hc)^2 
   	typedef const std::vector< double >& vector_ref;
   	typedef const std::vector< std::complex< double > >& cvector_ref;
   	
class RHEHpt
{
    public:
    	LOBaur Exact_FO_ME;
        RHEHpt(double CME, const std::string& PDFname = "NNPDF_30_nnlo_as_118", double MH = 125.09,
        	   double MT = 173.3, double MB=4.18, bool verbose = false);
		double core_hpt_infinite(double N) const;
        double hpt_infinite(double pt,double N) const;
        double hpt_finite(double pt,double N) const;
        double ptd_FO(double pt);
        double C(unsigned i,unsigned j,double pt) const; // coefficient of the series expansion of the resummed result

		std::vector< double > Integral_coeffs(unsigned order,double xp,bool heavy_quark = false) const;
		std::vector< double > Xp_prefactor_coeffs(unsigned order,double xp) const;
		
		/* 	Given a vector produced by pt_distr_series_terms, pt_distr_series returns the series up to a specific order. 
			In this way if you compute the NNNLO distr and at some point you need the NNLO distr you have to recompute nothing */
		std::vector< double > pt_distr_series_terms(unsigned order, double xp, double N, bool heavy_quark = false) const;

        template<class T> T _an_dim(T as_N) const{ // LO expansion of the LLx anomalous dimension (gamma_s)
			const T CA = 3.;
			const T lambda = CA * as_N ;
			return lambda/M_PI;
		}

		template<class T> std::vector< T > pt_distr_series_terms(T N, vector_ref integral_coeff_list,
																	vector_ref xp_coeff_list, unsigned order) const {
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

			T M = _an_dim(_as/N);

			std::vector<T> terms;
			terms.push_back(M*integral_coeff_list[1]*xp_coeff_list[0]); // integral_coeff_list[i] * xp_coeff_list[j] --> M^(i+j)
			if(order > 1){
				for(unsigned i = 2 ; i < order + 1 ; ++i){
		
					// get couples of int (x,y) s.t. x + y = i
					std::vector< std::array<unsigned,2> > indices = symmetric_partition_of(i);
					std::vector< T > tmp( indices.size() );
		
					// compute all coefficients of the i-th order
					std::transform(indices.cbegin(),indices.cend(),tmp.begin(),[&](const std::array<unsigned,2>& index_pair){
						//N.B. by definition of partition_of(), index_pair[0] >= index_pair[1]
						return std::pow(M,i)*(  integral_coeff_list[index_pair[0]]*xp_coeff_list[ index_pair[1] ]) ; 
			
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
			std::cout << "useful_number_of_terms = " << useful_number_of_terms << std::endl;
			for (int i = 0 ; i < useful_number_of_terms ; ++i){
				std::cout << terms[i] << " " ;
			}
			std::cout << std::endl;	
			const double sigma0 = _as * _as * Gf * std::sqrt(2.) / (576. * M_PI);
			return sigma0 * std::accumulate(terms.cbegin(),terms.cbegin() + useful_number_of_terms,0.);
		}


		double pt_distr_series(unsigned order, double xp, double N, bool heavy_quark = false) const;

		inline  double get_yt() const{ return (_mT*_mT)/(_mH*_mH); }
		inline  double get_alphas() const{ return _as; }
		inline double get_yb() const{ return (_mB*_mB)/(_mH*_mH); }
		void set_mt(double mt){ _mT = mt;}
		void set_mb(double mb){ _mB = mb;}
		void set_choice(int CHOICE){
		  if ((CHOICE!=1)||(CHOICE!=2)||(CHOICE!=3)||(CHOICE!=4)){
		    std::cout << "Error invalid choice; set a number from 1 to 4" << std::endl;  
		  }
		  else {
		    choice = CHOICE;
		  }
		}

        virtual ~RHEHpt();
        
        void set_scale(double Q){
            _CME = Q;
            _s = Q*Q;
            //_as = alphas(Q);
            Exact_FO_ME.SetCME(Q);
        }
        
    private:
		/* utility function */
    	std::vector< std::array<unsigned,2> > partition_of(unsigned n) const;
    	std::vector< std::array<unsigned,2> > symmetric_partition_of(unsigned n) const;
    	
        double R(double as_N) const;

        double _CME; // center of mass energy
        double _s;  // hadronic mandelstam variable
        double _mH; // Higgs mass
        double _mT; // Top mass
        double _mB; // Bottom mass
        double _as; // Current value of alphas, default is alphas(MH)
        bool _verbose;
	unsigned int choice=4;
        std::shared_ptr<LHAPDF::PDF> _PDF;
    //    Luminosity _lumi;   // That provides luminosity functions
    //    std::unique_ptr<LHAPDF::PDFSet> _PDFset;
    //    std::vector< std::unique_ptr< LHAPDF::PDFSet > > _PDFmembers;
};



}

#endif /* __RHEHpt_H__ */
