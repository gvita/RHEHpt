
#include "NLO_PL.h"
#include "Luminosity.h"
#include <functional>

class _parameters{
	public:
		typedef std::function < long double(long double, long double, long double, long double)> HardF1;
		typedef std::function < long double(long double, long double, long double)> HardF2;
		_parameters(){};
		
		void set(unsigned J,unsigned K,long double XP, HardF1 F_0F, HardF1 D_0F, HardF2 F_INF){
			j = J;
			k = K;
			xp = XP;
			F_0f = F_0F;
			D_0f = D_0F;
			F_inf = F_INF;
		}
		
		void switch_jk(){
			unsigned tmp = j;
			j = k;
			k = tmp;
		}
		
		unsigned j;
		unsigned k;
		long double xp;
		HardF1 F_0f;
		HardF1 D_0f;
		HardF2 F_inf;
};

class _par_part{
public:
  _par_part(NLOPL& Point, long double xp):_Pointlike_int(Point),_xp_int(xp){}
  NLOPL& _Pointlike_int;
  long double _xp_int;
  
};

class _par_hadro{
	public:
		_par_hadro(Luminosity& Lum, LOBaur& Fin, NLOPL& Point, long double xp, long double tau, long double mH, unsigned int channel): _Lumi_int(Lum),_Finite_int(Fin),_Pointlike_int(Point),_xp_int(xp),_tau_int(tau),_mH_int(mH),_channel_int(channel){};
		Luminosity& _Lumi_int;
		LOBaur& _Finite_int;
		NLOPL& _Pointlike_int;
		long double _xp_int;
		long double _tau_int;
		long double _mH_int;
		unsigned int _channel_int;
};

class _par_expansion{
	public:
		typedef std::function < long double(double)> Partonic_Distr;
		_par_expansion(Luminosity& Lum, Partonic_Distr part_f, long double xp, long double tau, long double mH):
						_Lumi_int(Lum), ds_xp_part(part_f),_xp_int(xp),_tau_int(tau),_mH_int(mH){};
		Luminosity& _Lumi_int;
		Partonic_Distr ds_xp_part;
		long double _xp_int;
		long double _tau_int;
		long double _mH_int;	
};
