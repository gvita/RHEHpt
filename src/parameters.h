#include <functional>

class _parameters{
	public:
		typedef std::function < long double(long double, long double, long double, long double)> HardF1;
		typedef std::function < long double(long double, long double, long double)> HardF2;
		_parameters(){};
		
		void set(unsigned J,unsigned K,long double XP, HardF1 F_0F, HardF1 D_0F, HardF2 F_INF){
//			std::cout << "Setting _parameters" << std::endl;
			j = J;
			k = K;
			xp = XP;
			F_0f = F_0F;
			D_0f = D_0F;
			F_inf = F_INF;
//			std::cout << "_parameters setted" << std::endl;
//			std::cout << "test:" << F_0f(1,1,1,1) << std::endl;
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
