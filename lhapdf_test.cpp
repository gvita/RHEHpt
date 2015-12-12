#include "LHAPDF/LHAPDF.h"
#include <memory>
#include <vector>
#include <iostream>

int main(){
	std::vector< LHAPDF::PDF* > ps = LHAPDF::mkPDFs("NNPDF30_nnlo_as_0118");
	for (auto pdf: ps){
		std::cout << pdf -> xfxQ(21, 1e-4, 100.) << std::endl;
	}
	return 0;
}

