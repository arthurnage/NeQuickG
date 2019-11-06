#include <assert.h>
#include <math.h>

#include "neQuick.h"

#include "ccir_f2.h"
#include "ccir_fm3.h"
#include "modip.h"
#include "validation_datasets.h"

void make_test() {
	const double eps = 3e-2;
	const int lat = 0;
	const int lon = 1;
	const int h = 2;
	MODIP_st modip = { 0 };
	modip.pdModip = &gstModip;

	CCIR_st ccir = {0};
	ccir.pdF2 = &gF2;
	ccir.pdM3000 = &gFm3;
	
	NeQuickInputData_st data = {0};
	data.siMaxRecurse = 100;
	data.siNumCoeff = 3;
	data.pstModip = &modip;
	data.pstCCIR = &ccir;
	int datasetId = 0;
	int rowId = 0;
	double ans = 0.0;
	double err = 0.0;
	for (datasetId = 0; datasetId < 3; datasetId++) {
		data.pdCoeff[0] = datasets[datasetId].broadcast_coefficients[0];
		data.pdCoeff[1] = datasets[datasetId].broadcast_coefficients[1];
		data.pdCoeff[2] = datasets[datasetId].broadcast_coefficients[2];
		for (rowId = 0; rowId < 36; rowId++) {
			data.siMonth = datasets[datasetId].inst[rowId].month;
			data.dUT = datasets[datasetId].inst[rowId].utc;
			
			data.pdGssPosLLH[lat] = datasets[datasetId].inst[rowId].stlat;
			data.pdGssPosLLH[lon] = datasets[datasetId].inst[rowId].stlon;
			data.pdGssPosLLH[h] = datasets[datasetId].inst[rowId].stheight;
			
			data.pdSatPosLLH[lat] = datasets[datasetId].inst[rowId].satlat;
			data.pdSatPosLLH[lon] = datasets[datasetId].inst[rowId].satlon;
			data.pdSatPosLLH[h] = datasets[datasetId].inst[rowId].satheight;
			ans = NeQuick(&data) / 1e16;
			err = fabs(ans - datasets[datasetId].inst[rowId].stec) / datasets[datasetId].inst[rowId].stec;
			assert(err  < eps);
		}
	}
}

int main() {
	make_test();
    return 0;
}
