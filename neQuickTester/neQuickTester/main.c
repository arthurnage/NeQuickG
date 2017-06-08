#include <stdio.h>
#include "neQuick.h"

int main()
{
    MODIP_st* modip = (MODIP_st*)malloc(sizeof(MODIP_st));
    GetModipFromFile(modip);
    CCIR_st* ccir = (CCIR_st*)malloc(sizeof(CCIR_st));
    GetCCIRFromFile(ccir);
    double pdKronrodTol[2] = {0.0, 0.0};
    int siMaxRecurse = 100;
    double    pdGssPosLLH[3] = {82.49, 297.66 - 360, 78.11};  /*// receiver position  (lat/lon/h)*/               //this one
    double    pdSatPosLLH[3] = {24.05, -158.03, 20275295.43};  /*// satellite position (lat/lon/h)*/              //this one
    double    pdSatPosLLH1[3] =  {54.29, 8.23, 20281546.18};
    int       siMonth = 1;         /*// month during which STEC value is required*/     //this one
    double    dUT = 0.0;             /*// time (UTC) at which STEC value is required*/    //this one
    int       siNumCoeff = 3;      /*// number of Az coeffs*/                               //3
    double    pdCoeff[3] = {236.831641, -0.39362878, 0.00402826613};      /*// Az coeffs*/                                 //these one
    double    dAzBase = 0.0;         /*// Az value at receiver locations*/
    NeQuickInputData_st *data = (NeQuickInputData_st*)malloc(sizeof(NeQuickInputData_st));
    data->pstModip = modip;
    data->pstCCIR = ccir;
    data->pdKronrodTol[0] = pdKronrodTol[0];
    data->pdKronrodTol[1] = pdKronrodTol[1];
    data->siMaxRecurse = siMaxRecurse;
    data->pdGssPosLLH[0] = pdGssPosLLH[0];
    data->pdGssPosLLH[1] = pdGssPosLLH[1];
    data->pdGssPosLLH[2] = pdGssPosLLH[2];
    data->pdSatPosLLH[0] = pdSatPosLLH1[0];
    data->pdSatPosLLH[1] = pdSatPosLLH1[1];
    data->pdSatPosLLH[2] = pdSatPosLLH1[2];
    data->siMonth = siMonth;
    data->dUT = dUT;
    data->siNumCoeff = siNumCoeff;
    data->pdCoeff[0] = pdCoeff[0];
    data->pdCoeff[1] = pdCoeff[1];
    data->pdCoeff[2] = pdCoeff[2];
    data->dAzBase = dAzBase;
    printf("keke\n");
    /* 3782752979130000.000000 */
    double res = NeQuick(data) * 40.3 / SQR(1575.42 * 1e6);
    printf("%f\n", res);
    
    return 0;
}
