/*// NeQuickG Header File*/

#ifndef NEQUICK_H
#define NEQUICK_H

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define RADIUS_OF_EARTH      6371.2

#ifndef PI
#define PI                   3.1415926536
#endif

#define DEGTORAD             0.01745329251
#define RADTODEG             57.2957795131

#ifndef SQR(x)
#define SQR(x) ( (x) * (x) )
#endif

typedef struct {
    
    double    pdModip[39][39];       /*// grid of modified dip latitude values*/
    
} MODIP_st;

typedef struct {
    
    double    pdF2[12][2][76][13];   /*// CCIR coeffs for computing FoF2*/
    double    pdM3000[12][2][49][9]; /*// CCIR coeffs for computing M(3000)F2*/
    
} CCIR_st;

typedef struct {
    
    double    dLat;        /*// latitude of point*/
    double    dLng;        /*// longitude of point*/
    double    dH;          /*// height of point*/
    double    dR;          /*// radius of point*/
    double    dS;          /*// distance of point to ray parigee*/
    double    dSinLat;     /*// sin of latitude*/
    double    dCosLat;     /*// cos of latitude*/
    
} SPoint_st;


typedef struct {
    
    MODIP_st* pstModip;        /*// structure containing grid of*/
    /*// modified dip latitude values*/                //this one
    CCIR_st*  pstCCIR;         /*// structure containing CCIR coeffs for
                                // computing FoF2 amd M(3000)F2*/
    double    pdKronrodTol[2]; /*// tolerances for Kronrod integration*/
    int       siMaxRecurse;    /*// max level of recursion allowed in               //?? какой он? где он? что он?
                                // Kronrod integration*/
    double    pdGssPosLLH[3];  /*// receiver position  (lat/lon/h)*/               //this one
    double    pdSatPosLLH[3];  /*// satellite position (lat/lon/h)*/                //this one
    int       siMonth;         /*// month during which STEC value is required*/     //this one
    double    dUT;             /*// time (UTC) at which STEC value is required*/    //this one
    int       siNumCoeff;      /*// number of Az coeffs*/                               //3
    double    pdCoeff[3];      /*// Az coeffs*/                                 //these one
    double    dAzBase;         /*// Az value at receiver locations*/
    
} NeQuickInputData_st;

typedef struct {
    
    double    dLat;            /*// latitude of point*/
    int       siMonth;         /*// month during which current STEC value
                                // has been computed*/
    double    dR12;            /*// current R12 index*/
    double    pdF0F2[988];     /*// interpolated coeffs for computing FoF2 for
                                // current month and R12 conditions*/
    /* 975 instead of 988 -bug*/
    double    pdM3000F2[441];  /*// interpolated coeffs for computing M(3000)F2 for
                                // current month and R12 conditions*/
    double    dUT;             /*// time (UTC) at which current STEC value
                                // has been computed*/
    double    pdLegCoeffs_F0[76];/*// spherical Legendre coeffs for calculating FoF2 for
                                  // current month and R12 conditions*/
    double    pdLegCoeffs_M3000[49];/*// spherical Legendre coeffs for calculating
                                     // M(3000)F2 for current month and R12 conditions*/
    
} CurrentCCIR_st;

typedef struct {
    
    double     pdAmp[3];       /*// Epstein amplitude param*/
    double     pdPeakHeight[3];/*// Epstein peak height param*/
    double     pdBotThick[3];  /*// Epstein bottom half-layer thickness param*/
    double     pdTopThick[3];  /*// Epstein top half-layer thickness param*/
    double     dM3000;         /*// current M(3000)F2 value*/
    double     pdF0[3];        /*// current Fo for the F2, F1, and E layers respectively*/
    
} LayersProperties_st;

typedef struct {
    
    SPoint_st  stP1;           /*// info for point 1 (receiver)*/
    SPoint_st  stP2;           /*// info for point 2 (satellite)*/
    SPoint_st  stRay;          /*// info for ray*/
    SPoint_st  stPactual;      /*// info for current integration point*/
    double     dZeta;          /*// zenith angle of point 2 seen from point 1*/
    double     dSinDelta;      /*// sin of angle of declination of sun*/
    double     dCosDelta;      /*// cos of angle of declination of sun*/
    double     dSinSig;        /*// sin of ray azimuth*/
    double     dCosSig;        /*// cos of ray azimuth*/
    
} Geometry_st;

typedef struct {
    
    NeQuickInputData_st*   pstNeQuickInputData; /*// input data for NeQuick func*/
    Geometry_st*           pstGeom;             /*// geometry data for ray*/
    CurrentCCIR_st*        pstCurrCCIR;         /*// FoF2 and M(3000)F2 info for
                                                 // current month and R12*/
    double                 dTolerance;          /*// tolerance for Kronrod integration*/
    int                    bVert;               /*// flag indicating whether ray
                                                 // is vertical or not*/
} IntegrateData_st;


void ConvertPosFromInputDataToGeom( NeQuickInputData_st* pstNeQuickInputData, Geometry_st* pstGeom );
void CalcSunDeclinationAngle( int mth, double ut, double* SinDelta, double* CosDelta );
extern void GetModipFromFile( MODIP_st* pstModip );
extern void GetCCIRFromFile( CCIR_st* pstCCIR );


extern double NeQuick( NeQuickInputData_st* pstNeQuickInputData );
int NeqCheckInputs( NeQuickInputData_st* pstNeQuickInputData );
double DoTECIntegration( IntegrateData_st* pstIntegrateData, int bVert, SPoint_st stP0,
                        double* pdNmax, LayersProperties_st* pstLayers);


double NeqCalcModip(double dLat, double dLng, MODIP_st* pstModip);
double NeqInterpolate(double* pdZ, double dDeltaX);
double NeqModipToAz(double dModip, int siNumCoeff, double* pdCoeff);


int NeqGetRayProperties( SPoint_st* pstP1, SPoint_st* pstP2, SPoint_st* pstRay,
                        double* pdZeta, double* pdSinSig, double* pdCosSig );
void NeqCalcRayProperties1( SPoint_st* pstP1, SPoint_st* pstP2, SPoint_st* pstRay,
                           double* pdZeta );
void NeqCalcRayProperties2( SPoint_st* pstP1, SPoint_st* pstP2, SPoint_st* pstRay,
                           double* pdSinSig, double* pdCosSig );


double NeqIntegrate( IntegrateData_st* pstIntegrateData, double dH1, double dH2,
                    int siCurrentLevel, SPoint_st* pstPactual, double* pdNmax,
                    LayersProperties_st* pstLayers, double AbsTol );


double NeqGetNeOnVertRay( double dH, LayersProperties_st* pstLayers, double* pdNmax );


double NeqGetNeOnSlantRay( double dS, IntegrateData_st* pstIntegrateData, double* pdNmax,
                          SPoint_st* pstPactual, LayersProperties_st* pstLayers);
void NeqCalcLLHOnRay( double dS, SPoint_st* pstRay, SPoint_st* pstP1,
                     double dSinSig, double dCosSig, SPoint_st* pstPactual );


void NeqCalcEpstParams( NeQuickInputData_st* pstNeQuickInputData, SPoint_st* pstPactual,
                       double dSinDelta, double dCosDelta, double* pdNmax,
                       LayersProperties_st* pstLayers, CurrentCCIR_st* pstCurrCCIR );
void NeqCalcSphLegCoeffs( double dUT, CurrentCCIR_st* pstCurrCCIR );
double NeqGetF2FreqFromCCIR( double dCosLat, double dLng, double* pdLegCoeffs,
                            double* pdSinModipToN, int siMode );
double NeqCriticalFreqToNe( double dF0 );
double NeqCalcF2PeakHeight( double dM3000, double dF0E, double dF0F2);
double NeqCalcEpstein( double dNmax, double dHmax, double dB, double dH );


double NeqCalcTopsideNe( double* pdNmax, LayersProperties_st* pstLayers, double dH );


double NeqCalcBottomsideNe( double dHH, LayersProperties_st* pstLayers );


double NeqJoin( double dF1, double dF2, double dAlpha, double dX );
double NeqClipExp( double dPower );

#endif
