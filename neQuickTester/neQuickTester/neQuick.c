/*// NeQuickG Source File*/

#include "neQuick.h"

#include <assert.h>
#include <math.h>

/*
 /////////////////////////////////////////////////////////////////////////////////
 // Auxiliary Funcs
 /////////////////////////////////////////////////////////////////////////////////
 */
void ConvertPosFromInputDataToGeom(NeQuickInputData_st* pstNeQuickInputData, Geometry_st* pstGeom)
{
	pstGeom->stP1.dLat = pstNeQuickInputData->pdGssPosLLH[0];
	pstGeom->stP1.dLng = pstNeQuickInputData->pdGssPosLLH[1];
	pstGeom->stP1.dH = pstNeQuickInputData->pdGssPosLLH[2] / 1000;
	pstGeom->stP1.dSinLat = sin(pstGeom->stP1.dLat * DEGTORAD); /*degtorad correction -bug*/
	pstGeom->stP1.dCosLat = cos(pstGeom->stP1.dLat * DEGTORAD);

	pstGeom->stP2.dLat = pstNeQuickInputData->pdSatPosLLH[0];
	pstGeom->stP2.dLng = pstNeQuickInputData->pdSatPosLLH[1];
	pstGeom->stP2.dH = pstNeQuickInputData->pdSatPosLLH[2] / 1000;
	pstGeom->stP2.dSinLat = sin(pstGeom->stP2.dLat * DEGTORAD); /*degtorad correction -bug*/
	pstGeom->stP2.dCosLat = cos(pstGeom->stP2.dLat * DEGTORAD);
}

void CalcSunDeclinationAngle(int mth, double ut, double* SinDelta, double* CosDelta)
{
	double dy = 30.5 * mth - 15;
	double t = dy + (double)(18 - ut) / 24;
	double am = (0.9856 * t - 3.289) * DEGTORAD;
	double al = am + (1.916 * sin(am) + 0.02 * sin(2 * am) + 282.634) * DEGTORAD;

	*SinDelta = 0.39782 * sin(al);
	*CosDelta = sqrt(1 - SQR(*SinDelta));
}

const int F2 = 0;
const int F1 = 1;
const int E = 2;

extern int numOfIter = 0;

///*
// /////////////////////////////////////////////////////////////////////////////////
// // NeQuick.c
// /////////////////////////////////////////////////////////////////////////////////
// */
/*checked*/
extern double NeQuick(NeQuickInputData_st* pstNeQuickInputData)
{
	int        mth = pstNeQuickInputData->siMonth;
	double     ut = pstNeQuickInputData->dUT;
	double*    ai = pstNeQuickInputData->pdCoeff;
	int        ncoeff = pstNeQuickInputData->siNumCoeff;

	double ModipRx, Az_U;
	/*double amrad, sdelta, cdelta;*/
	double dTEC = 0;
	double dNmax = 0;

	Geometry_st stGeom = { 0 };
	IntegrateData_st stIntegrateData = { 0 };
	CurrentCCIR_st stCurrCCIR = { 0 };
	LayersProperties_st stLayers = { 0 };

	SPoint_st* stP1 = &(stGeom.stP1);
	SPoint_st* stP2 = &(stGeom.stP2);
	SPoint_st  stP0;
	SPoint_st* stRay = &(stGeom.stRay);
	SPoint_st* stPactual = &(stGeom.stPactual);

	double* Zeta = &(stGeom.dZeta);
	double* SinSig = &(stGeom.dSinSig);
	double* CosSig = &(stGeom.dCosSig);

	ConvertPosFromInputDataToGeom(pstNeQuickInputData, &stGeom);
	CalcSunDeclinationAngle(mth, ut, &(stGeom.dSinDelta), &(stGeom.dCosDelta));

	stIntegrateData.pstNeQuickInputData = pstNeQuickInputData;
	stIntegrateData.pstGeom = &stGeom;
	stIntegrateData.pstCurrCCIR = &stCurrCCIR;                                 /*// not correct ??*/
	stIntegrateData.bVert = 0;

	/*// abs tolerance of integration below 1000 km*/
	pstNeQuickInputData->pdKronrodTol[0] = 0.001;
	/*// abs tolerance of integration above 1000 km/*/
	pstNeQuickInputData->pdKronrodTol[1] = 0.01;

	if (NeqCheckInputs(pstNeQuickInputData)) {

		//ModipRx = NeqCalcModip(81.3408488001651,   319.9066461160811, pstNeQuickInputData->pstModip);
		ModipRx = NeqCalcModip(stP1->dLat, stP1->dLng, pstNeQuickInputData->pstModip);
		Az_U = NeqModipToAz(ModipRx, ncoeff, ai);

		pstNeQuickInputData->dAzBase = Az_U;

		/*  why? */
		if (ncoeff == 3)
			if ((ai[0] == 0) && (ai[1] == 0) && (ai[2] == 0))
				Az_U = 63.7;
		if (ncoeff == 2)
			if ((ai[0] == 0) && (ai[1] == 0))
				Az_U = 63.7;
		if (ncoeff == 1)
			if (ai[0] == 0)
				Az_U = 63.7;

		if (Az_U < 0)
			Az_U = 0;

		if (Az_U > 400)
			Az_U = 400;

		if (!NeqGetRayProperties(stP1, stP2, stRay, Zeta, SinSig, CosSig)) /* returns true instead of false -bug*/
		{
			stP1->dS = sqrt(SQR(stP1->dR) - SQR(stRay->dR));
			stP2->dS = sqrt(SQR(stP2->dR) - SQR(stRay->dR));

			if (stP1->dH < 0)
				stP0.dH = 0;
			else
				stP0.dH = stP1->dH;

			stP0.dR = stP0.dH + RADIUS_OF_EARTH;
			stP0.dS = sqrt(SQR(stP0.dR) - SQR(stRay->dR)); /*it was a plus -bug*/

			stPactual->dH = stP2->dH; /*P2*/
			stPactual->dLat = stP1->dLat;
			stPactual->dLng = stP1->dLng;

			if (stRay->dR < 0.1) {
				NeqCalcEpstParams(pstNeQuickInputData, stPactual,
					stGeom.dSinDelta, stGeom.dCosDelta,
					&dNmax, &stLayers, &stCurrCCIR);
				stIntegrateData.bVert = 1;
			}
			/*update pstIntegrateData*/

			dTEC = DoTECIntegration(&stIntegrateData, stIntegrateData.bVert, stP0, &dNmax, &stLayers);
			return 1000 * dTEC;
		}
		else
			/*// todo: return error*/
			return 0;
	}
	else
		/*// todo: return error*/
		return 0;
}


int NeqCheckInputs(NeQuickInputData_st* pstNeQuickInputData)
{
	int check1 = fabs(pstNeQuickInputData->pdGssPosLLH[0]) <= 90;
	int check2 = fabs(pstNeQuickInputData->pdSatPosLLH[0]) <= 90;
	int check3 = (pstNeQuickInputData->siMonth >= 1) &&
		(pstNeQuickInputData->siMonth <= 12);
	int check4 = (pstNeQuickInputData->dUT >= 0) &&
		(pstNeQuickInputData->dUT < 24);
	int check5 = pstNeQuickInputData->siNumCoeff >= 1;
	int check6 = &(pstNeQuickInputData->pdCoeff) != 0;

	return (check1 && check2 && check3 && check4 && check5 && check6);
}
/*checked*/
double DoTECIntegration(IntegrateData_st* pstIntegrateData, int bVert, SPoint_st stP0,
	double* pdNmax, LayersProperties_st* pstLayers)
{

	double H1 = 1000;
	double H2 = 2000;
	double R = pstIntegrateData->pstGeom->stRay.dR;
	double S1a = sqrt(SQR(RADIUS_OF_EARTH + 1000) - SQR(R));
	double S1b = sqrt(SQR(RADIUS_OF_EARTH + 2000) - SQR(R));

	double AbsTolLo = pstIntegrateData->pstNeQuickInputData->pdKronrodTol[0];
	double AbsTolHi = pstIntegrateData->pstNeQuickInputData->pdKronrodTol[1];

	SPoint_st stP1 = pstIntegrateData->pstGeom->stP1;
	SPoint_st stP2 = pstIntegrateData->pstGeom->stP2;

	double HeightP0, HeightP1, HeightP2;
	double dTEC;

	if (bVert) {
		HeightP0 = stP0.dH;
		HeightP1 = stP1.dH;
		HeightP2 = stP2.dH;
	}
	else {
		HeightP0 = stP0.dS;
		HeightP1 = stP1.dS;
		HeightP2 = stP2.dS;
	}

	if (stP2.dH <= H1) {
		//both below 1000km
		dTEC = NeqIntegrate(pstIntegrateData, HeightP0, HeightP2, 1,
			&(pstIntegrateData->pstGeom->stPactual), pdNmax,
			pstLayers, AbsTolLo);
	}
	else {
		if (stP2.dH <= H2) {
			if (stP1.dH >= H1) {
				//satelete and receiver are both between 1000 and 2000
				dTEC = NeqIntegrate(pstIntegrateData, HeightP1, HeightP2, 1,
					&(pstIntegrateData->pstGeom->stPactual), pdNmax,
					pstLayers, AbsTolHi);
			}
			else {
				//ray path crosses 1000
				dTEC = NeqIntegrate(pstIntegrateData, HeightP0, S1a, 1,
					&(pstIntegrateData->pstGeom->stPactual), pdNmax,
					pstLayers, AbsTolLo)
					+ NeqIntegrate(pstIntegrateData, S1a, HeightP2, 1,
						&(pstIntegrateData->pstGeom->stPactual), pdNmax,
						pstLayers, AbsTolHi);

			}
		}
		else {
			if (stP1.dH >= H2) {
				//receiver and satelete are both above 2000km
				dTEC = NeqIntegrate(pstIntegrateData, HeightP1, HeightP2, 1,
					&(pstIntegrateData->pstGeom->stPactual), pdNmax,
					pstLayers, AbsTolHi);
			}
			else {
				if (stP1.dH >= H1) {
					//ray path crosses 2000
					dTEC = NeqIntegrate(pstIntegrateData, HeightP1, S1b, 1,
						&(pstIntegrateData->pstGeom->stPactual), pdNmax,
						pstLayers, AbsTolHi)
						+ NeqIntegrate(pstIntegrateData, S1b, HeightP2, 1,
							&(pstIntegrateData->pstGeom->stPactual), pdNmax,
							pstLayers, AbsTolHi);
				}
				else {
					//ray path crosses 1000 and 2000 km

					dTEC = NeqIntegrate(pstIntegrateData, HeightP0, S1a, 1,
						&(pstIntegrateData->pstGeom->stPactual), pdNmax,
						pstLayers, AbsTolLo)
						+ NeqIntegrate(pstIntegrateData, S1a, S1b, 1,
						&(pstIntegrateData->pstGeom->stPactual), pdNmax,
						pstLayers, AbsTolHi)
						+ NeqIntegrate(pstIntegrateData, S1b, HeightP2, 1,
						&(pstIntegrateData->pstGeom->stPactual), pdNmax,
						pstLayers, AbsTolHi);
					
				}
			}
		}
	}
	return dTEC;
}
/*
 ///////////////////////////////////////////////////////////////////////////////
 // NEQCALCMODIPAz.c
 ///////////////////////////////////////////////////////////////////////////////
 */
 /* checked */
double NeqCalcModip(double dLat, double dLng, MODIP_st* pstModip)
{
	int    lngp = 36;
	double dlatp = 5;
	double dlngp = 10;

	double lat1, lng1, dj, di;
	int sj, si;

	int k, j;

	double z1[4], z[4];

	if (dLat <= -90)
		return -90;
	else if (dLat >= 90)
		return 90;
	else {

		lng1 = (dLng + 180.0) / dlngp;
		sj = (int)floor(lng1) - 2;
		dj = lng1 - floor(lng1);

		if (sj < 0)
			sj += lngp;

		if (sj > (lngp - 3))
			sj -= lngp;

		lat1 = (dLat + 90.0) / dlatp + 1;
		si = (int)floor(lat1 - 1e-6) - 2;
		di = lat1 - si - 2;

		for (k = 1; k <= 4; k++) {
			for (j = 1; j <= 4; j++) {
				z1[j - 1] = (*pstModip->pdModip)[si + j][sj + k + 1];
			}
			z[k - 1] = NeqInterpolate(z1, di);
		}
		return NeqInterpolate(z, dj);
	}
}

/* checked */
double NeqInterpolate(double* pdZ, double dDeltaX)
{
	double g0, g1, g2, g3;
	double a0, a1, a2, a3;
	double dx;


	if (dDeltaX <= 1e-15)
		return pdZ[1];
	else {

		g0 = pdZ[2] + pdZ[1];
		g1 = pdZ[2] - pdZ[1];
		g2 = pdZ[3] + pdZ[0];
		g3 = (pdZ[3] - pdZ[0]) / 3;

		a0 = 9 * g0 - g2;
		a1 = 9 * g1 - g3;
		a2 = g2 - g0;
		a3 = g3 - g1;

		dx = 2 * dDeltaX - 1;
		/*added pow(d, i)*/
		/*interpolation changed*/
		return (a0 + a1 * pow(dx, 1) + a2 * pow(dx, 2) + a3 * pow(dx, 3)) / 16.0;
	}
}
/* checked */
double NeqModipToAz(double dModip, int siNumCoeff, double* pdCoeff)
{
	int i;
	double Az = 0;
	for (i = 0; i < siNumCoeff; i++)
		Az += pdCoeff[i] * pow(dModip, i);                        /*// int power*/

	if (Az < 0)
		Az = 0;

	return Az;
}
/*
 ///////////////////////////////////////////////////////////////////////////////
 // NEQGETRAYPROPERTIES.c
 ///////////////////////////////////////////////////////////////////////////////
 */
 /* checked */
int NeqGetRayProperties(SPoint_st* pstP1, SPoint_st* pstP2, SPoint_st* pstRay,
	double* pdZeta, double* pdSinSig, double* pdCosSig)
{

	if ((fabs(pstP2->dLat - pstP1->dLat) < 1e-5) && /* || instead of && -bug */
		(fabs(pstP2->dLng - pstP1->dLng) < 1e-5))
	{
		pstP2->dLng = pstP1->dLng;
	}

	pstP1->dR = pstP1->dH + RADIUS_OF_EARTH;
	pstP2->dR = pstP2->dH + RADIUS_OF_EARTH;

	NeqCalcRayProperties1(pstP1, pstP2, pstRay, pdZeta);

	if ((fabs(*pdZeta) > 90) && (pstRay->dR < RADIUS_OF_EARTH))
		return 1;

	if (pstRay->dR >= 0.1)
		NeqCalcRayProperties2(pstP1, pstP2, pstRay, pdSinSig, pdCosSig);

	return 0;
}

/* checked */
void NeqCalcRayProperties1(SPoint_st* pstP1, SPoint_st* pstP2, SPoint_st* pstRay,
	double* pdZeta)
{
	double CosDl12, SinDl12;
	double CosDel, SinDel;
	double SinSigp, CosSigp;
	double Delp, SinDelp, CosDelp;
	double SinPhp, CosPhp;
	double SinLamp, CosLamp;

	if ((fabs(pstP2->dLat - pstP1->dLat) < 1e-5) && /* || instead && -bug */
		(fabs(pstP2->dLng - pstP1->dLng) < 1e-5))
	{
		pstRay->dLat = pstP1->dLat;
		pstRay->dLng = pstP1->dLng;
		pstRay->dR = 0;
		*pdZeta = 0;
	}
	else {
		/* degtorad added for all of four statements -bug*/
		pstP1->dSinLat = sin(pstP1->dLat * DEGTORAD);
		pstP1->dCosLat = cos(pstP1->dLat * DEGTORAD);
		pstP2->dSinLat = sin(pstP2->dLat * DEGTORAD);
		pstP2->dCosLat = cos(pstP2->dLat * DEGTORAD);

		CosDl12 = cos((pstP2->dLng - pstP1->dLng) * DEGTORAD);
		SinDl12 = sin((pstP2->dLng - pstP1->dLng) * DEGTORAD);

		CosDel = pstP1->dSinLat * pstP2->dSinLat + pstP1->dCosLat * pstP2->dCosLat * CosDl12;
		SinDel = sqrt(1 - SQR(CosDel));

		*pdZeta = atan2(SinDel, CosDel - (pstP1->dR / pstP2->dR));

		SinSigp = (SinDl12 * pstP2->dCosLat / SinDel); /* no divisor -bug */ /* pstP2 instead of P1 -bug*/
		CosSigp = (pstP2->dSinLat - CosDel * pstP1->dSinLat) / (SinDel * pstP1->dCosLat); /* / no divisor -bug */
		Delp = 0.5 * PI - *pdZeta;
		SinDelp = sin(Delp);
		CosDelp = cos(Delp);
		SinPhp = pstP1->dSinLat * CosDelp - pstP1->dCosLat * SinDelp * CosSigp;
		CosPhp = sqrt(1 - SQR(SinPhp));

		pstRay->dLat = atan2(SinPhp, CosPhp) * RADTODEG;

		SinLamp = -(SinSigp * SinDelp) / CosPhp;
		CosLamp = (CosDelp - pstP1->dSinLat * SinPhp) / (pstP1->dCosLat * CosPhp);

		if (fabs(fabs(pstP1->dLat) - 90.0) < 1e-10) {
			if (*pdZeta >= 0.0) {
				pstRay->dLng = pstP2->dLng + 180.0;
			}
			else {
				pstRay->dLng = pstP2->dLng;
			}
		}
		else {
			pstRay->dLng = atan2(SinLamp, CosLamp) * RADTODEG + pstP1->dLng; //2.5.8.2.3
		}
		//pstRay->dLng = atan2(SinLamp, CosLamp) * RADTODEG;
		pstRay->dR = pstP1->dR * sin(*pdZeta);

		*pdZeta *= RADTODEG;
	}
}

/* checked */
void NeqCalcRayProperties2(SPoint_st* pstP1, SPoint_st* pstP2, SPoint_st* pstRay,
	double* pdSinSig, double* pdCosSig)
{
	double DeltaLong;
	double CosPsi, SinPsi;

	pstP1->dSinLat = sin(pstRay->dLat * DEGTORAD);
	pstP1->dCosLat = cos(pstRay->dLat * DEGTORAD);
	pstP2->dSinLat = sin(pstP2->dLat * DEGTORAD); /* pstRay instead of pstP2 -bug */
	pstP2->dCosLat = cos(pstP2->dLat * DEGTORAD); /* pstRay instead of pstP2 -bug */ /* ->dSinlat instead of dCosLat -bug*/

	DeltaLong = (pstP2->dLng - pstRay->dLng) * DEGTORAD; /* ?? P1 - ray???*/

	if (fabs(fabs(pstRay->dLat) - 90) < 1e-10) {

		*pdSinSig = 0;

		if (pstRay->dLat > 0)
			*pdCosSig = -1;
		else
			*pdCosSig = 1;
	}
	else {

		CosPsi = pstP1->dSinLat * pstP2->dSinLat + pstP1->dCosLat * pstP2->dCosLat * cos(DeltaLong);

		SinPsi = sqrt(1 - SQR(CosPsi));

		*pdSinSig = (pstP2->dCosLat * sin(DeltaLong)) / SinPsi;
		*pdCosSig = (pstP2->dSinLat - pstP1->dSinLat * CosPsi) /
			(SinPsi * pstP1->dCosLat);
	}
}
/*
 ///////////////////////////////////////////////////////////////////////////////
 // NEQINTEGRATE.c
 ///////////////////////////////////////////////////////////////////////////////
 */
 /*checked*/
double NeqIntegrate(IntegrateData_st* pstIntegrateData, double dH1, double dH2,
	int siCurrentLevel, SPoint_st* pstPactual, double* pdNmax,
	LayersProperties_st* pstLayers, double AbsTol)
{
	const double wi[15] =
	{ 0.022935322010529224963732008058970,
		0.063092092629978553290700663189204,
		0.104790010322250183839876322541518,
		0.140653259715525918745189590510238,
		0.169004726639267902826583426598550,
		0.190350578064785409913256402421014,
		0.204432940075298892414161999234649,
		0.209482141084727828012999174891714,
		0.204432940075298892414161999234649,
		0.190350578064785409913256402421014,
		0.169004726639267902826583426598550,
		0.140653259715525918745189590510238,
		0.104790010322250183839876322541518,
		0.063092092629978553290700663189204,
		0.022935322010529224963732008058970 };
	const double wig[7] =
	{ 0.129484966168869693270611432679082,
		0.279705391489276667901467771423780,
		0.381830050505118944950369775488975,
		0.417959183673469387755102040816327,
		0.381830050505118944950369775488975,
		0.279705391489276667901467771423780,
		0.129484966168869693270611432679082 };
	const double xi[15] =
	{ -0.991455371120812639206854697526329,
		-0.949107912342758524526189684047851,
		-0.864864423359769072789712788640926,
		-0.741531185599394439863864773280788,
		-0.586087235467691130294144838258730,
		-0.405845151377397166906606412076961,
		-0.207784955007898467600689403773245,
		0,
		0.207784955007898467600689403773245,
		0.405845151377397166906606412076961,
		0.586087235467691130294144838258730,
		0.741531185599394439863864773280788,
		0.864864423359769072789712788640926,
		0.949107912342758524526189684047851,
		0.991455371120812639206854697526329 };

	double h2, hh;
	int i, Gind;
	double x, y;
	double intk = 0;
	double intg = 0;

	h2 = 0.5 * (dH2 - dH1);
	hh = 0.5 * (dH2 + dH1); /* * instead of + -bug*/

	for (i = 0, Gind = 0; i < 15; i++) {

		x = h2 * xi[i] + hh;

		if (pstIntegrateData->bVert)
			y = NeqGetNeOnVertRay(x, pstLayers, pdNmax);
		else
			y = NeqGetNeOnSlantRay(x, pstIntegrateData, pdNmax,
				pstPactual, pstLayers);
		intk += y * wi[i];

		if ((i % 2) == 1) {

			intg += y * wig[Gind];
			Gind++;
		}
	}

	intk *= h2;
	intg *= h2;

	if (fabs(intk - intg) <= AbsTol || fabs((intk - intg) / intk) <= AbsTol)         /* подумать.. fabs ( ( - )/ ) */ /*// or RelTol ??*/
		return intk;
	else if (siCurrentLevel == pstIntegrateData->pstNeQuickInputData->siMaxRecurse)
		return intk;
	else
		return
		(NeqIntegrate(pstIntegrateData, dH1, (dH1 + h2), siCurrentLevel + 1,
			pstPactual, pdNmax, pstLayers, AbsTol)
			+ NeqIntegrate(pstIntegrateData, (dH1 + h2), dH2, siCurrentLevel + 1,
				pstPactual, pdNmax, pstLayers, AbsTol));
}
/*
 ///////////////////////////////////////////////////////////////////////////////
 // NEQGETNEONVERTRAY.c
 ///////////////////////////////////////////////////////////////////////////////
 */
 /* checked */
double NeqGetNeOnVertRay(double dH, LayersProperties_st* pstLayers, double* pdNmax) {

	if (dH > pstLayers->pdPeakHeight[F2])
		return NeqCalcTopsideNe(pdNmax, pstLayers, dH);
	else
		return NeqCalcBottomsideNe(dH, pstLayers);

}
/*
 ///////////////////////////////////////////////////////////////////////////////
 // NEQGETNEONSLANTRAY.c
 ///////////////////////////////////////////////////////////////////////////////
 */
 /* checked */
double NeqGetNeOnSlantRay(double dS, IntegrateData_st* pstIntegrateData, double* pdNmax,
	SPoint_st* pstPactual, LayersProperties_st* pstLayers)
{
	Geometry_st* pstGeom = pstIntegrateData->pstGeom;

	NeqCalcLLHOnRay(dS, &(pstGeom->stRay), &(pstGeom->stP1),
		pstGeom->dSinSig, pstGeom->dCosSig, pstPactual);
	NeqCalcEpstParams(pstIntegrateData->pstNeQuickInputData, pstPactual,
		pstGeom->dSinDelta, pstGeom->dCosDelta, pdNmax,
		pstLayers, pstIntegrateData->pstCurrCCIR);

	if (pstPactual->dH > pstLayers->pdPeakHeight[F2])
		return NeqCalcTopsideNe(pdNmax, pstLayers, pstPactual->dH);
	else
		return NeqCalcBottomsideNe(pstPactual->dH, pstLayers);
}

/* checked */
void NeqCalcLLHOnRay(double dS, SPoint_st* pstRay, SPoint_st* pstP1,
	double dSinSig, double dCosSig, SPoint_st* pstPactual)
{
	double TanDel, CosDel, SinDel;
	double arg, CLong;

	pstPactual->dCosLat = cos(pstPactual->dLat * DEGTORAD);
	pstPactual->dSinLat = sin(pstPactual->dLat * DEGTORAD);

	TanDel = dS / (pstRay->dR);
	CosDel = 1 / sqrt(1 + SQR(TanDel));
	SinDel = TanDel * CosDel;

	arg = pstP1->dSinLat * CosDel + pstP1->dCosLat * SinDel * dCosSig;

	pstPactual->dLat = atan2(arg, sqrt(1 - SQR(arg))) * RADTODEG;
	CLong = atan2(SinDel * dSinSig * pstP1->dCosLat, CosDel - pstP1->dSinLat *
		arg) * RADTODEG;
	pstPactual->dLng = CLong + pstRay->dLng;
	pstPactual->dH = sqrt(SQR(dS) + SQR(pstRay->dR)) - RADIUS_OF_EARTH; /* sqr) - Re, no)g -bug*/
	pstPactual->dCosLat = cos(pstPactual->dLat * DEGTORAD); /*added */
	pstPactual->dSinLat = sin(pstPactual->dLat * DEGTORAD);
	/*pstPactual->dS = sqrt( SQR(pstPactual->dR) - SQR(pstRay->dR) )*/
}
/*
 ///////////////////////////////////////////////////////////////////////////////
 // NEQCALCEPSTPARAMS.c
 ///////////////////////////////////////////////////////////////////////////////
 */
 /* checked */
void NeqCalcEpstParams(NeQuickInputData_st* pstNeQuickInputData, SPoint_st* pstPactual,
	double dSinDelta, double dCosDelta, double* pdNmax,
	LayersProperties_st* pstLayers, CurrentCCIR_st* pstCurrCCIR)
{
	double modip, Az, R12, RR1, RR2;

	CCIR_st* pstCCIR = pstNeQuickInputData->pstCCIR;
	int      mth = pstNeQuickInputData->siMonth;
	int      mth1 = mth - 1;
	double   ut = pstNeQuickInputData->dUT;

	double   SinNModip[12];
	double   xlt;
	double   CosChi, chi, chi0;
	int      seas = 1;
	double   ee, chin, sfac, fa;
	double   b2k, x, v, NdHmx;

	double   Nm[3];
	double*   PeakHeight = pstLayers->pdPeakHeight;
	double*   BotThick = pstLayers->pdBotThick;
	double*   TopThick = pstLayers->pdTopThick;
	double*   Amp = pstLayers->pdAmp;
	double*   F0 = pstLayers->pdF0;
	double*  M3000 = &pstLayers->dM3000;

	int i, j;
	modip = NeqCalcModip(pstPactual->dLat, pstPactual->dLng, pstNeQuickInputData->pstModip);
	Az = pstNeQuickInputData->dAzBase;
	R12 = sqrt(1123.6 * (Az - 63.7) + 167273) - 408.99;
	*pdNmax = -1;

	/*// todo: change stiff condition for R12 to nonstiff ??*/
	if ((mth != pstCurrCCIR->siMonth) || (R12 != pstCurrCCIR->dR12)) {

		RR2 = 0.01 * R12;
		RR1 = 1 - RR2;

		/*// ! index allignment changed*/
		for (i = 0; i <= 12; i++) {
			for (j = 0; j <= 75; j++) {
				pstCurrCCIR->pdF0F2[j * 13 + i] = (*pstCCIR->pdF2)[mth1][0][j][i] * RR1
					+ (*pstCCIR->pdF2)[mth1][1][j][i] * RR2;
			}
		}

		for (i = 0; i <= 8; i++) {
			for (j = 0; j <= 48; j++) {
				pstCurrCCIR->pdM3000F2[j * 9 + i] = (*pstCCIR->pdM3000)[mth1][0][j][i] * RR1
					+ (*pstCCIR->pdM3000)[mth1][1][j][i] * RR2;
			}
		}
		pstCurrCCIR->siMonth = mth;
		pstCurrCCIR->dR12 = R12;

		NeqCalcSphLegCoeffs(ut, pstCurrCCIR);
	}
	else if (ut != pstCurrCCIR->dUT)
		NeqCalcSphLegCoeffs(ut, pstCurrCCIR);

	SinNModip[0] = 1;
	SinNModip[1] = sin(modip * DEGTORAD); /* DEGTORAD added -bug*/
	for (i = 2; i <= 11; i++)
		SinNModip[i] = SinNModip[1] * SinNModip[i - 1];
	/*pstPactual->dCosLat = cos(pstPactual->dLat * DEGTORAD ); added coslat (to LLH function) -bug*/
	F0[F2] = NeqGetF2FreqFromCCIR(pstPactual->dCosLat, pstPactual->dLng,
		pstCurrCCIR->pdLegCoeffs_F0, SinNModip, 0);

	*M3000 = NeqGetF2FreqFromCCIR(pstPactual->dCosLat, pstPactual->dLng,       /* m3000 < 0 !?*/
		pstCurrCCIR->pdLegCoeffs_M3000, SinNModip, 1);

	xlt = ut + pstPactual->dLng / 15.0;
	if (xlt < 0)
		xlt = xlt + 24.0;
	else if (xlt >= 24)
		xlt = xlt - 24.0;
	/* sinlat and coslat are already calculated and contained in pstactual data */
	CosChi = pstPactual->dSinLat * dSinDelta + pstPactual->dCosLat * dCosDelta * cos(PI * (12 - xlt) / 12);
	chi = atan2(sqrt(1 - SQR(CosChi)), CosChi) * RADTODEG;
	chi0 = 86.23292796211615;

	switch (mth) {

	case 1: seas = -1;
		break;
	case 2: seas = -1;
		break;
	case 3: seas = 0;
		break;
	case 4: seas = 0;
		break;
	case 5: seas = 1;
		break;
	case 6: seas = 1;
		break;
	case 7: seas = 1;
		break;
	case 8: seas = 1;
		break;
	case 9: seas = 0;
		break;
	case 10: seas = 0;
		break;
	case 11: seas = -1;
		break;
	case 12: seas = -1;
	}


	chin = NeqJoin(90.0 - 0.24 * NeqClipExp(20.0 - 0.2 * chi), chi, 12, chi - chi0);
	ee = NeqClipExp(0.3 * pstPactual->dLat);                                            /* ! */
	double const seasp = seas * (ee - 1) / (ee + 1);
	sfac = (1.112 - 0.019 * seasp) * sqrt(sqrt(Az));
	fa = sfac * NeqClipExp(log(cos(chin * DEGTORAD)) * 0.3);


	F0[E] = sqrt(SQR(fa) + 0.49);//foE

	F0[F1] = NeqJoin(1.4 * F0[E], 0, 1000.0, F0[E] - 2);
	F0[F1] = NeqJoin(0, F0[F1], 1000.0, F0[E] - F0[F1]);
	F0[F1] = NeqJoin(F0[F1], 0.85 * F0[F1], 60.0, 0.85 * F0[F2] - F0[F1]);//f0F1

	if (F0[F1] < 1e-6)
		F0[F1] = 0;

	Nm[F2] = NeqCriticalFreqToNe(F0[F2]);
	Nm[F1] = NeqCriticalFreqToNe(F0[F1]);
	Nm[E] = NeqCriticalFreqToNe(F0[E]);

	PeakHeight[F2] = NeqCalcF2PeakHeight(*M3000, F0[E], F0[F2]);
	PeakHeight[E] = 120;
	PeakHeight[F1] = 0.5 * (PeakHeight[F2] + PeakHeight[E]);

	NdHmx = 0.01 * exp(-3.467 + 0.857 * log(SQR(F0[F2])) + 2.02 * log(*M3000));

	BotThick[F2] = 0.385 * Nm[F2] / NdHmx;
	TopThick[F1] = 0.3 * (PeakHeight[F2] - PeakHeight[F1]);
	BotThick[F1] = 0.5 * (PeakHeight[F1] - PeakHeight[E]); /* .5 = 0.5? */
	TopThick[E] = BotThick[F1];

	if (TopThick[E] < 7)
		TopThick[E] = 7;

	BotThick[E] = 5;
	/* code (Amp[i] = ...) missed -bug*/
	Amp[F2] = 4 * Nm[F2];
	Amp[F1] = 4 * Nm[F1];
	Amp[E] = 4 * Nm[E];

	if (F0[F1] < 0.5) {

		Amp[F1] = 0;
		Amp[E] = 4 * (Nm[E] - NeqCalcEpstein(Amp[F2], PeakHeight[F2],
			BotThick[F2], PeakHeight[E]));
	}
	else {

		for (i = 0; i < 5; i++) {

			Amp[F1] = 4 * (Nm[F1]
				- NeqCalcEpstein(Amp[F2], PeakHeight[F2], BotThick[F2], PeakHeight[F1])
				- NeqCalcEpstein(Amp[E], PeakHeight[E], BotThick[E], PeakHeight[F1]));
			Amp[F1] = NeqJoin(Amp[F1], 0.8 * Nm[F1], 1, Amp[F1] - 0.8 * Nm[F1]);
			Amp[E] = 4 * (Nm[E]
				- NeqCalcEpstein(Amp[F1], PeakHeight[F1], BotThick[F1], PeakHeight[E])
				- NeqCalcEpstein(Amp[F2], PeakHeight[F2], BotThick[F2], PeakHeight[E]));
		}
	}

	Amp[E] = NeqJoin(Amp[E], 0.05, 60.0, Amp[E] - 0.005);
	if ((mth > 3) && (mth < 10))
		b2k = 6.705 - 0.014 * R12 - 0.008 * PeakHeight[F2];
	else
		b2k = -7.77 + 0.097 * SQR(PeakHeight[F2] / BotThick[F2]) + 0.153 * Nm[F2];

	b2k = NeqJoin(b2k, 2, 1, b2k - 2.0);
	b2k = NeqJoin(8, b2k, 1, b2k - 8.0);

	TopThick[F2] = b2k * BotThick[F2];
	x = 0.01 * (TopThick[F2] - 150.0);
	v = (0.041163 * x - 0.183981) * x + 1.424472;
	TopThick[F2] /= v;
}
/* checked */
void NeqCalcSphLegCoeffs(double dUT, CurrentCCIR_st* pstCurrCCIR)//__F2fouriertimeseries__
{
	double* LegCoeffs_F0 = pstCurrCCIR->pdLegCoeffs_F0;
	double* LegCoeffs_M3000 = pstCurrCCIR->pdLegCoeffs_M3000;
	double* CCIR_F0 = pstCurrCCIR->pdF0F2;
	double* CCIR_M3000 = pstCurrCCIR->pdM3000F2;

	double t = (dUT * 15 - 180) * DEGTORAD;

	double sinHarm[6];
	double cosHarm[6];

	int i, i1, k;

	sinHarm[0] = sin(t);
	cosHarm[0] = cos(t);

	for (i = 1; i <= 5; i++) {

		i1 = i - 1;
		sinHarm[i] = sinHarm[i1] * cosHarm[0] + cosHarm[i1] * sinHarm[0];
		cosHarm[i] = cosHarm[i1] * cosHarm[0] - sinHarm[i1] * sinHarm[0];
	}

	/*// ! indexes changed*/
	for (i = 0; i <= 75; i++) {

		LegCoeffs_F0[i] = CCIR_F0[i * 13];
		for (k = 0; k <= 5; k++)
			LegCoeffs_F0[i] += CCIR_F0[i * 13 + 2 * (k + 1) - 1] * sinHarm[k]
			+ CCIR_F0[i * 13 + 2 * (k + 1)] * cosHarm[k];
	}

	for (i = 0; i <= 48; i++) {

		LegCoeffs_M3000[i] = CCIR_M3000[i * 9];
		for (k = 0; k <= 3; k++)
			LegCoeffs_M3000[i] += CCIR_M3000[i * 9 + 2 * (k + 1) - 1] * sinHarm[k] /*is was 0; now it's - 1 -bug*/
			+ CCIR_M3000[i * 9 + 2 * (k + 1)] * cosHarm[k];                        /*it was +1; now it's 0*/
	}

}
/*checked*/
double NeqGetF2FreqFromCCIR(double dCosLat, double dLng, double* pdLegCoeffs,
	double* pdSinModipToN, int siMode)
{
	int R;
	double term1 = 0;
	double term2 = 0;

	int i, l, k;

	int k1;
	int* nq;

	int nq0[9] = { 11, 11, 8, 4, 1, 0, 0, 0, 0 };
	int nq1[7] = { 6, 7, 5, 2, 1, 0, 0 };

	if (siMode == F2) {
		nq = nq0;
		k1 = 9;
	}
	else {
		nq = nq1;
		k1 = 7;
	}

	for (i = 0; i <= 11; i++)
		if (fabs(pdSinModipToN[i]) <= 1e-30)
			pdSinModipToN[i] = 0;

	for (i = 0; i <= nq[0]; i++)
	{
		term1 += pdLegCoeffs[i] * pdSinModipToN[i];
	}

	for (k = 1; k <= k1 - 1; k++)
		for (i = 0; i <= nq[k]; i++) {

			R = 0;
			for (l = 0; l <= k - 1; l++)
				R += nq[l];
			R += k;
			R *= 2;
			R -= nq[0] + 1;

			/*// int pow*/
			term2 += pdSinModipToN[i] * pow(dCosLat, k) *
				(pdLegCoeffs[R + 2 * i] * cos(k * dLng * DEGTORAD) + /* degtorad added -bug */
					pdLegCoeffs[R + 2 * i + 1] * sin(k * dLng * DEGTORAD));
		}
	assert(term2 + term1);
	return term2 + term1;
}
/* checked */
double NeqCriticalFreqToNe(double dF0)
{
	return 0.124 * SQR(dF0);
}
/* checked */
double NeqCalcF2PeakHeight(double dM3000, double dF0E, double dF0F2)
{
	/*if (dM3000 < 0)
		dM3000 *= -1;
	if (abs(dM3000 < 1))
		dM3000 = 1;*/
	double MF = dM3000 * sqrt((0.0196 * SQR(dM3000) + 1) /
		(1.2967 * SQR(dM3000) - 1));

	double r, exp_term, r2, deltaM;

	if (dF0E >= 1e-30)
	{
		r = dF0F2 / dF0E;
		exp_term = exp(20 * (r - 1.75));
		r2 = (r * exp_term + 1.75) / (exp_term + 1);
		assert(r2 > 1.73);
		deltaM = 0.253 / (r2 - 1.215) - 0.012;
	}
	else
		deltaM = -0.012;
	assert(dM3000 >= 1.0);
	assert(deltaM >= -0.012);
	return 1490 * MF / (dM3000 + deltaM) - 176;
}
/* checked */
double NeqCalcEpstein(double dNmax, double dHmax, double dB, double dH)
{
	double ExpTerm = NeqClipExp((dH - dHmax) / dB);
	return dNmax * ExpTerm / SQR(1 + ExpTerm);
}
/*
 ///////////////////////////////////////////////////////////////////////////////
 // NEQCALCTOPSIDE.c
 ///////////////////////////////////////////////////////////////////////////////
 */
double NeqCalcTopsideNe(double* pdNmax, LayersProperties_st* pstLayers, double dH)
{
	double g = 0.125;
	double rfac = 100;

	double* PeakHeight = pstLayers->pdPeakHeight;
	double* TopThick = pstLayers->pdTopThick;

	double ee = NeqClipExp((dH - PeakHeight[F2]) / (TopThick[F2] *
		(1 + rfac * g * (dH - PeakHeight[F2]) /
		(rfac * TopThick[F2] + g * (dH - PeakHeight[F2])))));
	double ep;

	if (ee > 1e11)
		ep = 4 / ee;
	else
		ep = 4 * ee / SQR(1 + ee);

	/*// if pdNmax has not been calculated*/
	if (*pdNmax < 0)
		*pdNmax = NeqCalcBottomsideNe(PeakHeight[F2], pstLayers);

	return (*pdNmax) * ep;
}
/*
 ///////////////////////////////////////////////////////////////////////////////
 // NEQCALCBOTTOMSIDE.c
 ///////////////////////////////////////////////////////////////////////////////
 */
 /* checked */
double NeqCalcBottomsideNe(double dHH, LayersProperties_st* pstLayers)
{
	double h0 = 100;
	double Hd = 10;
	double f1 = 10;
	double f2 = 1;

	double* PeakHeight = pstLayers->pdPeakHeight;
	double* TopThick = pstLayers->pdTopThick;
	double* BotThick = pstLayers->pdBotThick;
	double* Amp = pstLayers->pdAmp;

	double  B[3], s[3], ds[3];
	double  arg, z, bf, aN0;
	double  sum_sj, sum_dsj_sj;
	int j;

	B[F2] = BotThick[F2];

	if (dHH > PeakHeight[F1])
		B[F1] = TopThick[F1];
	else
		B[F1] = BotThick[F1];

	if (dHH > PeakHeight[E])
		B[E] = TopThick[E];
	else
		B[E] = BotThick[E];

	if (dHH < 100) {
		for (j = 0; j <= 2; j++) {

			if (j == F2)
				arg = (h0 - PeakHeight[F2]) / B[F2];
			else
				arg = ((h0 - PeakHeight[j]) / B[j]) *
				exp(f1 / (1 + f2 * fabs(h0 - PeakHeight[F2])));

			if (fabs(arg) > 25) {
				s[j] = 0;
				ds[j] = 0;
			}
			else {
				s[j] = Amp[j] * exp(arg) / SQR(1 + exp(arg));
				ds[j] = (1 - exp(arg)) / (B[j] * (1 + exp(arg)));
			}
		}

		z = (dHH - h0) / Hd;
		sum_sj = (s[0] + s[1] + s[2]);
		sum_dsj_sj = (s[0] * ds[0] + s[1] * ds[1] + s[2] * ds[2]);

		bf = 1 - Hd * sum_dsj_sj / sum_sj; /* NOT sum_sj_sqr BUT sum(sj* dsj) -bug*/
		aN0 = 1e11 * sum_sj;

		double const ans = aN0 * NeqClipExp(1 - bf * z - NeqClipExp(-z));
		return ans;
	}
	else {
		for (j = 0; j <= 2; j++) {

			if (j == F2)
				arg = (dHH - PeakHeight[F2]) / B[j];
			else
				arg = ((dHH - PeakHeight[j]) / B[j]) *
				exp(f1 / (1 + f2 * fabs(dHH - PeakHeight[F2])));

			if (fabs(arg) > 25)
				s[j] = 0;
			else
				s[j] = Amp[j] * exp(arg) / SQR(1 + exp(arg));
		}

		sum_sj = (s[0] + s[1] + s[2]);
		return 1e11 * sum_sj;
	}
}
/*
 ///////////////////////////////////////////////////////////////////////////////
 // NEQUTILS.c
 ///////////////////////////////////////////////////////////////////////////////
 */
double NeqJoin(double dF1, double dF2, double dAlpha, double dX)
{
	double ee = NeqClipExp(dAlpha * dX);
	return (dF1 * ee + dF2) / (ee + 1);
}

double NeqClipExp(double dPower)
{
	if (dPower > 80)
		return 5.5406e34;

	if (dPower < -80)
		return 1.8049e-35;

	return exp(dPower);
}
