#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <random>
#include </home/jallen/root-6.08/include/TROOT.h>
#include </home/jallen/root-6.08/include/TH1.h>
#include </home/jallen/root-6.08/include/TH2.h>
#include </home/jallen/root-6.08/include/TMath.h>
#include </home/jallen/root-6.08/include/TFile.h>
#include </home/jallen/root-6.08/include/TTree.h>
#include </home/jallen/root-6.08/include/TF1.h>
#include </home/jallen/root-6.08/include/TVector3.h>
#include </home/jallen/root-6.08/include/TRandom3.h>
#include <fstream>

using namespace std;

static TRandom3 *rando = new TRandom3(0.);

class Common
{
public:
	int aISG[19], aINN[19];
	double aDEN[19], aTHK[19], aPRS[19], aXLN[19], aARDEN[19];
	double aZNUMB[19][4], aANUMB[19][4], aELNUM[19][4], aCONCN[19][4];
	double aPRDEN[19][4], aPRTHK[19][4];
	double loste[19];

	double desorb(int ianz, double zp, double ap, double ep);
	double ads(int * pI1, double * psign, double * pXN1, double * pEPS, double * pA, double * pZ, double * pE, double * pDEE, int * pISTAT);
	double setabs(int * pINW, double A[], double Z[], double AN[], double T[], double * pTH, double D[], double * pDN);
	double setabg(int * pINW, double A[], double Z[], double AN[], double CN[], double T[], double * pTH, double D[], double * pDN, double * pPR, double * pXL);
	double dedx(double * pZ1, double * pA1, double * pZ2, double * pA2, double * pRHO, double * pENER, double * pV, int * pIFG, double * pDEDXHI, double * pDEDXTO);

};
double Common::desorb(int ianz, double zp, double ap, double ep)
{
	double aZNUMBW[4], aANUMBW[4], aELNUMW[4], aCONCNW[4];
	double aPRDENW[4], aPRTHKW[4];
	double aE[19], aDE[19], aXmem[19];
	double aTOUT[19], aTOUTE[19];
	double aEptable[50][500][2], aEmintabz[50];
	int INNW, ISTAT, ISTORE, INS, iArrayValue, aIANZV[5];
	int io1, IO2, IO3, io0, iopt, ianzi, ianzide, ISN, I, J, izp, I1STPASS, ipunch;
	double DENW, THKW, PRSW, XLNW, EI, EPS, XUPDN, XNS, DEI, EIOLD;

	io1 = 9;
	IO2 = 11;
	IO3 = 12;
	io0 = 1;
	iopt = 1;
	ianzi = 2;
	ianzide = 1;
	ISTAT = 0;
	XNS = 2.0;
	DEI = 0.0;
	


	/* the mass table is to be used only for iopt = 5,6
	use atomic masses to average for isotipic composition.
	taken from Formulas Facts and Constants, H. J. Fischbeck 
	K. H. Fischbeck. Springer - Verlag 1987 2nd ed, pages 164-183.*/

	double amass[70] = { 1.01, 4.00, 6.94, 9.01, 10.81, 12.01, 14.01, 16.00, 19.00,
		20.18, 22.99, 24.31, 26.98, 28.09, 30.97, 32.07, 35.45, 39.95,
		39.10, 40.08, 44.96, 47.88, 50.94, 52.00, 54.94, 55.85, 58.93,
		58.69, 63.55, 65.39, 69.72, 72.59, 74.92, 78.96, 79.90, 83.80,
		85.47, 87.62, 88.91, 91.22, 92.91, 95.94, 98., 101.07, 102.91,
		106.42, 107.87, 112.41, 114.82, 118.71, 121.75, 127.60, 126.90,
		131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24, 147.,
		150.36, 151.96, 157.25, 158.93, 162.5, 164.93, 167.26, 168.93,
		173.04 };

	// For Z > 70, you are in trouble!

	

/*************************************************************
	  IOPT = 1 - SUPPLY ENERGY OF PARTICLE ENTERING
                 THE ABSORBER ARRAY AND GET LOSS AND
                 RANGES
      IOPT = 2 - SUPPLY TARGET, PROECTILE AND EJECTILE
                 INFO. AND THEN GO THROUGH ABSORBER
                 SANDWICH
      IOPT = 3 - CALCULATE ENERGY DEPOSITS IN DETECTOR
                 DETECTOR DIMENSIONS ARE STANDARD AND
                 THE VARIABLE -'IDET' - CHOOSES BETWEEN
                 VARIETY OF AVAILABLE DETECTORS
      IOPT = 4 - FINDS MAXIMUM ENERGY THAT CAN BE STOPPED IN
                 IANZ ELEMENTS OF THE SANDWICH FOR GIVEN
                 ZP, AP.
                 WHEN CALCULATION IS FINISHED, THE PROGRAM READS
                 IN NEW VALUES OF ZP, AP AND RESTARTS. TO END
                 THE PROGRAM, GIVE ZP < 0.
                 IN ORDER TO HELP THE SPEED OF THE PROGRAM,
                 GIVE THE PARTICLE'S "Z" IN increasing  ORDER.
      IOPT = 5 - STORES ARRAYS OF Edet AS A FUNCTION OF INCIDENT    
                 ENERGY AND THE PARTICLE'S ID (Z,A)
                 ARRAY LABELED  eptable(Z-Zth,Einc,ipunch)
                 ipunch = 1 stopped,  = 2 punched through
                 Einc = E(incident)/detable
                 Zth  = lowest Z considered - 1
*************************************************************/

	
/*    ianz  = number of elements in absorber "sandwich" - the 
              particle deposits all its energy in these layers.
      ianzi = index of last layer in which the energy of the 
              particle is not recorded - this unreecorded energy
              is used in the DST production in two modes:
              when making DST tape from data:
              Since the detector records only deposited energy
              the output tables are used to correct for this 
              deficinecy and for a given dharge and mass extrapolate
              the measured energy to the target exit energy
              when making DST tape from model calculations:
              The "lost" energy is a source of broadening of the 
              energy spectra due to straggling - this "smearing" 
              is estimated and superposed on the calculated spectra.
	   ianzide = element # for DE calculation*
	   
      zth= threshold Z for incident energy tble calc. (iopt=5,6)
      zmax = maximum z for table calculation
      detable = the energy step size used for array storage   
      emin = starting incident energy for table
      emax = mazximum incident energy for table calculation
      EP is ignored when iopt = 5 or 6*/

	for (ISN = 0; ISN < ianz; ISN++)
	{
		if (Common::aISG[ISN] != 1)
		{
			aTOUT[ISN] = Common::aTHK[ISN] / (Common::aDEN[ISN] * 1000.0);
			aTOUTE[ISN] = aTOUT[ISN] / 2.54;
		}

	}


	for (I = 0; I < ianz; I++)
	{
		INNW = Common::aINN[I];
		for (J = 0; J < INNW; J++)
		{
			aANUMBW[J] = Common::aANUMB[I][J];
			aZNUMBW[J] = Common::aZNUMB[I][J];
			aELNUMW[J] = Common::aELNUM[I][J];
			aCONCNW[J] = Common::aCONCN[I][J];


		}
		DENW = Common::aDEN[I];
		XLNW = Common::aXLN[I];
		PRSW = Common::aPRS[I];
		THKW = Common::aTHK[I];

		if (Common::aISG[I] == 1)
		{
			Common::setabg( &INNW, aANUMBW, aZNUMBW, aELNUMW, aPRTHKW, aCONCNW, &THKW, aPRDENW, &DENW, &PRSW, &XLNW);
			Common::aDEN[I] = 0.0;
			Common::aTHK[I] = 0.0;

			for (J = 0; J < INNW; J++)
			{
				Common::aPRDEN[I][J] = aPRDENW[J];
				Common::aPRTHK[I][J] = aPRTHKW[J];
				Common::aDEN[I] = Common::aDEN[I] + Common::aPRDEN[I][J];
				Common::aTHK[I] = Common::aTHK[I] + Common::aPRTHK[I][J];
			}
		}
		else if (Common::aISG[I] != 1)
		{
			Common::setabs(&INNW, aANUMBW, aZNUMBW, aELNUMW, aPRTHKW, &THKW, aPRDENW, &DENW);
			for (J = 0; J < INNW; J++)
			{
				Common::aPRDEN[I][J] = aPRDENW[J];
				Common::aPRTHK[I][J] = aPRTHKW[J];
			}

		}
	}
	zp = zp + 0.5;
	izp = (int)zp;
	zp = zp - 0.5;

	// EI energy in
	EI = ep;
	XUPDN = -1.0;
	EPS = 0.0001;
	I1STPASS = 1;

	ipunch = 2;


	for (I = 0; I < ianz; I++)
	{

		if (Common::aISG[I] == 0)
		{

		}
		else if (Common::aISG[I] == 1)
		{ 		

		}
		
		for (J = 0; J < Common::aINN[I]; J++)
		{

		}
	}
	/*XNS Initial nuber of intervals for integration of DE
	DEI energy out- energy in(<0 for energy loss)
	E[I] Energy left after I'th element (EP-DE[0] - DE[1] + .....)
	if particle stopped in detector this is equal to energy lost in
	remaining layers */

	I = 0;
	Common::ads(&I, &XUPDN, &XNS, &EPS, &ap, &zp, &EI, &DEI, &ISTAT);
	EIOLD = EI;
	aDE[I] = DEI;
	aE[I] = EI + DEI;
	EI = aE[I];
	XNS = XNS + 0.1;
	INS = (int)XNS;


	Common::loste[I] = -1.0 * aDE[I];

	if (EI < EPS || ISTAT == -1)
	{
	}

	iArrayValue = ianz + 1;
	Common::loste[iArrayValue] = aE[ianz];

	ISTORE = I;
	
	if (EI < EPS)
	{
		ipunch = 1;
	}


	return 0;

}

double Common::setabs(int * pINW, double A[], double Z[], double AN[], double T[], double * pTH, double D[], double * pDN)
{
	//Function for setting up composite absorber data
	//Partial Densities and Thicknesses

	int INW;
	double DN, TH;

	INW = *pINW;
	DN = *pDN;
	TH = *pTH;

	int I;
	double AW;

	AW = 0;

	for (I = 0; I < INW; I++)
	{
		AW = AW + (A[I] * AN[I]);
	}
	for (I = 0; I < INW; I++)
	{
		AN[I] = A[I] * AN[I] / AW;
		T[I] = TH * AN[I];
		D[I] = DN*AN[I];

		
	}

	*pINW = INW;
	*pDN = DN;
	*pTH = TH;

	return 0;
}

double Common::setabg(int * pINW, double A[], double Z[], double AN[], double CN[], double T[], double * pTH, double D[], double * pDN, double * pPR, double * pXL)
{
	double P, X, AWW, AW, TH, DN, PR, XL;
	int I, INW;
	TH = *pTH;
	DN = *pDN;
	PR = *pPR;
	XL = *pXL;
	INW = *pINW;

	P = PR / 760.0;
	X = XL / 22.4;
	AWW = 0.0;
	AW = 0.0;

	for (I = 0; I < INW; I++)
	{
		AW = AW + (A[I] * AN[I]);
		AWW = AWW + (A[I] * AN[I] * CN[I]);
		T[I] = P * X * A[I] * AN[I] * CN[I];
		D[I] = T[I] / XL;
	}

	*pTH = TH;
	*pDN = DN;
	*pPR = PR;
	*pXL = XL;
	*pINW = INW;

	return 0;

}

double Common::ads(int * pI1, double * psign, double * pXN1, double * pEPS, double * pA, double * pZ, double * pE, double * pDEE, int * pISTAT)
{
	//Subroutine for energy loss calculations
	//Call DEDX for stopping power calculations
	double EH, DEDNEXT, AX, ZX, FX, DENST, VH, DE, DEX, sign, XN1, EPS, A, Z, E, DEE;
	double DDD, DED1ST, DDS, DDR;
	int N1, K, J1, J, I, ISGW, I1, ISTAT;
	I1 = *pI1;
	ISTAT = *pISTAT;
	sign = *psign;
	XN1 = *pXN1;
	EPS = *pEPS;
	A = *pA;
	Z = *pZ;
	E = *pE;
	DEE = *pDEE;
	DE = 0.0;
	DEX = 0.0;
	// NI number of integrations for energy loss
	EH = E;
	XN1 = XN1 + 0.001;
	N1 = (int)XN1;
	DEDNEXT = 0.0;

	for (K = 1; K <= N1; K++)
	{
		J1 = Common::aINN[I1];
		ISGW = Common::aISG[I1];
		I = I1;
		for (J = 0; J < J1; J++)
		{
			AX = Common::aANUMB[I][J];
			ZX = Common::aZNUMB[I][J];
			FX = Common::aPRTHK[I][J]/XN1;
			DENST = Common::aPRDEN[I][J];
			VH = sqrt(2.13e-3 * EH / A);

			Common::dedx(&Z, &A, &ZX, &AX, &DENST, &EH, &VH, &ISGW, &DEX, &DE);
			EH = EH + DE * sign * FX;


		
			if (EH <= 0.0)
			{
				if (K <= 2)
				{
					N1 = N1 * 2;
                                        XN1 = (double)N1;
					J = -1;
					K = 0;
					DEDNEXT = 0.0;
					EH = E;

				}
				ISTAT = -1;
				break;

			}
			if (K <= 2)
			{
				DEDNEXT = DEDNEXT + DE * FX;
			}

		}
		if (K == 1)
		{
			DED1ST = DEDNEXT;
			DEDNEXT = 0.0;
		}
		if (K == 2)
		{
			DDD = DED1ST - DEDNEXT;
			if (DDD < 0.0)
			{
				DDD = -DDD;
			}

			DDS = DED1ST + DEDNEXT;
			DDR = DDD / DDS;
			if (DDR > EPS)
			{
				N1 = N1 * 2;
                                XN1 = (double)N1;
				J = -1;
				K = 0;
				DEDNEXT = 0.0;
				EH = E;


			}
		}
	}

	ISTAT = 0;
	DEE = EH - E;

	*pI1 = I1;
	*pISTAT = ISTAT;
	*psign = sign;
	*pXN1 = XN1;
	*pEPS = EPS;
	*pA = A;
	*pZ = Z;
	*pE = E;
	*pDEE = DEE;
	return 0;
}

double Common::dedx(double * pZ1, double * pA1, double * pZ2, double * pA2, double * pRHO, double * pENER, double * pV, int * pIFG, double * pDEDXHI, double * pDEDXTO)
{

	double Z1, A1, Z2, A2, RHO, ENER, V, DEDXHI, DEDXTO;
	int IFG;

	Z1 = *pZ1;
	A1 = *pA1;
	Z2 = *pZ2;
	A2 = *pA2;
	RHO = *pRHO;
	ENER = *pENER;
	V = *pV;
	DEDXHI = *pDEDXHI;
	DEDXTO = *pDEDXTO;
	IFG = *pIFG;



	/*Program calculates the differential energy loss DE/DX in solid targets using a semiempirical formula
	deduced from experimental work
	
	This program is modified for gas absorbers
	REF.: K.BRAUNE,R.NOVOTNY,D.PELTE,D.HUSAR,D.SCHWALM,
	PROCEEDINGS - SPRING MEETING OF THE GERMAN PHYSICAL
	SOCIETY, VERHANDLUNGEN 4/1978
	K.BRAUNE, DIPLOM, HEIDELBERG 1979
	 	  H(Z2) IS A SUM OF FIVE GAUSSIAN FUNCTIONS.
	      A1        MASS NUMBER  - PROJECTILE
	      Z2        ATOMIC NUMBER ABSORBER
	      A1        MASS NUMBER   ABSORBER
	      RHO       DENSITY OF THE ABSORBER (GRAMM/CM**3)
	      (MEANLESS IF GAS ABSORBER )
	      ENER      ENERGY OF THE PROJECTILE (MEV)
	      V         VELOCITY OF THE PROJECTILE
	      IN MEV/(MG/CM**2)
		  Z1       ATOMIC NUMBER - PROJECTILE	*/
	double XI, A2SAV, Z2SAV, FY, G1, G2, G3, G4, G5, HZ2, Z2ZWD, FG;
	double ALEFG, GXI, SQXI, C, FG0, AL, Y, VV0, FV, AZ1, QQ, GHI, VZ1;
	double ZA, EPS, SIGMAN, DEDXNU;
	A2SAV = 0.0;
	Z2SAV = 0.0;

	if (IFG == 1)
	{
		RHO = 1.0;
	}

	XI = V * V / Z2;

	/*Absorber - Function
	G(XI) = Y(EXP) - Y(Theory) Is deduced from experimental energy loss measurements*/

	if (A2 != A2SAV && Z2 != Z2SAV)
	{
		A2SAV = A2;
		Z2SAV = Z2;

		//FY is function Y
		FY = 54721.0 * (1.0 + 5.15E-2 * sqrt(A2 / RHO) - exp(-0.23*Z2));
		
		if (IFG == 1)
		{
			FY = 54721.0 * (1.35 - exp(Z2*((-0.13 + 0.0014*Z2))));
		}
	}

	//G(XI) Is the derivation of a guassian with variable height H(Z2)
	if (Z2 <= 26.0)
	{
		G1 = 19.84 * exp(-0.17 * (Z2 - 4.25)*(Z2 - 4.25));
	}
	else if (Z2 > 26.0)
	{
		G1 = 0.000001;
	}

	if (Z2 <= 38.0)
	{
		G2 = 17.12 * exp(-0.12 * (Z2 - 11.63)* (Z2 - 11.63));
	}
	else if (Z2 > 38.0)
	{
		G2 = 0.0000001;
	}

	G3 = 7.95 * exp(-0.015 * (Z2 - 30.2) * (Z2 - 30.2));
	G4 = 5.84 * exp(-0.022 * (Z2 - 48.63) * (Z2 - 48.63));
	G5 = 7.27 * exp(-0.005 * (Z2 - 73.06) * (Z2 - 73.06));
	HZ2 = (9.0 - (G1 + G2 + G3 + G4 + G5)) * 1.32e-5;
	Z2ZWD = cbrt(Z2) * cbrt(Z2);

	//Multiplication factors of G(XI)

	FG = (1.2e-4)*Z2*Z2 + ((2.49e-2)* A2 / RHO);
	if (IFG == 1)
	{
		FG = 1.3 / (1.0 + exp(3.0 - (Z2 / 5.0)));
	}

	ALEFG = log((2.7e-5) / FG);

	//Calculation of G(XI)

	GXI = 0.0;

	if (XI >= 1.0e-9 && XI <= 5.0e-4)
	{
		SQXI = sqrt(XI);
		C = (2.0 / Z2) * (SQXI / (1.0 + 1.0e4 * SQXI));

		if (IFG == 1)
		{
			C = C / 2.0;
		}

		FG0 = 1.0 / (1.0 + (XI * 10000.0)*(XI * 10000.0)*(XI * 10000.0));
		AL = log(XI) - ALEFG;
		GXI = (C - HZ2 * AL * exp(-0.32 * AL *AL))*FG0;
	}

	// Calculation of Y(XI)

	Y = 3.3e-4 * log(1.0 + (XI*FY)) + GXI;

	//Energy Loss of heavy ions
	//Effective charge

	VV0 = V * 137.0;
	FV = 1.0;

	if (V >= 0.62)
	{
		FV = 1.0 - exp(-VV0);
	}

	AZ1 = log(1.035 - 0.4 * exp(-0.16*Z1));

	QQ = V / (pow(Z1, 0.509));
	GHI = Z1;
	VZ1 = (-116.79 - 3350.4 * QQ) * QQ;
	if (VZ1 > -85.2)
	{
		GHI = Z1 * (1.0 - exp( VZ1 ));
	}
	if (Z1 > 2.0)
	{
		GHI = Z1 * (1.0 - exp(FV * AZ1 - 0.879 * (VV0 / pow(Z1, 0.65))));
	}

	//Effective charge for protons and alpha particles

	//Electronic energy loss DEDXHI
	DEDXHI = GHI * GHI * Z2 * Y / (A2 * V * V);
	//Nuclear energy loss DEDXNU
	ZA = sqrt(cbrt(Z1) * cbrt(Z1) + Z2ZWD);

	EPS = 3.25e4 * A2 * ENER / (Z1 * Z2 * (A1 + A2) * ZA);
	SIGMAN = 1.7 * sqrt(EPS) * log(EPS + 2.1718282) / (1.0 + 6.8 * EPS + 3.4 * sqrt(EPS) *  sqrt(EPS) *  sqrt(EPS));

	DEDXNU = SIGMAN * 5.105 * Z1 * Z2 * A1 / (ZA * A2 * (A1 + A2));
	//Total energy loss
	DEDXTO = DEDXHI + DEDXNU;

	*pZ1 = Z1;
	*pA1 = A1;
	*pZ2 = Z2;
	*pA2 = A2;
	*pRHO = RHO;
	*pENER = ENER;
	*pV = V;
	*pDEDXHI = DEDXHI;
	*pDEDXTO = DEDXTO;
	*pIFG = IFG;

	return 0;
}

void defprta( double energy, double * elost, double depth)
{
	Common Target;

	Target.aPRS[0] = 10.0;
	Target.aXLN[0] = 10.0;
	Target.aCONCN[0][0] = 1.0;
	Target.aCONCN[0][1] = 1.0;


	Target.aISG[0] = 0;
	Target.aINN[0] = 2;
	Target.aDEN[0] = 3.18;
	Target.aZNUMB[0][0] = 20;
	Target.aZNUMB[0][1] = 9;
	Target.aANUMB[0][0] = 40;
	Target.aANUMB[0][1] = 19;
	Target.aELNUM[0][0] = 1;
	Target.aELNUM[0][1] = 2;
	Target.loste[0] = 0;
	Target.aTHK[0] = 0.5*depth;

	Target.desorb(1, 2.0, 3.0, energy);
        
        *elost = Target.loste[0];


}
void defprtb( double energy, double * elost, double depth, double angle )
{
	Common Target;

	Target.aPRS[0] = 10.0;
	Target.aXLN[0] = 10.0;
	Target.aCONCN[0][0] = 1.0;
	Target.aCONCN[0][1] = 1.0;


	Target.aISG[0] = 0;
	Target.aINN[0] = 2;
	Target.aDEN[0] = 3.18;
	Target.aZNUMB[0][0] = 20;
	Target.aZNUMB[0][1] = 9;
	Target.aANUMB[0][0] = 40;
	Target.aANUMB[0][1] = 19;
	Target.aELNUM[0][0] = 1;
	Target.aELNUM[0][1] = 2;
	Target.loste[0] = 0;
	Target.aTHK[0] = 0.5*(1.0-depth)/cos(angle);

	Target.desorb(1, 1.0, 3.0, energy);
        
        *elost = Target.loste[0];


}

class reaction {
public:

   reaction(double Z1, double A1, double Me1, double Z2, double A2, double Me2, double Z3, double A3, double Me3, double Z4, double A4, double Me4);
   ~reaction();

   void SetMass();

   void GenerateState();
   void GetExcitedE();
   void GenerateReaction(int j);
   void SSnaps();
   double B_CALC(double * BRSL1, double * BZSL1,double MA,double ERS,double ZZ,double CC,double LL,double RR );
   vector<double> CreateCOMVelocityFrame1Vector();
   vector<double> CreateCOMHe5Vector();
   vector<double> CreateHe4VelocityLabVector(vector<double> &vec);

   double LabthetaHe4Final();
   double LabthetaNFinal();

   double m1, m2, m3, m4,a1,a2,a3,a4,me1,me2,me3,me4, z1,z2,z3,z4, current, m_FinalZOut;
   double m_Q, m_ExciteM, m_QExcite, m_CenterOfMassTotalE, m_CenterOfMassEOut, m_CenterOfMassERec, m_RandCOMThetaOut; 
   double m_RandCOMThetaRec, m_COMVelocityOut, m_COMVelocityRec, m_COMVelocityBeam, m_COMVelocity, m_KOut, m_KRec;
   double m_LabEnergyOut, m_LabEnergyRec, m_LabVelocityOut, m_LabVelocityRec, m_LabthetaRec;
   double m_LabthetaOut, m_RandPhiOut, m_RandPhiRec;
   double T1, m_ExcitedEspread, m_InitExcitedE, m_Width;
	
private:
   double Pi;
   static constexpr double hbar = 6.58217e-22;       //Units of MeV*s
   static constexpr double lhtspd = 2.99792e08;
   vector<double> VelocityCOMHe5Vector;
   vector<double> VelocityVectorCOMFrame1;
   vector<double> VelocityHe4LabVector;
};
//   reaction *Reaction = new reaction(11, 20.174, 9, 11.3476, 8, 4.94167, 12, 25.076);
     // m2 is the light initial particle m3 is the detected partile m4 is the recoil
  reaction::reaction(double Z1, double A1, double Me1, double Z2, double A2, double Me2, double Z3, double A3, double Me3, double Z4, double A4, double Me4) 
     :z1(Z1), a1(A1), me1(Me1), z2(Z2), a2(A2), me2(Me2), z3(Z3), a3(A3), me3(Me3), z4(Z4), a4(A4), me4(Me4)
{
     VelocityVectorCOMFrame1.push_back(0.0);
     VelocityVectorCOMFrame1.push_back(0.0);
     VelocityVectorCOMFrame1.push_back(0.0);
}

  reaction::~reaction(){}

  void reaction::SetMass() {
     m1 = a1*931.494027+me1;                               //Units of MeV/c^2
     m2 = a2*931.494027+me2;                               //Units of MeV/c^2
     m3 = a3*931.494027+me3;                               //Units of MeV/c^2
     m4 = a4*931.494027+me4;                               //Units of MeV/c^2
     Pi = 4.0 * atan(1.0);
}
   
   void reaction::GenerateState() {
     //All in units of MeV
     double a = rando->Uniform(0.0,1.0);
     double b = (1.0/9.0);
//     a=0.1*b;
     if (0.0 <= a && a < b){
       m_InitExcitedE = 0;
       m_Width = hbar/(2*21.3e-3);
     }
     else if (b <= a && a < 2*b){
       m_InitExcitedE = 0.23827;
       m_Width = 0.0; 
     }
     else if (2*b <= a && a < 3*b){
       m_InitExcitedE = 1.536;
       m_Width = 0.0;
     }
     else if(3*b <= a && a < 4*b) {
       m_InitExcitedE = 4.0329;
       m_Width = 0.0;
     }
     else if(4*b <= a && a < 5*b) {
       m_InitExcitedE = 6.290;
       m_Width = 0.0;
     }
     else if (b <= a && a < 6*b){
       m_InitExcitedE = 6.416;
       m_Width = 0.0; 
     }
     else if (2*b <= a && a < 7*b){
       m_InitExcitedE = 6.440;
       m_Width = 0.0;
     }
     else if(3*b <= a && a < 8*b) {
       m_InitExcitedE = 6.459;
       m_Width = 0.0;
     }
     else if(4*b <= a && a < 9*b) {
       m_InitExcitedE = 6.741;
       m_Width = 0.0;
     }

     GetExcitedE();
   }

   void reaction::GetExcitedE() {  
     do{
     m_ExcitedEspread = rando->BreitWigner(m_InitExcitedE,m_Width);  
     }while(m_ExcitedEspread < 0 || m_ExcitedEspread > 14);
   }

void reaction::GenerateReaction(int T) {

     GenerateState();    

	double *kinematics = new double[9];
	
	double temp, ExciteE=0;
	TVector3 productDirection;
	
	double RandNum = rando->Uniform(0.0,1.0);
	
	double pi = 4*atan(1.0);	
	
	double T1 = T;	
	
  	m_QExcite = (m2 + m1 - m3 - m4-m_InitExcitedE);   

	double qValue = m1 + m2 - m3 - m4;

//-------------------------------------------------
        double depth = rando->Uniform(0.0,1.0);
        double energy, eloss;
//-------------------------------------------------

        energy = T1;
        eloss = 0.0;
//        cout << depth << " " << energy << endl;
        
        defprta(energy, &eloss, depth);
//        cout << depth << " = depth, eloss = " << eloss << endl;

        T1 = T1 - eloss;
        
//-------------------------------------------------


	m_CenterOfMassTotalE = (m1*T1)/(m1+m2);
	m_COMVelocityBeam= m2*sqrt(2*T1/m2)/(m1+m2)*lhtspd;
	
//	cout << m_COMVelocityBeam << "   " << m2 << "  " << m1 <<endl;
	
	m_CenterOfMassEOut = (m4*(m_CenterOfMassTotalE+m_QExcite))/(m3+m4);
   	m_CenterOfMassERec = (m3*(m_CenterOfMassTotalE+m_QExcite))/(m3+m4);   
//		cout << m_CenterOfMassTotalE << "  "  << m_CenterOfMassEOut << endl;	
	//Create a random theta
	//double m_RandCOMThetaOut = pi*rando->Uniform(0.0,0.9999);

	m_RandCOMThetaOut = acos( 2.0*rando->Uniform(0.1,1.0) - 1.0 )*70/180;	//lab angles < 90
//  	m_RandCOMThetaOut = 30.4*pi/180;
	m_RandCOMThetaRec = pi - m_RandCOMThetaOut;   	
	 m_RandPhiOut = 2*acos (2.0*rando->Uniform(0.0,1.0)-1.0);
	 m_RandPhiRec = abs(pi - m_RandPhiOut);
//	 cout << m_RandCOMThetaOut << "  " << m_RandPhiOut << endl;

	m_KOut = sqrt((m2*m3*m_CenterOfMassTotalE)/(m1*m4*(m_CenterOfMassTotalE+m_QExcite)));
  	m_KRec = sqrt((m2*m4*m_CenterOfMassTotalE)/(m1*m3*(m_CenterOfMassTotalE+m_QExcite)));   	
	m_LabEnergyOut = m_CenterOfMassEOut*(1. + m_KOut*m_KOut + 2.*m_KOut*cos(m_RandCOMThetaOut));
   	m_LabEnergyRec = m_CenterOfMassERec*(1. + m_KRec*m_KRec + 2.*m_KRec*cos(m_RandCOMThetaRec));   	

//cout << m_QExcite << "  " << m_CenterOfMassTotalE << "  " <<qValue<<"  "<<m4<<endl;
//     cout << m_KOut << "  " << m_KRec << endl;
  	
	 m_LabthetaOut = atan2((sin(m_RandCOMThetaOut)),(cos(m_RandCOMThetaOut)+m_KOut));
   	 m_LabthetaRec = atan2((sin(m_RandCOMThetaRec)),(cos(m_RandCOMThetaRec)+m_KRec));

//-----------------------------------------------------------------------------

        energy = m_LabEnergyOut;
        
        defprtb( energy, &eloss, depth, m_LabthetaOut);
        
//        cout << "energy = " << energy << " eloss = " << eloss << endl;

        m_LabEnergyOut = m_LabEnergyOut - eloss;

//        cout << m_LabEnergyOut << " = Final Energy, Theta = " << m_LabthetaOut << endl;
//----------------------------------------------------------------------------- 	

	
   double tempPhi =  rando->Uniform(0.0,1.0)* 2.0 * pi;
   double HeavyPhi; //= abs(tempPhi-pi);
   
	if( tempPhi < pi )
		HeavyPhi = pi + tempPhi;
	else if( tempPhi > pi )
		HeavyPhi = tempPhi - pi;	

   m_RandPhiOut = tempPhi;
   m_RandPhiRec = HeavyPhi;
 	
   productDirection[0] = sin(m_LabthetaOut)*cos(tempPhi);
   productDirection[1] = sin(m_LabthetaOut)*sin(tempPhi);
   productDirection[2] = cos(m_LabthetaOut);			
	
   double Theta = productDirection.Theta();	
   double productEnergy = m_LabEnergyOut;

   
   SSnaps();

	
;

	
	kinematics[0] = m_LabEnergyOut;
	kinematics[1] = m_LabthetaOut;
	kinematics[2] = tempPhi;
	kinematics[3] = m_CenterOfMassEOut;
	
	kinematics[4] = m_LabEnergyRec;
	kinematics[5] = m_LabthetaRec;		
	kinematics[6] = HeavyPhi;
	kinematics[7] = m_CenterOfMassERec;
	kinematics[8] = m_FinalZOut;	

//	cout << m_LabEnergyOut << "   " << m_LabthetaOut*57.3 << "   "   << m_FinalZOut <<endl; 
}     


  vector<double> reaction::CreateCOMVelocityFrame1Vector() {

   VelocityVectorCOMFrame1[0] = 0.;
   VelocityVectorCOMFrame1[1] = 0.;
   VelocityVectorCOMFrame1[2] = m_COMVelocity;
   return VelocityVectorCOMFrame1;
}

  vector<double> reaction::CreateCOMHe5Vector() {
     double mu1 = m3*m4/(m3+m4);
     double a = (2.*(m_QExcite + m_CenterOfMassTotalE)*mu1)/(m3*m3);
     double m_Vx = sqrt(a)*sin(m_RandCOMThetaOut);;
     double m_Vy = 0.;
     double m_Vz = sqrt(a)*cos(m_RandCOMThetaOut);
   
     VelocityCOMHe5Vector.clear();
     VelocityCOMHe5Vector.push_back(m_Vx);
     VelocityCOMHe5Vector.push_back(m_Vy);
     VelocityCOMHe5Vector.push_back(m_Vz);
     return VelocityCOMHe5Vector;
}

  vector<double> reaction::CreateHe4VelocityLabVector(vector<double> &vec) {
     CreateCOMHe5Vector();
     CreateCOMVelocityFrame1Vector();
     double m_Vx = VelocityVectorCOMFrame1[0]+VelocityCOMHe5Vector[0]+vec[0];//*m3/m4;
     double m_Vy = VelocityVectorCOMFrame1[1]+VelocityCOMHe5Vector[1]+vec[1];//*m3/m4;
     double m_Vz = VelocityVectorCOMFrame1[2]+VelocityCOMHe5Vector[2]+vec[2];//*m3/m4;
   
     VelocityHe4LabVector.clear();
     VelocityHe4LabVector.push_back(m_Vx);
     VelocityHe4LabVector.push_back(m_Vy);
     VelocityHe4LabVector.push_back(m_Vz);
     return VelocityHe4LabVector;
}

   double reaction::LabthetaHe4Final() {
     return (180./Pi)*atan2( sqrt(VelocityHe4LabVector[0]*VelocityHe4LabVector[0]+VelocityHe4LabVector[1]*VelocityHe4LabVector[1]),
          (VelocityHe4LabVector[2]) );
}

   double reaction::LabthetaNFinal() {
     return (180./Pi)*atan2( sqrt(VelocityHe4LabVector[0]*VelocityHe4LabVector[0]+VelocityHe4LabVector[1]*VelocityHe4LabVector[1]),
          (VelocityHe4LabVector[2]) );
}

   void reaction::SSnaps() {

       int XYZ;
       double ZYX, TEMP;

       double  QO,QI,ZO,AMU,EO,THETA,RN;
       int NR, SETCOUNT;
 
       double    PAR_SOLENOID;
       int SOLNAME1,SOLNAME2;
       double    S1,LL,LLS,ER,ERS,MA1,MA2,I1,I2  ;
       double    ZTGT,Z1,Z2,SEPAR,CC,CCC,BM1,BM2,LENGTH,BORE;
; 
       int     E_RAD;
       double    V_INNER,V_OUTER,R_INNER,R_OUTER,L_ELECTRODE,Z_E_1,Z_E_2;

       int ION, ABSORBER, POS_ABSORBER, POS_TARGET;
       double THICKNESS_ABSORBER,POSITION_ABSORBER,AMU_ABSORBER;
       double M1,M2,RHO,ATRHO,VFERMI,LFCTR,EPSIL, SE, SN, VELSQR;
       double POSITION_TARGET;
       double PCOEF[8];
       int MM1;      

       double  DS_S1, DS_LL, DS_ER, DS_ERS;
       double  DS_MA, DS_Z1, DS_CC, DS_BM, BM3, DS_LENGTH, DS_BORE;       

       double TORR;            
       double    N,M,E,T,ANG3,V0,VX,VY,VZ,V1,V2,V3,RR,RRS,INRR,INRRS,EL;
       double    X0,Y0,Z0,X,Y,Z,GX,GY,GZ,HGX,HGY,HGZ,F7,F8,F9,HF7,HF8,HF9;
       double    XX,YY,ZZ,XXS,YYS,AX,AY,AZ,T0;
       double    ZND,DZ,HDZ,DT,DTS,CQM,XP,YP;
       double    BZS,BRS,BR,BZ,BX,BY,E_FIELD,EX,EY,EZ,EFinal;  
       int I,STEP,STEPNUM,J, J_MAX;
       double    STEPNOM;
       
       
       current = 97;
       V0=0.0;
       HDZ=0.0;

      NR=1; //number of rays
      QO=z3;
      QI=z3;
      ZO=z3;
      AMU=a3+me3/931.494;
      EO=m_LabEnergyOut;//+45.0/1000.0;
      THETA=m_LabthetaOut;

//	cout << "EO   " << EO << "  " << m_LabEnergyOut<< "  " << THETA << "  " << m_LabthetaOut << endl;
//      ZTGT = 1.021;
      ZTGT = 0.35;
      POSITION_TARGET = 5.443;

      DZ=.0010000;
      HDZ=.0010000;
      ZND=1.5;

      STEPNOM = ZND/DZ;
      STEPNUM = int(STEPNOM);
      POS_TARGET = POSITION_TARGET/DZ;
      SOLNAME1 = 1;


     I1=current;
     S1=0.62;
     ER=0.2123;
     MA1=3.636*I1/100;
     LENGTH = 1.00;
     BORE = 0.15;
     Z1 = ZTGT - 0.5*S1;  

//    * effective length of solenoid    *

      LL = S1;         
//    * entrance position of solenoid    *

      CC = Z1;                
      CCC = Z2;

//    * generally an "S" indicates a squared variable     *

      LLS=LL*LL;
      ERS=ER*ER;
//    * BM = Max B at the center of the solenoid      *
      
      BM1=MA1*LL/sqrt(ERS+LLS/4);
      BM2=MA2*LL/sqrt(ERS+LLS/4);
      BM3=DS_MA;

/*************************************
 MAIN LOOP - done for each particle |
************************************/

         EL=ZO;
         N=QO;
         M=AMU;
         ANG3=THETA;
         E=EO;
         T=ANG3; //   EMERGING ANGLE
      
         V0=1.389e+07*sqrt(E/M);

         double randx = rando->Uniform(0.0,1.0);
         double randy = rando->Uniform(0.0,1.0);

         X0=randx * 0.0015;
         Y0=randy * 0.0015;
         Z0=0.0;
         T0=0.0;
         X=X0;
         Y=Y0;
         Z=Z0;
         VX=V0*sin(T);
         VY=0.0;
         VZ=V0*cos(T);
         GX=0.0;
         GY=0.0;
         GZ=0.0;
         HGX=0.5*GX;
         HGY=0.5*GY;
         HGZ=0.5*GZ;
         F7=0.0;
         F8=0.0;
         F9=0.0;
         HF7=0.5*F7;
         HF8=0.5*F8;
         HF9=0.5*F9;
          
/************************************************************************
!   CALCULATIONS OF POSITIONS AND FIELDS STARTS HERE - THIS HALF STEP   *
!   TECHNIQUE ALLOWS LARGE STEPS - Note that all the field calculations *
!   are done at the halfway point through the interval.                 *
!************************************************************************/
//	cout << ANG3*57.3 << "   " << AMU << "   " << EO << endl;
         for ( STEP=1;STEP<=STEPNUM; STEP++)  {
//	    cout << STEP << endl;
            V1=VX+HF7;
            V2=VY+HF8;
            V3=VZ+HF9;
            DT=DZ/V3;
            DTS=DT*DT;
            XX=X+HGX;
            YY=Y+HGY;
            ZZ=Z+HDZ;
            XXS=XX*XX;
            YYS=YY*YY;
            RRS=XXS+YYS;
            RR=sqrt(RRS);
            INRRS=1/(RRS+0.00000001);
            INRR=1/(RR+0.00000001);
           
/**************************************************
!  This next section uses rationalized B fields.  *
!**************************************************/

            BZS=0.;
            BRS=0.;
            B_CALC(&BRS,&BZS, MA1,ERS,ZZ,CC,LL,RR);
//	    cout<< "" <<endl;;//MA1 << "  " << ERS << "  " << ZZ << "  " << CC << "  " << LL << "  " << RR <<endl;
//	    cout << STEP << "   " << BZS << " after " << BRS << endl;
            BR=BRS;
//	    cout << "B   " << BZS << "   " << BRS << endl;
            BZ=BZS;
            BX=BR*INRR*XX;
            BY=BR*INRR*YY;
            CQM=9.648e+07*(N/M);


/*******************************************************
!  calculate the cross products V x B and also include the electrostatic field            *
!*******************************************************/
            AX=CQM*(V2*BZ-V3*BY+EX); 
            AY=CQM*(V3*BX-V1*BZ+EY);
            AZ=CQM*(V1*BY-V2*BX+EZ);

      
/**********************************************************************
!  calculate the X and Y displacements and new X and Y positions       
!**********************************************************************/
            GX=V1*DT+.5*AX*DTS;
            GY=V2*DT+.5*AY*DTS;
            X=X+GX;
            Y=Y+GY;
            Z=DZ*STEP;
            T0 = T0 + (sqrt(GX*GX+GY*GY+DZ*DZ))/(1.0e-09*V0); //Time in ns          

/**********************************************************************
!  calculate the increments of the velocity                            
!  calculate the half steps to be used in next iteration of loop       
!**********************************************************************/
            F7=AX*DT;
            F8=AY*DT;
            F9=AZ*DT;
            HF7=.5*F7;
            HF8=.5*F8;
            HF9=.5*F9;
            VX=VX+F7;
            VY=VY+F8;
            VZ=VZ+F9;
//	 cout << STEP << "   "<< RR << "   "  <<  THETA*57.3 <<"   "   << Z << endl;  
	if(RR > 0.15){
		STEP = STEPNUM; 
		m_FinalZOut = 0.0;
	} 
	if(STEP > 80 && RR < 0.02 && STEP!=STEPNUM){
//	    cout << STEP << "    " << THETA*57.3 <<"   "  <<RR << "  " <<Z <<endl;
	    STEP = STEPNUM;
	    m_FinalZOut = 100*Z;	//puts Z in cm
	 }
	
//	if(STEP % 50 == 0) cout << STEP << "   " << RR << "  " << Z << endl;
   }
//	cout << "EO   " << EO << "  " << THETA*57.3 <<endl;// "  " << BZS << endl;
 }   
      double reaction::B_CALC( double * BRSL, double * BZSL,double MA,double ERS,double ZZ,double CC,double LL,double RR){

      /******************************************************************
      ! calculate the radial and axial components of the magnetic field  
      ! using rationalized expressions from Liu's original code.         
      !******************************************************************/

          double MAERS,ZC,ZCS,ZL,ZLS,SCS,SC,SLS,SL,SC3,SC5,SC7;
          double SL3,SL5,SL7,HR,HRS,HR3,BZZ,B1,B3,BRR;
          MAERS=MA*ERS;
          ZC=ZZ-CC;
          ZCS=ZC*ZC;
          ZL=ZZ-CC-LL;
          ZLS=ZL*ZL;
          SCS=ERS+ZCS;
          SC=sqrt(SCS);
          SLS=ERS+ZLS;
          SL=sqrt(SLS);
          SC3=SCS*SC;
          SC5=SC3*SCS;
          SC7=SC5*SCS;
          SL3=SLS*SL;
          SL5=SL3*SLS;
          SL7=SL5*SLS;
          HR=.5*RR;
          HRS=HR*HR;
          HR3=HRS*HR;
          BZZ=MA*(ZC/SC-ZL/SL) - HRS*3*MAERS*(-ZC/SC5+ZL/SL5);
          B1=-HR*MAERS*(1/SC3-1/SL3);
          B3=HR3/2*3*MAERS*((4*ZCS-ERS)/SC7-(4*ZLS-ERS)/SL7);
          BRR=B1+B3;
          *BZSL=*BZSL + BZZ;
          *BRSL=*BRSL + BRR;
//	cout <<*BZSL << "  before  " <<*BRSL << endl;
	return 2.0;
/***************************************************************************
!  Differences between the first and second solenoids when calling          
!  this subroutine.:                                                        
!  - The variable on the left is the one found inside the subroutine.       
!  - The two variables on the right are the ones that are sent to the       
!     subroutine from the main program.                                     
!                                                                           
!                              1st Solenoid        2nd Solenoid             
!            MA            =       MA1                 MA2                  
!            CC            =       CC                  CCC                  
!            PAR_SOLENOID  =        1              PAR_SOLENOID             
!                                                                           
! Before the B_CALC subroutine is called for the first solenoid, make       
! sure that BZSL and BRSL have been initialzed to zero. Also when           
! calling the subroutine for the first solenoid, the variable PAR_SOLENOID  
! receives zero as its value.  PAR_SOLENOID specifies whether the two       
! solenoids are in parallel or antiparallel configurations. Since solenoids 
! are always focusing, the absolute orientations of the two magnetic fields 
! is unimportant, only their relative orientations. Note that having the    
! two solenoids parallel will reduce the strength of the fringe fields      
! in the space between the them.  Also, there is slightly greater focusing  
! power when the two solenoids are parallel.                                
!***************************************************************************/
}
             


int main()
{   

      // m2 is the beam. m3 is the detected partile m4 is the recoil  
   //19F(3He,t)19Ne
  


   double COMThetaOut, LabThetaOut, LabEnergyOut, COMThetaRec, LabThetaRec, LabEnergyRec;
   double COMEnergyOut, COMEnergyRec, FinalZOut;
   double KN, RandPhiN, period;
   double ExcitedE, BeamEnergy, ExcitedEDist, QVal;
   double PhiOut, PhiRec,heliosTheta, heliosZ;
   double Bfield = 6;
   double OutMass = 3;
   double charge = 1;
   reaction *Reaction = new reaction( 9,19, -1.487,2, 3, 14.931, 1,OutMass,14.950, 10, 19, 1.751);
   double Vcm;
   
   Reaction->SetMass();
   vector<double> COMHe4Frame1Vec;
   vector<double> COMNFrame1Vec;
       
// Constants and Conversions
	
   double Pi = 4.0 * atan(1.0);
   double MeVMass = 931.502; // MeV/c^2
   double ElCharge = 1.6e-19; // Coulombs
   double SpeedOfLight = 299792458.0 // m/s
   double SecondToNS = 1000000000; // ns/s
   double MeVToJoules = 1.6e-13; // J/MeV
   
   

	int InputBeamEnergy; //Units of MeV. total beam energy
	
	cout << "Input the total beam energy of 3He in units of MeV:";
	cin >> InputBeamEnergy;

	char buffer[124];
   sprintf(buffer,"results/19F_3Het_%dMeV_res.root",InputBeamEnergy);

   TFile *f = new TFile(buffer,"RECREATE");
   
   TTree *Tree = new TTree("react","MonteCarlo Simulation");   

	default_random_engine generator;

	
   Tree->Branch("ExcitedE",&ExcitedE,"ExcitedE/D");
   Tree->Branch("QVal",&QVal,"QVal/D");
   Tree->Branch("ExcitedEDist",&ExcitedEDist,"ExcitedEDist/D");
   Tree->Branch("BeamEnergy",&BeamEnergy,"BeamEnergy/D");
   Tree->Branch("PhiOut", &PhiOut, "PhiOut/D");
   Tree->Branch("PhiRec", &PhiRec, "PhiRec/D");
   
   Tree->Branch("COMThetaOut",&COMThetaOut,"COMThetaOut/D");
   Tree->Branch("COMEnergyOut",&COMEnergyOut,"COMEnergyOut/D");   
   Tree->Branch("LabThetaOut",&LabThetaOut,"LabThetaOut/D");
   Tree->Branch("LabEnergyOut",&LabEnergyOut,"LabEnergyOut/D");
   Tree->Branch("FinaZOut",&FinalZOut,"FinalZOut/D");
   Tree->Branch("heliosTheta",&heliosTheta,"heliosTheta/D");
      
   Tree->Branch("COMThetaRec",&COMThetaRec,"COMThetaRec/D");
   Tree->Branch("COMEnergyRec",&COMEnergyRec,"COMEnergyRec/D"); 
   Tree->Branch("LabThetaRec",&LabThetaRec,"LabThetaRec/D");
   Tree->Branch("LabEnergyRec",&LabEnergyRec,"LabEnergyRec/D");
   	ofstream outf("values.dat");
   for (int i=0;i<10000;i++)
   //for (int i=0;i<5;i++)
	{
        Reaction->GenerateReaction(InputBeamEnergy);
//Fill in Tree for the reaction
        QVal = Reaction->m_QExcite;
        //if( Reaction->m_InitExcitedE ==0)       
            //cout << Reaction->m_Q << endl;
        ExcitedE = Reaction -> m_InitExcitedE;
        ExcitedEDist = Reaction -> m_ExcitedEspread;
        BeamEnergy = Reaction->T1;
		period = 2*Pi*OutMass*931.502*0.00000000000016/charge/ElCharge/6/SpeedOfLight/SpeedOfLight*1000000000; //in ns
		Vcm = Reaction->m_COMVelocityBeam;
		COMThetaOut = (180./Pi)*(Reaction->m_RandCOMThetaOut);
		COMEnergyOut = Reaction->m_CenterOfMassEOut;
        LabThetaOut = (180./Pi)*(Reaction->m_LabthetaOut);
        LabEnergyOut = Reaction->m_LabEnergyOut;
        FinalZOut = Reaction->m_FinalZOut;
//	cout << FinalZOut <<endl;
	
		heliosTheta = acos(1/(2*Pi)*(charge*1.609e-19*Bfield*FinalZOut-2*Pi*OutMass*Vcm)/sqrt(2*OutMass*LabEnergyOut+OutMass*OutMass*Vcm*Vcm - OutMass*Vcm*charge*1.602e-19*Bfield*FinalZOut/Pi));
		heliosTheta = 180/Pi*acos(FinalZOut*sqrt(931.502*OutMass/(2*LabEnergyOut))/(SpeedOfLight*period/1000000000*100));
		heliosZ=sqrt(2*LabEnergyOut/OutMass/931.502)*cos(LabThetaOut*3.14159/180)*SpeedOfLight*period*0.000000001*100;
	
	//normal_distribution<double> distribution(LabEnergyOut,45.0/1000);	
	//LabEnergyOut = distribution(generator);	
		 
		COMThetaRec = (180./Pi)*(Reaction->m_RandCOMThetaRec);		  		  
		LabThetaRec = (180./Pi)*(Reaction->m_LabthetaRec);
		LabEnergyRec = Reaction->m_LabEnergyRec;
		COMEnergyRec = Reaction->m_CenterOfMassEOut;
	
	PhiOut = Reaction->m_RandPhiOut;
	PhiRec = Reaction->m_RandPhiRec;
	
	Tree->Fill();
        {
        outf << LabEnergyOut << " " << FinalZOut << " " << COMThetaOut << " " << COMEnergyOut << endl;
        }
	}
        
   f->Write();
   
   delete Reaction;

   delete rando;

   return 0;
}   
