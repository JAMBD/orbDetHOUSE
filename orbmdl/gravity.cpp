
#include "coordTrans.hpp"
#include "constants.hpp"
#include "satRefSys.hpp"
#include "gravity.hpp"
#include "jplEph.hpp"
#include "sofa.hpp"

#include <iostream>
#include <vector>
// typedef std::ostream Trace;

Vector3d CalcPolarAngles(Vector3d mVec)
{
	/* Norm of vector
	 */
	double mR = mVec.norm();

	/* Azimuth of vector
	 */
	double mPhi;
	if ((mVec(0) == 0) && (mVec(1) == 0))
	{
		mPhi = 0;
	}
	else
	{
		mPhi = atan2(mVec(1), mVec(0));
	}

	if (mPhi < 0)
	{
		mPhi += 2 * PI;
	}

	double rho = sqrt(mVec(0) * mVec(0) + mVec(1) * mVec(1)); // Length of projection in x - y - plane:

	/* Altitude of vector
	 */
	double mTheta;
	if ((mVec(2) == 0) && (rho == 0))
	{
		mTheta = 0;
	}
	else
	{
		mTheta = atan2(mVec(2), rho);
	}

	Vector3d vecRAE = Vector3d::Zero();
	vecRAE(0) = mR;
	vecRAE(1) = mPhi;
	vecRAE(2) = mTheta;

	return vecRAE;
}

/* Calculate normalized Legendre polynomial values
 *
 */
const std::pair<MatrixXd, MatrixXd> Legendre(
	int m,		///< Maximum degree
	int n,		///< Maximum order
	double phi) ///< Geocentric latitude in radian
{
	const double cos_phi = cos(phi);
	const double sin_phi = sin(phi);

	static std::vector<double> sqrt_table_a;
	static std::vector<double> sqrt_table_b;
	static std::vector<double> sqrt_table_c;
	static std::vector<double> sqrt_table_d;
	static int cache_m = 0;
	if (cache_m < m)
	{
		cache_m = m;
		sqrt_table_a.resize(m);
		sqrt_table_a[0] = sqrt(3);
		for (int i=1; i < m; i++)
		{
			sqrt_table_a[i] = sqrt(1.0 + (1.0 / (2.0 * (i + 1.0))));
		}
		sqrt_table_b.resize(m);
		for (int i=0; i < m; i++)
		{
			sqrt_table_b[i] = sqrt(2.0 * i + 3.0);
		}
		sqrt_table_c.resize((m+1) * (n+1));
		sqrt_table_d.resize((m+1) * (n+1));
		
		for (int j=0; j<=n; j++)
		{
			const int offset = j * (cache_m + 1);
			for (int i = j+2; i <= m; i++)
			{
				sqrt_table_c[offset + i] = sqrt((2.0 * i + 1.0) / ((i - j) * (i + j)));
				sqrt_table_d[offset + i] = sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0));
			}
		}
	}

	MatrixXd pnm = MatrixXd::Zero(m + 1, n + 1);
	// Copy faster than a Zero initializer;
	MatrixXd dpnm = pnm;
	pnm.coeffRef(0, 0) = 1;
	dpnm.coeffRef(0, 0) = 0.0;
	double prev_d = 0.0;
	double prev = 1.0;
	for (int i = 0; i < m; i++)
	{
		/* horizontal first step coefficients
		 */
		pnm.coeffRef(i + 1, i) = sqrt_table_b[i] * sin_phi * prev;

		/* horizontal first step coefficients
		 */
		dpnm.coeffRef(i + 1, i) = sqrt_table_b[i] * ((cos_phi * prev) + (sin_phi * prev_d));

		/* diagonal coefficients
		 */
		prev_d = sqrt_table_a[i] * (cos_phi * prev_d - sin_phi * prev);
		dpnm.coeffRef(i+1, i+1) = prev_d; 

		/* diagonal coefficients
		 */
		prev *= sqrt_table_a[i] * cos_phi;
		pnm.coeffRef(i + 1, i + 1) = prev;
	}

	/* horizontal second step coefficients
	 */
	for (int j=0; j<=n; j++)
	{
		const int offset = j * (cache_m + 1);
		double past_1 = pnm.coeff(j+1,j);
		double past_2 = pnm.coeff(j,j);

		double past_d_1 = dpnm.coeff(j+1,j);
		double past_d_2 = dpnm.coeff(j,j);

		for (int i = j+2; i <= m; i++)
		{
			const double tmp_a = sqrt_table_c[offset + i];
			const double tmp_b = sqrt_table_d[offset + i];
			const double tmp_c = sqrt_table_b[i-2];
			const double sqrt_sin_phi = tmp_c * sin_phi;
			const double next_value = tmp_a * (sqrt_sin_phi * past_1 - tmp_b* past_2);
			const double next_value_d = tmp_a * ((sqrt_sin_phi * past_d_1) + tmp_c * cos_phi * past_1 - tmp_b * past_d_2);

			pnm.coeffRef(i, j) = next_value; 
			past_2 = past_1;
			past_1 = next_value;

			dpnm.coeffRef(i, j) = next_value_d;
			past_d_2 = past_d_1;
			past_d_1 = next_value_d;
		}
	}


	return {std::move(pnm), std::move(dpnm)};
}

GravityModel::GravityModel(
	EarthGravMdlOpt gravMdlOpt,
	EarthTidesOpt tidesOpt,
	EGMCoef egmCoef,
	void *pJPLEph)
{
	mGravModelType = gravMdlOpt.earthGravMdl;
	mEarthGravAccDeg = gravMdlOpt.earthGravAccDeg;
	mEarthGravSTMDeg = gravMdlOpt.earthGravSTMDeg;
	mEarthTidesOpt = tidesOpt;
	mEGMCoef = egmCoef;
	mJPLEphemeris = new JPLEphemeris(pJPLEph);
}

void GravityModel::solidEarthTidesCorrection(
	double mjdUTC,				///< UTC in modified julian date format
	EGMCoef &egmCoef,			///< Struct of Earth gravity coefficients
	IERS iersIns,				///< Instance of IERS class
	const Vector3d &vecRAESun,	///< Rho, azimuth and altitude information of Sun
	const Vector3d &vecRAEMoon) ///< Rho, azimuth and altitude information of Moon
{
	/* rho, azimuthi and elevation (altitude) of Sun
	 */
	double rhoSun = vecRAESun(0);
	double azSun = vecRAESun(1);
	double elSun = vecRAESun(2);

	/* rho, azimuthi and elevation (altitude) of Sun
	 */
	double rhoMoon = vecRAEMoon(0);
	double azMoon = vecRAEMoon(1);
	double elMoon = vecRAEMoon(2);

	/* Effect of Solid Earth Tides(elastic Earth)
	 * For dC21and dS21
	 * The coefficients we choose are in - phase(ip) amplitudes and out - of - phase amplitudes of the
	 * corrections for frequency dependence, and multipliers of the Delaunay variables
	 * Refer to Table 6.5a in IERS2010
	 */
	double coeff0[48][7] = {
		{2, 0, 2, 0, 2, -0.1, 0},
		{0, 0, 2, 2, 2, -0.1, 0},
		{1, 0, 2, 0, 1, -0.1, 0},
		{1, 0, 2, 0, 2, -0.7, 0.1},
		{-1, 0, 2, 2, 2, -0.1, 0},
		{0, 0, 2, 0, 1, -1.3, 0.1},
		{0, 0, 2, 0, 2, -6.8, 0.6},
		{0, 0, 0, 2, 0, 0.1, 0},
		{1, 0, 2, -2, 2, 0.1, 0},
		{-1, 0, 2, 0, 1, 0.1, 0},
		{-1, 0, 2, 0, 2, 0.4, 0},
		{1, 0, 0, 0, 0, 1.3, -0.1},
		{1, 0, 0, 0, 1, 0.3, 0},
		{-1, 0, 0, 2, 0, 0.3, 0},
		{-1, 0, 0, 2, 1, 0.1, 0},
		{0, 1, 2, -2, 2, -1.9, 0.1},
		{0, 0, 2, -2, 1, 0.5, 0},
		{0, 0, 2, -2, 2, -43.4, 2.9},
		{0, -1, 2, -2, 2, 0.6, 0},
		{0, 1, 0, 0, 0, 1.6, -0.1},
		{-2, 0, 2, 0, 1, 0.1, 0},
		{0, 0, 0, 0, -2, 0.1, 0},
		{0, 0, 0, 0, -1, -8.8, 0.5},
		{0, 0, 0, 0, 0, 470.9, -30.2},
		{0, 0, 0, 0, 1, 68.1, -4.6},
		{0, 0, 0, 0, 2, -1.6, 0.1},
		{-1, 0, 0, 1, 0, 0.1, 0},
		{0, -1, 0, 0, -1, -0.1, 0},
		{0, -1, 0, 0, 0, -20.6, -0.3},
		{0, 1, -2, 2, -2, 0.3, 0},
		{0, -1, 0, 0, 1, -0.3, 0},
		{-2, 0, 0, 2, 0, -0.2, 0},
		{-2, 0, 0, 2, 1, -0.1, 0},
		{0, 0, -2, 2, -2, -5.0, 0.3},
		{0, 0, -2, 2, -1, 0.2, 0},
		{0, -1, -2, 2, -2, -0.2, 0},
		{1, 0, 0, -2, 0, -0.5, 0},
		{1, 0, 0, -2, 1, -0.1, 0},
		{-1, 0, 0, 0, -1, 0.1, 0},
		{-1, 0, 0, 0, 0, -2.1, 0.1},
		{-1, 0, 0, 0, 1, -0.4, 0},
		{0, 0, 0, -2, 0, -0.2, 0},
		{-2, 0, 0, 0, 0, -0.1, 0},
		{0, 0, -2, 0, -2, -0.6, 0},
		{0, 0, -2, 0, -1, -0.4, 0},
		{0, 0, -2, 0, 0, -0.1, 0},
		{-1, 0, -2, 0, -2, -0.1, 0},
		{-1, 0, -2, 0, -1, -0.1, 0}};

	/* For dC22and dS22, Refer to Table 6.5c in IERS2010
	 */
	double coeff2[2][6] =
		{
			{1, 0, 2, 0, 2, -0.3},
			{0, 0, 2, 0, 2, -1.2}};

	/* STEP1 CORRECTIONS
	 */
	const auto tmp_lgM = Legendre(4, 4, elMoon);
	MatrixXd lgM = tmp_lgM.first;
	MatrixXd dlgM = tmp_lgM.second;

	const auto tmp_lgS = Legendre(4, 4, elSun);
	MatrixXd lgS = tmp_lgS.first;
	MatrixXd dlgS = tmp_lgS.second;

	double dCmn20 = (0.29525 / 5) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(2, 0) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(2, 0));
	double dCmn21 = (0.29470 / 5) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(2, 1) * cos(azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(2, 1) * cos(1 * azSun));
	double dSmn21 = (0.29470 / 5) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(2, 1) * sin(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(2, 1) * sin(1 * azSun));
	double dCmn22 = (0.29801 / 5) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(2, 2) * cos(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(2, 2) * cos(2 * azSun));
	double dSmn22 = (0.29801 / 5) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(2, 2) * sin(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(2, 2) * sin(2 * azSun));

	double dCmn30 = (0.093 / 7) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 0) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 0));
	double dCmn31 = (0.093 / 7) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 1) * cos(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 1) * cos(1 * azSun));
	double dSmn31 = (0.093 / 7) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 1) * sin(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 1) * sin(1 * azSun));
	double dCmn32 = (0.093 / 7) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 2) * cos(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 2) * cos(2 * azSun));
	double dSmn32 = (0.093 / 7) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 2) * sin(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 2) * sin(2 * azSun));
	double dCmn33 = (0.094 / 7) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 3) * cos(3 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 3) * cos(3 * azSun));
	double dSmn33 = (0.094 / 7) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 3) * sin(3 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 3) * sin(3 * azSun));

	double dCmn40 = (-0.00087 / 5) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(4, 0) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(4, 0));
	double dCmn41 = (-0.00079 / 5) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(4, 1) * cos(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(4, 1) * cos(1 * azSun));
	double dSmn41 = (-0.00079 / 5) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(4, 1) * sin(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(4, 1) * sin(1 * azSun));
	double dCmn42 = (-0.00057 / 5) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(4, 2) * cos(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(4, 2) * cos(2 * azSun));
	double dSmn42 = (-0.00057 / 5) * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(4, 2) * sin(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(4, 2) * sin(2 * azSun));

	/* STEP2 CORRECTIONS
	 */
	double mjdUT1 = iersIns.UT1_UTC(mjdUTC) / 86400.0 + mjdUTC;
	double mjdTT = iersIns.TT_UTC(mjdUTC) / 86400.0 + mjdUTC;

	double theta_g = iauGmst06(JD2MJD, mjdUT1, JD2MJD, mjdTT);

	double dC21 = 0;
	double dS21 = 0;
	for (int jj = 0; jj < 48; jj++)
	{
		dC21 += 1e-12 * coeff0[jj][5] * sin(theta_g + PI);
		dS21 += 1e-12 * coeff0[jj][5] * cos(theta_g + PI);
	}
	dCmn21 += dC21;
	dSmn21 += dS21;

	double dC22 = 0;
	double dS22 = 0;
	for (int jj = 0; jj <= 1; jj++)
	{
		dC22 += 1e-12 * coeff2[jj][5] * sin(theta_g + PI);
		dS22 += 1e-12 * coeff2[jj][5] * cos(theta_g + PI);
	}
	dCmn22 += dC22;
	dSmn22 += dS22;

	/* Treatment of the Permanent Tide(elastic Earth)
	 */
	double dC20 = 4.4228e-8 * (-0.31460) * 0.29525;
	dCmn20 -= dC20;

	/* Effect of Solid Earth Pole Tide(elastic Earth)
	 */
	dC21 = -1.290e-9 * (iersIns.x_pole(mjdUTC));
	dS21 = 1.290e-9 * (iersIns.y_pole(mjdUTC));
	dCmn21 = dCmn21 + dC21;
	dSmn21 = dSmn21 + dS21;

	egmCoef.cmn(2, 0) += dCmn20;
	egmCoef.cmn(2, 1) += dCmn21;
	egmCoef.cmn(2, 2) += dCmn22;
	egmCoef.smn(2, 1) += dSmn21;
	egmCoef.smn(2, 2) += dSmn22;

	egmCoef.cmn(3, 0) += dCmn30;
	egmCoef.cmn(3, 1) += dCmn31;
	egmCoef.cmn(3, 2) += dCmn32;
	egmCoef.cmn(3, 3) += dCmn33;
	egmCoef.smn(3, 1) += dSmn31;
	egmCoef.smn(3, 2) += dSmn32;
	egmCoef.smn(3, 3) += dSmn33;

	egmCoef.cmn(4, 0) += dCmn40;
	egmCoef.cmn(4, 1) += dCmn41;
	egmCoef.cmn(4, 2) += dCmn42;
	egmCoef.smn(4, 1) += dSmn41;
	egmCoef.smn(4, 2) += dSmn42;
}

void GravityModel::oceanTidesCorrection(
	double mjdUTC,				///< UTC in modified julian date format
	EGMCoef &egmCoef,			///< Struct of Earth gravity coefficients
	const Vector3d &vecRAESun,	///< Rho, azimuth and altitude information of Sun
	const Vector3d &vecRAEMoon) ///< Rho, azimuth and altitude information of Moon
{

	/* rho, azimuthi and elevation (altitude) of Sun
	 */
	double rhoSun = vecRAESun(0);
	double azSun = vecRAESun(1);
	double elSun = vecRAESun(2);

	/* rho, azimuthi and elevation (altitude) of Sun
	 */
	double rhoMoon = vecRAEMoon(0);
	double azMoon = vecRAEMoon(1);
	double elMoon = vecRAEMoon(2);

	const auto tmp_lgM = Legendre(6, 6, elMoon);
	const MatrixXd lgM = tmp_lgM.first;
	const MatrixXd dlgM = tmp_lgM.second;

	const auto tmp_lgS = Legendre(6, 6, elSun);
	MatrixXd lgS = tmp_lgS.first;
	MatrixXd dlgS = tmp_lgS.second;

	double dCmn20 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.3075) / 5 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(2, 0) * cos(0 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(2, 0) * cos(0 * elSun));
	double dCmn21 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.3075) / 5 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(2, 1) * cos(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(2, 1) * cos(1 * elSun));
	double dCmn22 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.3075) / 5 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(2, 2) * cos(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(2, 2) * cos(2 * elSun));

	double dCmn30 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.195) / 7 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 0) * cos(0 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 0) * cos(0 * elSun));
	double dCmn31 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.195) / 7 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 1) * cos(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 1) * cos(1 * elSun));
	double dCmn32 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.195) / 7 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 2) * cos(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 2) * cos(2 * elSun));
	double dCmn33 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.195) / 7 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 3) * cos(3 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 3) * cos(3 * elSun));

	double dCmn40 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.132) / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 5)) * lgM(4, 0) * cos(0 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 5)) * lgS(4, 0) * cos(0 * elSun));
	double dCmn41 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.132) / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 5)) * lgM(4, 1) * cos(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 5)) * lgS(4, 1) * cos(1 * elSun));
	double dCmn42 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.132) / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 5)) * lgM(4, 2) * cos(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 5)) * lgS(4, 2) * cos(2 * elSun));
	double dCmn43 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.132) / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 5)) * lgM(4, 3) * cos(3 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 5)) * lgS(4, 3) * cos(3 * elSun));
	double dCmn44 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.132) / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 5)) * lgM(4, 4) * cos(4 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 5)) * lgS(4, 4) * cos(4 * elSun));

	double dCmn50 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.1032) / 11 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 6)) * lgM(5, 0) * cos(0 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 6)) * lgS(5, 0) * cos(0 * elSun));
	double dCmn51 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.1032) / 11 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 6)) * lgM(5, 1) * cos(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 6)) * lgS(5, 1) * cos(1 * elSun));
	double dCmn52 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.1032) / 11 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 6)) * lgM(5, 2) * cos(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 6)) * lgS(5, 2) * cos(2 * elSun));
	double dCmn53 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.1032) / 11 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 6)) * lgM(5, 3) * cos(3 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 6)) * lgS(5, 3) * cos(3 * elSun));
	double dCmn54 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.1032) / 11 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 6)) * lgM(5, 4) * cos(4 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 6)) * lgS(5, 4) * cos(4 * elSun));
	double dCmn55 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.1032) / 11 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 6)) * lgM(5, 5) * cos(5 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 6)) * lgS(5, 5) * cos(5 * elSun));

	double dCmn60 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.0892) / 13 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 7)) * lgM(6, 0) * cos(0 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 7)) * lgS(6, 0));
	double dCmn61 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.0892) / 13 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 7)) * lgM(6, 1) * cos(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 7)) * lgS(6, 1) * cos(1.0 * elSun));
	double dCmn62 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.0892) / 13 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 7)) * lgM(6, 2) * cos(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 7)) * lgS(6, 2) * cos(2.0 * elSun));
	double dCmn63 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.0892) / 13 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 7)) * lgM(6, 3) * cos(3 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 7)) * lgS(6, 3) * cos(3.0 * elSun));
	double dCmn64 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.0892) / 13 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 7)) * lgM(6, 4) * cos(4 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 7)) * lgS(6, 4) * cos(4.0 * elSun));
	double dCmn65 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.0892) / 13 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 7)) * lgM(6, 5) * cos(5 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 7)) * lgS(6, 5) * cos(5.0 * elSun));
	double dCmn66 = 4 * PI * pow(RE_WGS84, 2) * 1025 / (5.9722e24) * (1 - 0.0892) / 13 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 7)) * lgM(6, 6) * cos(6 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 7)) * lgS(6, 6) * cos(6.0 * elSun));

	double dSmn21 = -0.3075 / 5 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(2, 1) * sin(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(2, 1) * sin(1 * elSun));
	double dSmn22 = -0.3075 / 5 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 3)) * lgM(2, 2) * sin(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 3)) * lgS(2, 2) * sin(2 * elSun));

	double dSmn31 = -0.195 / 7 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 1) * sin(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 1) * sin(1 * elSun));
	double dSmn32 = -0.195 / 7 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 2) * sin(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 2) * sin(2 * elSun));
	double dSmn33 = -0.195 / 7 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 4)) * lgM(3, 3) * sin(3 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 4)) * lgS(3, 3) * sin(3 * elSun));

	double dSmn41 = -0.132 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 5)) * lgM(4, 1) * sin(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 5)) * lgS(4, 1) * sin(1 * elSun));
	double dSmn42 = -0.132 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 5)) * lgM(4, 2) * sin(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 5)) * lgS(4, 2) * sin(2 * elSun));
	double dSmn43 = -0.132 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 5)) * lgM(4, 3) * sin(3 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 5)) * lgS(4, 3) * sin(3 * elSun));
	double dSmn44 = -0.132 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 5)) * lgM(4, 4) * sin(4 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 5)) * lgS(4, 4) * sin(4 * elSun));

	double dSmn51 = -0.1032 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 6)) * lgM(5, 1) * sin(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 6)) * lgS(5, 1) * sin(1 * elSun));
	double dSmn52 = -0.1032 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 6)) * lgM(5, 2) * sin(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 6)) * lgS(5, 2) * sin(2 * elSun));
	double dSmn53 = -0.1032 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 6)) * lgM(5, 3) * sin(3 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 6)) * lgS(5, 3) * sin(3 * elSun));
	double dSmn54 = -0.1032 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 6)) * lgM(5, 4) * sin(4 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 6)) * lgS(5, 4) * sin(4 * elSun));
	double dSmn55 = -0.1032 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 6)) * lgM(5, 5) * sin(5 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 6)) * lgS(5, 5) * sin(5 * elSun));

	double dSmn61 = -0.0892 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 7)) * lgM(6, 1) * sin(1 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 7)) * lgS(6, 1) * sin(1 * elSun));
	double dSmn62 = -0.0892 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 7)) * lgM(6, 2) * sin(2 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 7)) * lgS(6, 2) * sin(2 * elSun));
	double dSmn63 = -0.0892 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 7)) * lgM(6, 3) * sin(3 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 7)) * lgS(6, 3) * sin(3 * elSun));
	double dSmn64 = -0.0892 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 7)) * lgM(6, 4) * sin(4 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 7)) * lgS(6, 4) * sin(4 * elSun));
	double dSmn65 = -0.0892 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 7)) * lgM(6, 5) * sin(5 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 7)) * lgS(6, 5) * sin(5 * elSun));
	double dSmn66 = -0.0892 / 9 * ((GM_Moon / GM_Earth) * (pow((RE_WGS84 / rhoMoon), 7)) * lgM(6, 6) * sin(6 * azMoon) + (GM_Sun / GM_Earth) * (pow((RE_WGS84 / rhoSun), 7)) * lgS(6, 6) * sin(6 * elSun));

	egmCoef.cmn(2, 0) += dCmn20;
	egmCoef.cmn(2, 1) += dCmn21;
	egmCoef.cmn(2, 2) += dCmn22;
	egmCoef.smn(2, 1) += dSmn21;
	egmCoef.smn(2, 2) += dSmn22;

	egmCoef.cmn(3, 0) += dCmn30;
	egmCoef.cmn(3, 1) += dCmn31;
	egmCoef.cmn(3, 2) += dCmn32;
	egmCoef.cmn(3, 3) += dCmn33;
	egmCoef.smn(3, 1) += dSmn31;
	egmCoef.smn(3, 2) += dSmn32;
	egmCoef.smn(3, 3) += dSmn33;

	egmCoef.cmn(4, 0) += dCmn40;
	egmCoef.cmn(4, 1) += dCmn41;
	egmCoef.cmn(4, 2) += dCmn42;
	egmCoef.cmn(4, 3) += dCmn43;
	egmCoef.cmn(4, 4) += dCmn44;
	egmCoef.smn(4, 1) += dSmn41;
	egmCoef.smn(4, 2) += dSmn42;
	egmCoef.smn(4, 3) += dSmn43;
	egmCoef.smn(4, 4) += dSmn44;

	egmCoef.cmn(5, 0) += dCmn50;
	egmCoef.cmn(5, 1) += dCmn51;
	egmCoef.cmn(5, 2) += dCmn52;
	egmCoef.cmn(5, 3) += dCmn53;
	egmCoef.cmn(5, 4) += dCmn54;
	egmCoef.cmn(5, 5) += dCmn55;
	egmCoef.smn(5, 1) += dSmn51;
	egmCoef.smn(5, 2) += dSmn52;
	egmCoef.smn(5, 3) += dSmn53;
	egmCoef.smn(5, 4) += dSmn54;
	egmCoef.smn(5, 5) += dSmn55;

	egmCoef.cmn(6, 0) += dCmn60;
	egmCoef.cmn(6, 1) += dCmn61;
	egmCoef.cmn(6, 2) += dCmn62;
	egmCoef.cmn(6, 3) += dCmn63;
	egmCoef.cmn(6, 4) += dCmn64;
	egmCoef.cmn(6, 5) += dCmn65;
	egmCoef.cmn(6, 6) += dCmn66;
	egmCoef.smn(6, 1) += dSmn61;
	egmCoef.smn(6, 2) += dSmn62;
	egmCoef.smn(6, 3) += dSmn63;
	egmCoef.smn(6, 4) += dSmn64;
	egmCoef.smn(6, 5) += dSmn65;
	egmCoef.smn(6, 6) += dSmn66;
}

Vector3d GravityModel::centralBodyGravityAcc(
	const double mjdUTC,
	const double *erpv,		 ///< xp, yp, ut1_utc, lod and leapsecond
	const Vector3d &rSat,	 ///< Satellite position vector in the inertial system
	const Matrix3d &mECI2BF, ///< Transformation matrix from ECI to central body fixed system
	bool bVarEq)			 ///< bVarEq = 1 if for the variational equation
{
	double q1 = 0, q2 = 0, q3 = 0;
	double b1, b2, b3;
	int nd;
	int mMax = 0;
	int nMax = 0;
	EGMCoef egmCoef = mEGMCoef;
	IERS iersIns;
	Vector3d vecRAEMoon;
	Vector3d vecRAESun;

	if (!bVarEq)
	{
		double dUT1_UTC = erpv[2];
		double dUTC_TAI = -(19 + erpv[4]);
		double xp = erpv[0];
		double yp = erpv[1];
		double lod = erpv[3];
		iersIns.Set(dUT1_UTC, dUTC_TAI, xp, yp, lod);
		double jdTT = iersIns.TT_UTC(mjdUTC) / 86400.0 + mjdUTC + JD2MJD;

		Vector3d rSun = Vector3d::Zero();
		mJPLEphemeris->getJPLEphemeris(jdTT, eSun, eEarth, rSun);
		rSun = mECI2BF * rSun; // from inertial coordinate to earth centred fixed coordiante
		Vector3d rMoon = Vector3d::Zero();
		mJPLEphemeris->getJPLEphemeris(jdTT, eMoon, eEarth, rMoon);
		rMoon = mECI2BF * rMoon;

		vecRAESun = CalcPolarAngles(rSun); // calculating the range, azimuth and altitude
		vecRAEMoon = CalcPolarAngles(rMoon);

		mMax = mEarthGravAccDeg.mMax;
		nMax = mEarthGravAccDeg.nMax;
	}
	else
	{
		mEarthTidesOpt.flagSolidEarthTides = 0;
		mEarthTidesOpt.flagOceanTides = 0;
		mMax = mEarthGravSTMDeg.mMax;
		nMax = mEarthGravSTMDeg.nMax;
	}

	if (mEarthTidesOpt.flagSolidEarthTides)
	{
		solidEarthTidesCorrection(mjdUTC, egmCoef, iersIns, vecRAESun, vecRAEMoon);
	}

	if (mEarthTidesOpt.flagOceanTides)
	{
		oceanTidesCorrection(mjdUTC, egmCoef, vecRAESun, vecRAEMoon);
	}

	const Vector3d rSat_bf = mECI2BF * rSat; // Body-fixed position

	double rSat_latgc = asin(rSat_bf(2) / rSat_bf.norm()); // Geocentric latitude of satellite (n)
	double rSat_longc = atan2(rSat_bf(1), rSat_bf(0));	   // Geocentric longitude of satellite (n)

	const auto tmp_m = Legendre(mMax, nMax, rSat_latgc);		// Legendre matrix given order/degree

	const MatrixXd& pnm = tmp_m.first;
	const MatrixXd& dpnm = tmp_m.second; 
	double dUdr = 0;
	double dUdlatgc = 0;
	double dUdlongc = 0;

	const double tmp_a =  (GM_Earth / rSat_bf.norm());
	const double tmp_b = (RE_WGS84 / rSat_bf.norm());
	const double tmp_c = (-GM_Earth / rSat_bf.squaredNorm());

	std::vector<double> sin_table(mMax + 1);
	std::vector<double> cos_table(mMax + 1);
	for (int n = 0; n <= mMax; n++) {
		sin_table[n] = sin(n * rSat_longc);
		cos_table[n] = cos(n * rSat_longc);
	}
	for (int m = 0; m <= nMax; m++)
	{

		for (int n = 0; n <= m; n++)
		{
			const double tmp_cmn = egmCoef.cmn.coeff(m, n);
			const double tmp_smn = egmCoef.smn.coeff(m, n);
			const double tmp_sin = sin_table[n];
			const double tmp_cos = cos_table[n];
			const double tmp_pnm = pnm.coeff(m,n);
			q1 += tmp_pnm * (tmp_cmn * tmp_cos + tmp_smn * tmp_sin);
			q2 += dpnm.coeff(m, n) * (tmp_cmn * tmp_cos + tmp_smn * tmp_sin);
			q3 += n * tmp_pnm * (tmp_smn * tmp_cos - tmp_cmn * tmp_sin);
		}

		const double tmp_d = pow(tmp_b,m);
		const double b1 = tmp_c * tmp_d * (m + 1);
		const double b2 = tmp_a * tmp_d;
		const double b3 = tmp_a * tmp_d;
		dUdr += q1 * b1;
		dUdlatgc += q2 * b2;
		dUdlongc += q3 * b3;

		q3 = 0;
		q2 = 0;
		q1 = 0;
	}

	double r2xy = SQR(rSat_bf.x()) + SQR(rSat_bf.y()); // Body-fixed acceleration

	Vector3d acc_bf = Vector3d::Zero();
	acc_bf.x() = (1 / rSat_bf.norm() * dUdr - rSat_bf(2) / (rSat_bf.squaredNorm() * sqrt(r2xy)) * dUdlatgc) * rSat_bf(0) - (1 / r2xy * dUdlongc) * rSat_bf(1);
	acc_bf.y() = (1 / rSat_bf.norm() * dUdr - rSat_bf(2) / (rSat_bf.squaredNorm() * sqrt(r2xy)) * dUdlatgc) * rSat_bf(1) + (1 / r2xy * dUdlongc) * rSat_bf(0);
	acc_bf.z() = (1 / rSat_bf.norm() * dUdr * rSat_bf(2)) + sqrt(r2xy) / rSat_bf.squaredNorm() * dUdlatgc;

	return mECI2BF.transpose() * acc_bf; // Inertial acceleration
}

/* Relativistic Effects
 *
 */
Vector3d GravityModel::relativityEffectsAcc(
	const Vector3d &rSat, ///< Inertial position of satellite (m)
	const Vector3d &vSat) ///< Inertial velocity of satellite (m/s)
{

	double rSatNorm = rSat.norm();
	double vSatNorm = vSat.norm();

	return EARTH_G_CONST / (pow(CLIGHT, 2) * pow(rSatNorm, 3)) * ((4 * EARTH_G_CONST / rSatNorm - pow(vSatNorm, 2)) * rSat + 4 * rSat.dot(vSat) * vSat);
}
