#pragma once

// header inclusion
#include "math.h"
#include "Matrix.h"
#include <vector>
using namespace std;

#include "2DPointMap.h"
#include "3DPointMap.h"

// some constants
#define CALIBRATION_FUZZINESS 16.0
#define ALIGN_FUZZINESS	64.0

///////////////////////////////////
// CCamera class's declaration
///////////////////////////////////

///////////////////////////////////
// IMPORTANT!
//
// Don't modify the declaration of
// CCamera!
//
///////////////////////////////////
class __declspec(dllexport) CCamera
{
public:
	///////////////////////////////////////////
	// For calibration from Calibration Grid //
	///////////////////////////////////////////

	virtual BOOL Calibrate(const vector<C2DPoint*>& src2D,
		                     const vector<C3DPoint*>& src3D,
							 const vector<C2DPoint*>& corners,
							 CMatrix<double>& matPrj);


	// decompose a projection matrix into K and Rt
	virtual void Decompose(const CMatrix<double>& prjMatrix,
										   CMatrix<double>& prjK,
										   CMatrix<double>& prjRt);

	// For Reconstruction
	virtual void Triangulate(const vector<CMatrix<double>*>& prjMats, 
							 const vector<vector<C2DPoint*>*>& src2Ds,
							 vector<C3DPoint*>& res3D);
};
