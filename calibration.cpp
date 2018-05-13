#include "stdAfx.h"
#include "Assignment4.h"

// Assignment4 Source File 

////////////////////////////////////////////////////////////////////////////////
// A brief description of C2DPoint and C3DPoint
//
// class C2DPoint
// {
// public:
//		double x; // stores the x coordinate
//		double y; // stores the y coordinate
// };
//
// class C3DPoint
// {
// public:
//		double x; // stores the x coordinate
//		double y; // stores the y coordinate
//		double z; // stores the y coordinate
// };
//
// Note:
// Although C2DPoint and C3DPoint do contain other data members,
// don't bother with them
//

BOOL CCamera::Calibrate(const vector<C2DPoint*>& src2D, const vector<C3DPoint*>& src3D,
	const vector<C2DPoint*>& corners, CMatrix<double>& matPrj)
{
	// INPUT:
	//     vector<C2DPoint*>& src2D      This contains a list of 2D coordinates for the image points specified 
	//                                   by the user in the image. Each point in this list is in 1-to-1 
	//                                   correspondence with the point having the same index in src3D.
	//
	//     vector<C3DPoint*>& src3D      This contains a list of 3D coordinates for the image points specified
	//                                   by the user in the image. Each point in this list is in 1-to-1 
	//                                   correspondence with the point having the same index in src2D.
	//
	//     vector<C2DPoint*>& corners    This contains a list of 2D coordinates for the detected corners.
	//
	// OUTPUT:
	//     CMatrix<double>& pPrjMatrix   A 3x4 camera projection matrix computed from the detected corners.
	//
	// Please refer to the tutorial for the usage of the related libraries.


	//////////////////////////
	// Begin your code here

	// Step 1: Classify the input 3D points into points on the x-z planes, points on
	//         the y-z plane, and points not on any calibration plane
	//Initialize src2D_xz, src3D_xz, src2D_yz, src3D_yz
	vector<C2DPoint*> src2D_xz;
	vector<C3DPoint*> src3D_xz;
	vector<C2DPoint*> src2D_yz;
	vector<C3DPoint*> src3D_yz;
	vector<C2DPoint*> src2D_other;
	vector<C3DPoint*> src3D_other;
	//Classify
	for (int i = 0; i < src3D.size(); i++) {
		C3DPoint* pPnt = src3D[i];
		C2DPoint* pPnt2 = src2D[i];
		if (pPnt->y == 0) {
			//x-z plane
			src3D_xz.push_back(pPnt);
			src2D_xz.push_back(pPnt2);
		}
		else if (pPnt->x == 0) {
			//y-z plane
			src3D_yz.push_back(pPnt);
			src2D_yz.push_back(pPnt2);
		}
		else {
			src3D_other.push_back(pPnt);
			src2D_other.push_back(pPnt2);
		}
	}
	// Step 2: Estimate a plane-to-plane projectivity for each of the calibration planes
	//         using the input 2D/3D point pairs
	//1. For x-z plane
	//Initialize the D_xz matrix [X1 Y1 1 0 0 0 -u1X1 -u1Y1, 0 0 0 ...]
	CMatrix<double> D_xz(8, 8);
	int D_xz_rowNum = 0;
	for (int i = 0; i < src3D_xz.size(); i++) {
		C3DPoint* pPnt = src3D_xz[i];
		C2DPoint* pPnt2 = src2D_xz[i];
		D_xz(D_xz_rowNum, 0) = pPnt->x;
		D_xz(D_xz_rowNum, 1) = pPnt->z;
		D_xz(D_xz_rowNum, 2) = 1;
		D_xz(D_xz_rowNum, 3) = 0; D_xz(D_xz_rowNum, 4) = 0; D_xz(D_xz_rowNum, 5) = 0;
		D_xz(D_xz_rowNum, 6) = -1*(pPnt2->x * pPnt->x);
		D_xz(D_xz_rowNum, 7) = -1*(pPnt2->x * pPnt->z);
		D_xz_rowNum++;
		D_xz(D_xz_rowNum, 0) = 0; D_xz(D_xz_rowNum, 1) = 0; D_xz(D_xz_rowNum, 2) = 0;
		D_xz(D_xz_rowNum, 3) = pPnt->x;
		D_xz(D_xz_rowNum, 4) = pPnt->z;
		D_xz(D_xz_rowNum, 5) = 1;
		D_xz(D_xz_rowNum, 6) = -1*(pPnt2->y * pPnt->x);
		D_xz(D_xz_rowNum, 7) = -1*(pPnt2->y * pPnt->z);
		D_xz_rowNum++;
	}
	//Calculate Pseudoinverse
	CMatrix<double> D_xz_pi(8, 8);
	D_xz_pi = ((D_xz.Transpose() * D_xz).Inverse())*D_xz.Transpose();
	//Create f_xz
	CMatrix<double> f_xz(8, 1);
	int f_xz_rowNum = 0;
	for (int i = 0; i < 4; i++) {
		C2DPoint* pPnt2 = src2D_xz[i];
		f_xz(f_xz_rowNum, 0) = pPnt2->x;
		f_xz_rowNum++;
		f_xz(f_xz_rowNum, 0) = pPnt2->y;
		f_xz_rowNum++;
	}
	CMatrix<double> p_xz(8,1);
	p_xz = D_xz_pi * f_xz;
	//1. For y-z plane
	//Initialize the D_yz matrix [X1 Y1 1 0 0 0 -u1X1 -u1Y1, 0 0 0 ...]
	CMatrix<double> D_yz(8, 8);
	int D_yz_rowNum = 0;
	for (int i = 0; i < src3D_yz.size(); i++) {
		C3DPoint* pPnt = src3D_yz[i];
		C2DPoint* pPnt2 = src2D_yz[i];
		D_yz(D_yz_rowNum, 0) = pPnt->y;
		D_yz(D_yz_rowNum, 1) = pPnt->z;
		D_yz(D_yz_rowNum, 2) = 1;
		D_yz(D_yz_rowNum, 3) = 0; D_yz(D_yz_rowNum, 4) = 0; D_yz(D_yz_rowNum, 5) = 0;
		D_yz(D_yz_rowNum, 6) = -(pPnt2->x * pPnt->y);
		D_yz(D_yz_rowNum, 7) = -(pPnt2->x * pPnt->z);
		D_yz_rowNum++;
		D_yz(D_yz_rowNum, 0) = 0; D_yz(D_yz_rowNum, 1) = 0; D_yz(D_yz_rowNum, 2) = 0;
		D_yz(D_yz_rowNum, 3) = pPnt->y;
		D_yz(D_yz_rowNum, 4) = pPnt->z;
		D_yz(D_yz_rowNum, 5) = 1;
		D_yz(D_yz_rowNum, 6) = -(pPnt2->y * pPnt->y);
		D_yz(D_yz_rowNum, 7) = -(pPnt2->y * pPnt->z);
		D_yz_rowNum++;
	}
	//Calculate Pseudoinverse of D_yz
	CMatrix<double> D_yz_pi(8, 8);
	D_yz_pi = ((D_yz.Transpose() * D_yz).Inverse())*D_yz.Transpose();
	//Create f_yz
	CMatrix<double> f_yz(8, 1);
	int f_yz_rowNum = 0;
	for (int i = 0; i < 4; i++) {
		C2DPoint* pPnt2 = src2D_yz[i];
		f_yz(f_yz_rowNum, 0) = pPnt2->x;
		f_yz_rowNum++;
		f_yz(f_yz_rowNum, 0) = pPnt2->y;
		f_yz_rowNum++;
	}
	//Construct a CMatrix of p values in shape (8,1)
	CMatrix<double> p_yz(8,1);
	p_yz = D_yz_pi * f_yz;
	//Convert the p_xz and p_yz that are in the vector format to a 3*3 matrix format
	CMatrix<double> p_xz_33(3, 3);
	int pxz33Index = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			p_xz_33(i, j) = p_xz(pxz33Index, 0);
			pxz33Index++;
		}
	}
	p_xz_33(2, 2) = 1;
	CMatrix<double> p_yz_33(3, 3);
	int pyz33Index = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			p_yz_33(i, j) = p_yz(pyz33Index, 0);
			pyz33Index++;
		}
	}
	p_yz_33(2, 2) = 1;
	// Step 3: Using the estimated plane-to-plane projectivities, assign 3D coordinates
	//         to all the detected corners on the calibration pattern
	// Figure out the 3D points of corners on xz plane and store them on corners3D_xz
	vector<C2DPoint*> corners2D;
	vector<C3DPoint*> corners3D;
	for (double Z = 0.5; Z < 8; Z++) {
		for (double X = 0.5; X < 10; X++) {
			corners3D.push_back(new C3DPoint(X, 0, Z));
			CMatrix<double> pPnt3DMat_xz(3, 1);
			pPnt3DMat_xz(0, 0) = X;
			pPnt3DMat_xz(1, 0) = Z;
			pPnt3DMat_xz(2, 0) = 1;
			CMatrix<double> pPnt2DMat_xz(3, 1);
			pPnt2DMat_xz = (p_xz_33 * pPnt3DMat_xz);
			double u = pPnt2DMat_xz(0, 0) / pPnt2DMat_xz(2, 0);
			double v = pPnt2DMat_xz(1, 0) / pPnt2DMat_xz(2, 0);
			corners2D.push_back(new C2DPoint(u, v));
		}
	}
	// Figure out the 3D points of corners on yz plane and store them on corners3D_yz
	for (double Z = 0.5; Z < 8; Z++) {
		for (double Y = 0.5; Y < 10; Y++) {
			corners3D.push_back(new C3DPoint(0, Y, Z));
			CMatrix<double> pPnt3DMat_yz(3, 1);
			pPnt3DMat_yz(0, 0) = Y;
			pPnt3DMat_yz(1, 0) = Z;
			pPnt3DMat_yz(2, 0) = 1;
			CMatrix<double> pPnt2DMat_yz(3,1);
			pPnt2DMat_yz = (p_yz_33 * pPnt3DMat_yz);
			double u = pPnt2DMat_yz(0, 0) / pPnt2DMat_yz(2, 0);
			double v = pPnt2DMat_yz(1, 0) / pPnt2DMat_yz(2, 0);
			corners2D.push_back(new C2DPoint(u, v));
		}
	}
	
	//Use distance to find closest 2D corners
	vector<C2DPoint*> closest2D;
	vector<C3DPoint*> closest3D;
	for (int i = 0; i < corners2D.size(); i++) {
		double xCalc = corners2D[i]->x;
		double yCalc = corners2D[i]->y;
		double minD = 3;
		int minDIndex = -1;
		for (int j = 0; j < corners.size(); j++) {
			double xGiven = corners[j]->x;
			double yGiven = corners[j]->y;
			double distance = sqrt((xCalc - xGiven)*(xCalc - xGiven) + (yCalc - yGiven)*(yCalc - yGiven));
			if (distance < minD) {
				minD = distance;
				minDIndex = j;
			}
		}
		if (minDIndex != -1) {
			int ind = minDIndex;
			closest2D.push_back(new C2DPoint(corners[minDIndex]->x, corners[minDIndex]->y));
			closest3D.push_back(new C3DPoint(corners3D[i]->x, corners3D[i]->y, corners3D[i]->z));
		}
		

	}
	// Step 4: Estimate a 3x4 camera projection matrix from all the detected corners on
	//         the calibration pattern using linear least squares
	
	CMatrix<double> A(2 * closest2D.size(), 12), U, D, V;
	int A_rowNum = 0;
	for (int i = 0; i < closest2D.size(); i++) {
		C3DPoint* pPnt = closest3D[i];
		C2DPoint* pPnt2 = closest2D[i];
		A(A_rowNum, 0) = pPnt->x;
		A(A_rowNum, 1) = pPnt->y;
		A(A_rowNum, 2) = pPnt->z;
		A(A_rowNum, 3) = 1;
		A(A_rowNum, 4) = 0;  A(A_rowNum, 5) = 0; A(A_rowNum, 6) = 0; A(A_rowNum, 7) = 0;
		A(A_rowNum, 8) = -(pPnt2->x * pPnt->x);
		A(A_rowNum, 9) = -(pPnt2->x * pPnt->y);
		A(A_rowNum, 10) = -(pPnt2->x * pPnt->z);
		A(A_rowNum, 11) = -(pPnt2->x);
		A_rowNum++;
		//second row NEEDS MODIFICATION
		A(A_rowNum, 0) = 0;  A(A_rowNum, 1) = 0; A(A_rowNum, 2) = 0; A(A_rowNum, 3) = 0;
		A(A_rowNum, 4) = pPnt->x;
		A(A_rowNum, 5) = pPnt->y;
		A(A_rowNum, 6) = pPnt->z;
		A(A_rowNum, 7) = 1;
		A(A_rowNum, 8) = -(pPnt2->y * pPnt->x);
		A(A_rowNum, 9) = -(pPnt2->y * pPnt->y);
		A(A_rowNum, 10) = -(pPnt2->y * pPnt->z);
		A(A_rowNum, 11) = -(pPnt2->y);
		A_rowNum++;
	}
	A.SVD2(U, D, V);
	//	virtual CMatrix<Type> GetRow(int r) const;
	CMatrix<double> sol(12, 1);
	sol = V.GetCol(11);
	/*
	for (int i = 0; i < 12; i++) {
		double help = sol(i, 0);
	}
	*/
	double lastVal = sol(11, 0);
	sol = sol / lastVal;
	matPrj.SetSize(3, 4);
	int matPrjCounter = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			matPrj(i, j) = sol(matPrjCounter, 0);
			matPrjCounter++;
		}
	}
	sol.CleanUp();

	V.CleanUp();
	D.CleanUp();
	U.CleanUp();
	A.CleanUp();
	p_xz.CleanUp();
	p_xz_33.CleanUp();
	p_yz_33.CleanUp();
	p_yz.CleanUp();
	D_xz.CleanUp();
	D_xz_pi.CleanUp();
	D_yz.CleanUp();
	D_yz_pi.CleanUp();
	f_xz.CleanUp();
	f_yz.CleanUp();
	





	return TRUE;
}

void CCamera::Decompose(const CMatrix<double>& prjMatrix, CMatrix<double>& prjK, CMatrix<double>& prjRt)
{
	// INPUT:
	//     CMatrix<double>& prjMatrix    This is a 3x4 camera projection matrix to be decomposed.
	//
	// OUTPUT:
	//     CMatrix<double>& prjK         This is the 3x3 camera calibration matrix K.
	//
	//     CMatrix<double>& prjRt        This is the 3x4 matrix composed of the rigid body motion of the camera.
	//
	// Please refer to the tutorial for the usage of the related libraries.

	prjK.SetSize(3, 3, 0);
	prjRt.SetSize(3, 4, 0);

	//////////////////////////
	// Begin your code here
	
	// Step 1: Decompose the 3x3 sub-matrix composed of the first 3 columns of
	//         prjMatrix into the product of K and R using QR decomposition
	CMatrix<double>Q, R, K, Rotation;
	CMatrix<double>Temp(3, 3);
	Temp = prjMatrix.SubMat(0, 2, 0, 2);
	Temp.Inverse().QR(Q, R);
	K = R.Inverse();
	Rotation = Q.Transpose();
	// Step 2: Normalize the 3x3 camera calibration matrix K
	double a = K(2, 2);
	K = K / a;
	if (K(0, 0) < 0){
		for (int i = 0; i < 3; i++){
			K(i, 0) = -K(i, 0);
			Rotation(0, i) = -Rotation(0, i);
		}
	}
	if (K(1, 1) < 0){
		for (int i = 0; i < 3; i++){
			K(i, 1) = -K(i, 1);
			Rotation(1, i) = -Rotation(1, i);
		}
	}
	// Step 3: Compute the translation vector T from the last column of prjMatrix
	CMatrix<double>T, t;
	T = (1 / a)* (K.Inverse()*prjMatrix.GetCol(3));
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 4; j++){
			if (j != 3)
				prjRt(i, j) = Rotation(i, j);
			else
				prjRt(i, j) = T(i, 0);
		}
	}
	prjK = K;
	return;
}

void CCamera::Triangulate(const vector<CMatrix<double>*>& prjMats, const vector<vector<C2DPoint*>*>& src2Ds,
	vector<C3DPoint*>& res3D)
{
	// INPUT:
	//     vector<CMatrix<double>*> prjMats 	A list of projection matrices
	//
	//     vector<vector<C2DPoint*>*> src2Ds	A list of image point lists, each image point list is in 1-to-1
	//                                          correspondence with the projection matrix having the same index in prjMats.
	//
	// OUTPUT:
	//     vector<C3DPoint*> res3D				A list of 3D coordinates for the triangulated points.
	//
	// Note:
	//    - src2Ds can be considered as a 2D array with each 'column' containing the image positions 
	//      for the same 3D point in different images. If any of the image does not contain the image for a particular
	//      point, the corresponding element in src2Ds will be a Null vector. For example, if there are two images, 
	//      and we know 8 pairs of corresponding points, then
	//      
	//			prjMats.size() = 2
	//			src2Ds.size() = 2
	//			src2Ds[k]->size() = 8           // k >= 0 and k < no. of images - 1
	//    
	//    - If for any reason the 3D coordinates corresponding to a 'column' in src2Ds cannot be computed,
	//      please push in a NULL as a place holder. That is, you have to make sure that each point in res3D
	//      must be in 1-to-1 correspondence with a column in src2Ds, i.e., 
	//      
	//			src2Ds[k]->size() == src3D.size()	// k >= 0 and k < no. of images - 1
	//
	// Please refer to the tutorial for the usage of related libraries.

	//////////////////////////
	// Begin your code here
	//Step1: find numOfPoints and numOfImages
	int numOfPoints = src2Ds[0]->size();
	int numOfImages = prjMats.size();

	//Step2: for each point, calculate matrix A and solve for corresponding 3D coordinates.
	for (int i = 0; i < numOfPoints; i++) {
		//iterate through each point
		//set up matrix A & b
		CMatrix<double> A(numOfImages*2, 4,0);
		int AIndex = 0;
		CMatrix<double> b(numOfImages * 2);
		int bIndex = 0;
		for (int j = 0; j < numOfImages; j++) {
			if (src2Ds[j]->at(i) != NULL) {
				CMatrix<double> *imgPrjMat = prjMats[j];
				double x = src2Ds[j]->at(i)->x;
				double y = src2Ds[j]->at(i)->y;
				A(AIndex, 0) = (*imgPrjMat)(0, 0) - x*(*imgPrjMat)(2, 0);
				A(AIndex, 1) = (*imgPrjMat)(0, 1) - x*(*imgPrjMat)(2, 1);
				A(AIndex, 2) = (*imgPrjMat)(0, 2) - x*(*imgPrjMat)(2, 2);
				AIndex++;
				A(AIndex, 0) = (*imgPrjMat)(1, 0) - y*(*imgPrjMat)(2, 0);
				A(AIndex, 1) = (*imgPrjMat)(1, 1) - y*(*imgPrjMat)(2, 1);
				A(AIndex, 2) = (*imgPrjMat)(1, 2) - y*(*imgPrjMat)(2, 2);
				AIndex++;

				b(bIndex, 0) = -((*imgPrjMat)(0, 3) - x*(*imgPrjMat)(2, 3));
				bIndex++;
				b(bIndex, 0) = -((*imgPrjMat)(1, 3) - y*(*imgPrjMat)(2, 3));
				bIndex++;
			}
			
		}
		if (AIndex >= 4) {
			//Calculate Pseudoinverse
			CMatrix<double> A_pi(3, 4);
			A_pi = ((A.Transpose()*A).Inverse())*A.Transpose();
			CMatrix<double> sol(3, 1);
			sol = A_pi*b;
			res3D.push_back(new C3DPoint(sol(0, 0), sol(1, 0), sol(2, 0)));
		}
		else {
			res3D.push_back(NULL);
		}
		
	}
	
	return;
}


