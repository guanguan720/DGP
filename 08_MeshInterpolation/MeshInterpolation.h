#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "MeshViewer/MeshDefinition.h"

class MeshInterpolation
{
public:
	MeshInterpolation(Mesh& m0, Mesh& m1);
	~MeshInterpolation();
	void DoMorphing(double t);

private:
	void Init(void);
	void LocalStep(double t);
	void GlobalStep(void);

private:
	Mesh& mesh0;
	Mesh& mesh1;
	OpenMesh::FPropHandleT<Eigen::Matrix2d> Jacobi0To1;
	OpenMesh::FPropHandleT<Eigen::Matrix2d> DesiredJacobi;
	OpenMesh::FPropHandleT<Eigen::Matrix2d> Rotation;
	OpenMesh::FPropHandleT<Eigen::Matrix2d> Scale;

	Eigen::SparseMatrix<double> CoefMatA;
	Eigen::SparseMatrix<double> CoefMat;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
};