#pragma once
#include "MeshViewer\MeshDefinition.h"
#include "HW\TutteEmbedding.h"
#include <Eigen\Dense>
#include <Eigen\Sparse>

class ARAP
{
public:
	ARAP(Mesh& m);
	~ARAP();
	void DoARAP(void);
private:
	void InitPara(void);
	void InitFaceProp(void);
	void CalcCot(void);
	double cos(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2);
	double sin(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2);
	double cot(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2);
	void LocalStep(void);
	void GlobalStep(void);
	void InitGlobalMat(void);
	void GenRighthandTerm(void);
	double getEnergy(void);
private:
	Mesh& mesh;
	OpenMesh::FPropHandleT<Eigen::Matrix2d> jacobi;
	OpenMesh::FPropHandleT<double> face_area;
	OpenMesh::FPropHandleT<Eigen::Matrix2d> rotationMat;
	OpenMesh::FPropHandleT<Eigen::Matrix<double, 3, 2>> ori_face_uv;  // put original mesh face onto a plane
	OpenMesh::FPropHandleT<Eigen::Matrix<double, 3, 2>> gradientMat;

	Eigen::SparseMatrix<double> globalMat;
	Eigen::MatrixXd b;
	//Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	Eigen::Vector2d b0;
	std::vector<double> cot_value;
};

inline double ARAP::cos(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2)
{
	float oppoE = (Point1 - Point2).norm();
	float adjE1 = (Point1 - AnglePoint).norm();
	float adjE2 = (Point2 - AnglePoint).norm();
	return ((adjE1 * adjE1 + adjE2 * adjE2 - oppoE * oppoE) / (2 * adjE1 * adjE2));
}

inline double ARAP::sin(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2)
{
	float c = cos(AnglePoint, Point1, Point2);
	return (sqrt(1 - c * c));
}


inline double ARAP::cot(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2)
{
	return (((Point1 - AnglePoint) | (Point2 - AnglePoint)) / ((Point1 - AnglePoint) % (Point2 - AnglePoint)).norm());
}

