#pragma once
#include"src/MeshViewer/MeshDefinition.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

class Deformation
{
public:
	Deformation(Mesh& m);
	~Deformation();
	void DoDeformation(int iterNum, std::vector<int> Fixed, OpenMesh::Vec3d Active);
	
private:
	void InitProp(void);
	void Calc_Cot(void);
	void Calc_CoefMat(void);
	void UpdateRotation(void);
	void UpadatePos(void);

	double cot(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2);

private:
	Mesh& mesh;
	OpenMesh::HPropHandleT<double> cot_value;
	OpenMesh::VPropHandleT<Eigen::Matrix3d> R;
	Eigen::SparseMatrix<double> CoefMat;
	//Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver; //the result is terrible!!!
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

	std::vector<int> FixedPoints;
	OpenMesh::Vec3d ActivePos;
	std::vector<OpenMesh::Vec3d> LastFramePos;
};

inline double Deformation::cot(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2)
{
	return (((Point1 - AnglePoint) | (Point2 - AnglePoint)) / ((Point1 - AnglePoint) % (Point2 - AnglePoint)).norm());
}