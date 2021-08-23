#pragma once
#include "MeshViewer/MeshDefinition.h"
#include <Eigen/Sparse>

class MeanValuePara
{
public:
	MeanValuePara(Mesh& m): mesh(m) {};
	~MeanValuePara();
	void DoPara(void);

private:
	void CalcBoundary(void);
	void CalcInteriorPos(void);
	double cos(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2);
	double sin(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2);
	double tan_half(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2); //tan(a/2)

private:
	Mesh& mesh;
	OpenMesh::SmartHalfedgeHandle boundary_start;
	int boundaryNum;
	Eigen::MatrixXd b;
};

inline double MeanValuePara::cos(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2)
{
	return ((Point1 - AnglePoint) | (Point2 - AnglePoint) / ((Point1 - AnglePoint).norm() * (Point2 - AnglePoint).norm()));
}

inline double MeanValuePara::sin(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2)
{
	return (((Point1 - AnglePoint) % (Point2 - AnglePoint)).norm() / ((Point1 - AnglePoint).norm() * (Point2 - AnglePoint).norm()));
}

inline double MeanValuePara::tan_half(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2)
{
	return(sin(AnglePoint, Point1, Point2) / (1 + cos(AnglePoint, Point1, Point2)));
}