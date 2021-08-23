#pragma once
#include "MeshViewer/MeshDefinition.h"

class DelaunayTriangulation
{
public:
	DelaunayTriangulation(Mesh& m);
	~DelaunayTriangulation();
	void DoDelaunay(void);

private:
	void UpdateTriangulation(void);
	void UpdateVertexPos(void);
	void UpdateCircumcenter(void);
	void UpdateFaceArea(void);
	void Update1ringArea(void);
	double cos(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2);
	double sin(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2);

private:
	Mesh& mesh;
	OpenMesh::FPropHandleT<OpenMesh::Vec3d> Circumcenter;
	OpenMesh::FPropHandleT<double> FaceArea;
	OpenMesh::VPropHandleT <double> OneRingArea;
};

inline double DelaunayTriangulation::cos(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2)
{
	return ((Point1 - AnglePoint) | (Point2 - AnglePoint)) / ((Point1 - AnglePoint).norm() * (Point2 - AnglePoint).norm());
}

inline double DelaunayTriangulation::sin(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2)
{
	return ((Point1 - AnglePoint) % (Point2 - AnglePoint)).norm() / ((Point1 - AnglePoint).norm() * (Point2 - AnglePoint).norm());
}