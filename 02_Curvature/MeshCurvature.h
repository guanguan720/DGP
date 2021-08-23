#pragma once
#include "MeshViewer\MeshDefinition.h"
#include <OpenMesh/Core/Utils/PropertyManager.hh>
class Curvature
{
public:
	Curvature(Mesh m);
	~Curvature();

public:
	OpenMesh::Vec3d getCircumcenter(int id);
	OpenMesh::Vec3d getMixedVoronoiCellPoint(int id);
	float getMeanCurvature(int id);
	float getAbsoluteMeanCurvature(int id);
	float getGaussianCurvature(int id);

private:
	OpenMesh::Vec3d Laplacian(int id);

	// Initiate property of the mesh
	void InitArea(void);
	void InitCircumcenter(void);
	void InitMixedVoronoiCellPoint(void);
	void InitAngleWeightNormal(void);

	float sin(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2);
	float cos(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2);
	float cot(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2);
private:
	Mesh mesh;

	OpenMesh::VPropHandleT<double> Area;
	OpenMesh::VPropHandleT<OpenMesh::Vec3d> AngleWeightNormal;
	OpenMesh::FPropHandleT<OpenMesh::Vec3d> Circumcenter;
	OpenMesh::FPropHandleT<OpenMesh::Vec3d> MixedVoronoiCellPoint;

};

