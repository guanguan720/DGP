#pragma once
#include"MeshViewer/MeshDefinition.h"

class MeshSmoothing
{
public:
	MeshSmoothing(Mesh& m);
	~MeshSmoothing();
	void SmoothMesh(void);

private:
	void InitFaceProp(void);
	void UpdateNormal(void);
	void UpdateVertex(void);

	double CalculateVolume(void);
	void UpdateVolume(void);
	
private:
	Mesh& mesh;
	OpenMesh::FPropHandleT<double> FArea;
	OpenMesh::FPropHandleT<OpenMesh::Vec3d> FCentroid;

	int K1;
	int K2;
	double original_volume;
};