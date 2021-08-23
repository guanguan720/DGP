#pragma once
#include "MeshViewer/MeshDefinition.h"
#include <queue>

class LloydIteration
{
public:
	LloydIteration(Mesh& m, int k, std::vector<int> InitSeed);
	~LloydIteration();
	void GetApproximation(void);
	
public:
	void Init(void);
	void UpdatePartition(void);
	void UpdateProxyFit(void);
	void UpdateSeedTri(void);

public:
	OpenMesh::FPropHandleT<int> FaceCluster;  //0,1,2,...,k-1
	Mesh& mesh;

public:
	int kNum;
	int faceNum;
	double CurEnergy;
	std::vector<OpenMesh::Vec3d> CurProxyNormal;
	std::vector<int> CurSeedIdx;	
	OpenMesh::FPropHandleT<double> FaceArea;
};