#pragma once
#include "MeshViewer/MeshDefinition.h"
#include "AABBTree/AABB_Tree.h"
#include "AABBTree/TinyVector.h"

class Remeshing
{
public:
	Remeshing(Mesh& m);
	~Remeshing();
	void DoRemeshing(void);

private:
	void Init(void);
	void SplitLong(void);
	void CollapseShort(void);
	void Flip(void);
	void Smoothing(void);
	void Project(void);

private:
	Mesh& mesh;
	double TargetLenth;
	double ThreasholdMax;
	double ThreasholdMin;
	std::vector<int> Degree;
	OpenMesh::EPropHandleT<double> EdgeLenth;
	AABB_Tree* ABtree;
};