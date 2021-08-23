#pragma once
#include "MeshViewer\MeshDefinition.h"
#include <Eigen/Sparse>

class TutteEmbedding
{
public:
	TutteEmbedding(Mesh &m);
	~TutteEmbedding();
	void DoTutte(void);

private:
	void BindBoundToCircle(void);
	void InnerVEmbedding(void);
private:
	Mesh& mesh;
};