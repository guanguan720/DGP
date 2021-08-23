#pragma once
#include <complex>
#include "MeshViewer/MeshDefinition.h"

class CrossField
{
public:
	CrossField(Mesh& m);
	~CrossField();
	std::vector<OpenMesh::Vec3d> GetCrossField(void);

private:
	void Init(void);
	Eigen::VectorXd Calc_uf4(void);

private:
	Mesh& mesh;
	OpenMesh::FPropHandleT<OpenMesh::Vec3d> centroid;
	OpenMesh::FPropHandleT<OpenMesh::Vec3d> local_x_axis;
	OpenMesh::FPropHandleT<OpenMesh::Vec3d> local_y_axis;
	OpenMesh::HPropHandleT<std::complex<double>> local_halfedge_4;
	std::vector<std::complex<double>> uf_4;	
	
};
