#pragma once

#include "MeshViewer\MeshDefinition.h"
#include <Eigen/Dense>
//#include <Eigen/LU> 
#include <queue>
#include <vector>
//#include <algorithm>

class Simplification
{
public:
	Simplification(Mesh& m);
	~Simplification();

	void SimplifyMesh();

private:
	void Init_FaceProp(void);
	void Init_VertexProp(void);
	void Init_EdgeProp(void);
	void Init_EOptimalPos(void);
	void Init_EError(void);

	void Update_FaceProp(OpenMesh::SmartVertexHandle v);	//update v's adjacent faces property
	void Update_VertexProp(OpenMesh::SmartVertexHandle v);	//update v's adjacent vertices property
	void Update_EdgeProp(OpenMesh::SmartEdgeHandle e);		//update e's edge property
	void Update_Error(OpenMesh::SmartEdgeHandle e);			//update e's edge property

	void Set_EOptimalPos(OpenMesh::SmartEdgeHandle e);


private:
	Mesh& mesh;
	OpenMesh::FPropHandleT<Eigen::Matrix4d> F_ErrorMat;

	OpenMesh::VPropHandleT<Eigen::Matrix4d> V_ErrorMat;

	OpenMesh::EPropHandleT<Eigen::Matrix4d> E_ErrorMat;
	OpenMesh::EPropHandleT<Eigen::Vector3d> E_OptimalPos;
	OpenMesh::EPropHandleT<double> E_Error;
	
};
