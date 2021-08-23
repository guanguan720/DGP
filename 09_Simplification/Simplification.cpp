#include "HW\Simplification.h"
#include <OpenMesh/Core/Mesh/TriConnectivity.hh>

Simplification::Simplification(Mesh& m):mesh(m)
{
	Init_FaceProp();
	Init_VertexProp();
	Init_EdgeProp();
	Init_EOptimalPos();
	Init_EError();
}

Simplification::~Simplification()
{
	mesh.remove_property(F_ErrorMat);
	mesh.remove_property(V_ErrorMat);
	mesh.remove_property(E_ErrorMat);
	mesh.remove_property(E_OptimalPos);
	mesh.remove_property(E_Error);
}

void Simplification::SimplifyMesh()
{
	int nNum = mesh.n_vertices();
	struct Edge
	{
		OpenMesh::SmartEdgeHandle eh;
		int state;
	};

	typedef std::pair<double, Edge> pair;
	std::map<OpenMesh::SmartEdgeHandle, int> state; //record latest edge number


	struct compare
	{
		bool operator()(pair a, pair b)
		{
			return a.first > b.first;
		}
	};
	std::priority_queue<pair, std::vector<pair>, compare> queue;

	//std::priority_queue<pair, std::vector<pair>, std::greater<pair>> queue;
	for (auto eh : mesh.edges())
	{
		//if (!mesh.status(eh).tagged())
		//{
		pair tmp(mesh.property(E_Error, eh), Edge{ eh, 0 });
		queue.push(tmp);
		state.insert(std::make_pair(eh, 0));
		//}
	}

	pair min;
	int cout = nNum;
	int target_num = 0.5 * nNum;
	while (cout > target_num)
	{	
		/*double test_err = queue.top().first;*/
		min = queue.top();
		auto e = min.second;
		/*if (mesh.status(e).deleted())
		{
			queue.pop();
			continue;
		}*/	
		auto he = e.eh.h0();
		if (mesh.is_collapse_ok(he))
		{
			if (!mesh.status(e.eh).tagged() && e.state == state.find(e.eh)->second) //collapse ok && edge state number == latest edge state number
			{
				auto v = he.to();
				Eigen::Vector3d OptimalPos = mesh.property(E_OptimalPos, e.eh);
				OpenMesh::Vec3d tmp(OptimalPos[0], OptimalPos[1], OptimalPos[2]);
				mesh.set_point(v, tmp);
				mesh.collapse(he);
				cout--;
				mesh.property(V_ErrorMat, v) = mesh.property(E_ErrorMat, e.eh);
				//update 2-ring edge error 
				Update_FaceProp(v);
				Update_VertexProp(v);
				for (auto vv : v.vertices())
				{
					for (auto vve : v.edges())
					{
						Update_EdgeProp(vve);
						pair tmp(mesh.property(E_Error, vve), Edge{ vve, ++state[vve] });
						queue.push(tmp);
					}
				}
				mesh.status(e.eh).set_tagged(true);
			}
			queue.pop();
		}
		else
		{
			mesh.status(e.eh).set_tagged(true);
			queue.pop();
		}
	}
	mesh.garbage_collection();
	}
	
void Simplification::Init_FaceProp(void)
{
	mesh.add_property(F_ErrorMat, "F_ErrorMat");
	for (auto fh : mesh.faces())
	{
		auto n = mesh.normal(fh);
		auto x = mesh.point(fh.halfedge().to());
		double d = n | x;
		Eigen::Vector4d n_bar(n[0], n[1], n[2], -d);
		mesh.property(F_ErrorMat, fh) = n_bar * n_bar.transpose();
	}
}

void Simplification::Init_VertexProp(void)
{
	mesh.add_property(V_ErrorMat, "V_ErrorMat");	
	for (auto vh : mesh.vertices())
	{
		mesh.property(V_ErrorMat, vh) = Eigen::Matrix4d::Zero();
		for (auto vhf : vh.faces())
		{
			mesh.property(V_ErrorMat, vh) += mesh.property(F_ErrorMat, vhf);
		}
	}
}

void Simplification::Init_EdgeProp(void)
{
	mesh.add_property(E_ErrorMat, "E_ErrorMat");
	for (auto eh : mesh.edges())
	{
		mesh.property(E_ErrorMat, eh) = mesh.property(V_ErrorMat, eh.v0()) + mesh.property(V_ErrorMat, eh.v1());
	}
}

void Simplification::Init_EOptimalPos(void)
{
	mesh.add_property(E_OptimalPos, "E_OptimalPos");
	for (auto eh : mesh.edges())
	{
		Set_EOptimalPos(eh);
	}
}

void Simplification::Init_EError(void)
{
	mesh.add_property(E_Error, "E_Error");
	for (auto eh : mesh.edges())
	{
		Eigen::Matrix4d Qe = mesh.property(E_ErrorMat, eh);
		Eigen::Vector4d v;
		v << mesh.property(E_OptimalPos, eh), 1.0;
		mesh.property(E_Error, eh) = v.transpose() * Qe * v;
	}
}

void Simplification::Update_FaceProp(OpenMesh::SmartVertexHandle v)
{
	for (auto vf : v.faces())
	{
		mesh.update_normal(vf);
		auto n = mesh.normal(vf);
		auto x = mesh.point(vf.halfedge().to());
		double d = n | x;
		Eigen::Vector4d n_bar(n[0], n[1], n[2], -d);
		mesh.property(F_ErrorMat, vf) = n_bar * n_bar.transpose();
	}
}

void Simplification::Update_VertexProp(OpenMesh::SmartVertexHandle v)
{
	for (auto vv : v.vertices())
	{
		mesh.property(V_ErrorMat, vv) = Eigen::Matrix4d::Zero();
		for (auto vvf : vv.faces())
		{
			mesh.property(V_ErrorMat, vv) += mesh.property(F_ErrorMat, vvf);
		}
	}	
}

void Simplification::Update_EdgeProp(OpenMesh::SmartEdgeHandle e)
{	
	mesh.property(E_ErrorMat, e) = mesh.property(V_ErrorMat, e.v0()) + mesh.property(V_ErrorMat, e.v1());
	Set_EOptimalPos(e);
	Update_Error(e);	
}

void Simplification::Update_Error(OpenMesh::SmartEdgeHandle e)
{
	Eigen::Matrix4d Qe = mesh.property(E_ErrorMat, e);
	Eigen::Vector4d v;
	v << mesh.property(E_OptimalPos, e), 1.0;
	double error = v.transpose() * Qe * v;
	mesh.property(E_Error, e) = error;
}

void Simplification::Set_EOptimalPos(OpenMesh::SmartEdgeHandle e)
{
	Eigen::Matrix4d Qe = mesh.property(E_ErrorMat, e);
	Eigen::Matrix3d A = Qe.block<3, 3>(0, 0);
	Eigen::Vector3d b = Qe.block<3, 1>(0, 3);
	bool invertible; double determinant; Eigen::Matrix3d A_inv;
	A.computeInverseAndDetWithCheck(A_inv, determinant, invertible);
	if (invertible) //if (A.determinant() > 1e-5)
	{
		mesh.property(E_OptimalPos, e) = A_inv * -b;
	}
	else
	{
		OpenMesh::Vec3d mid_pos = (mesh.point(e.v0()) + mesh.point(e.v1())) / 2;
		Eigen::Vector3d tmp(mid_pos[0], mid_pos[1], mid_pos[2]);
		mesh.property(E_OptimalPos, e) = tmp;
	}
}

