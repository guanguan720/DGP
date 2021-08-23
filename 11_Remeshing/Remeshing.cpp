#include "HW/Remeshing.h"


Remeshing::Remeshing(Mesh& m):mesh(m)
{
	Init();
}

Remeshing::~Remeshing()
{
	mesh.remove_property(EdgeLenth);
	delete ABtree;
}

void Remeshing::DoRemeshing(void)
{
	for (size_t i = 0; i < 5; i++)
	{
		SplitLong();
		CollapseShort();
		Flip();
		Smoothing();
		Project();
	}
	mesh.update_normals();
}

void Remeshing::Init(void)
{
	//Init target lenth & edge lenth
	mesh.add_property(EdgeLenth);
	TargetLenth = 0.0;
	double tmp;
	for (auto& eh : mesh.edges())
	{
		tmp = mesh.calc_edge_length(eh);
		TargetLenth += tmp;
		mesh.property(EdgeLenth, eh) = tmp;
	}
	TargetLenth /= mesh.n_edges();
	TargetLenth *= 0.8;
	ThreasholdMax = 4 * TargetLenth / 3;
	ThreasholdMin = 4 * TargetLenth / 5;
	//Init AABB tree
	std::vector<Vector3f> List;
	Vector3f p;
	OpenMesh::Vec3d pos;
	for (auto& fh : mesh.faces())
	{
		for (auto& fv : fh.vertices())
		{
			pos = mesh.point(fv);
			p[0] = pos[0];
			p[1] = pos[1];
			p[2] = pos[2];
			List.push_back(p);
		}
	}
	ABtree = new AABB_Tree(List);
	
}

void Remeshing::SplitLong(void)
{
	OpenMesh::Vec3d p1, p2;
	for (auto& eh : mesh.edges())
	{
		if (mesh.property(EdgeLenth, eh) > ThreasholdMax)
		{
			p1 = mesh.point(eh.h0().to());
			p2 = mesh.point(eh.h0().from());
			auto new_vh = mesh.add_vertex((p1 + p2) / 2);
			mesh.split(eh, new_vh);
			for (auto& ve : new_vh.edges())
			{
				mesh.property(EdgeLenth, ve) = mesh.calc_edge_length(ve);
			}
		}
	}
}

void Remeshing::CollapseShort(void)
{
	OpenMesh::Vec3d p0, p1, mid;
	double len;
	bool IsLongEdge;
	for (auto& eh : mesh.edges())
	{
		if (mesh.property(EdgeLenth, eh) < ThreasholdMin && !mesh.is_boundary(eh))
		{
			IsLongEdge = false;
			auto he = eh.h0();
			auto v0 = he.from();
			auto v1 = he.to();
			p0 = mesh.point(v0);
			p1 = mesh.point(v1);
			if (!mesh.is_boundary(v0) && !mesh.is_boundary(v1))
			{
				mid = (p0 + p1) / 2;
			}
			if(mesh.is_boundary(v1))
			{
				mid = p1;				//avoid moving boundary vertex when collapsing
			}	
			if (mesh.is_boundary(v0))
			{
				mid = p0;
			}
			for (auto& v0v : v0.vertices())
			{
				len = (mid - mesh.point(v0v)).norm();
				if (len > ThreasholdMax)
				{
					IsLongEdge = true;
					break;
				}
			}
			if (!IsLongEdge)
			{
				for (auto& v1v : v1.vertices())
				{
					len = (mid - mesh.point(v1v)).norm();
					if (len > ThreasholdMax)
					{
						IsLongEdge = true;
						break;
					}
				}
			}
			if (!IsLongEdge && mesh.is_collapse_ok(he))
			{
				mesh.set_point(v1, mid);			
				mesh.collapse(he);
				//update edge length
				for (auto& ve : v1.edges())
				{
					mesh.property(EdgeLenth, ve) = mesh.calc_edge_length(ve);
				}
			}
		}
	}
	mesh.garbage_collection();
}

void Remeshing::Flip(void)
{
	Degree.resize(mesh.n_vertices());
	std::fill(Degree.begin(), Degree.end(), 0);
	for (auto& vh : mesh.vertices())
	{
		for (auto& vv : vh.vertices())
		{
			Degree[vh.idx()] += 1;
		}
	}
	double preEnergy, Energy;
	std::vector<OpenMesh::SmartVertexHandle> v(4);
	for (auto& eh : mesh.edges())
	{
		if (!mesh.is_boundary(eh) && mesh.is_flip_ok(eh))
		{
			auto h0 = eh.h0();
			auto h1 = eh.h1();
			v[0] = h0.to();
			v[1] = h1.to();
			v[2] = h0.next().to();
			v[3] = h1.next().to();
			preEnergy = Energy = 0.0;
			for (size_t i = 0; i < 4; i++)
			{
				if (mesh.is_boundary(v[i]))
				{
					preEnergy += std::abs(Degree[v[i].idx()] - 4);
				}
				else
				{
					preEnergy += std::abs(Degree[v[i].idx()] - 6);
				}
			}
			for (size_t i = 0; i < 2; i++)
			{
				if (mesh.is_boundary(v[i]))
				{
					Energy += std::abs(Degree[v[i].idx()] - 1 - 4);
				}
				else
				{
					Energy += std::abs(Degree[v[i].idx()] - 1 - 6);
				}
			}
			for (size_t i = 2; i < 4; i++)
			{
				if (mesh.is_boundary(v[i]))
				{
					Energy += std::abs(Degree[v[i].idx()] + 1 - 4);
				}
				else
				{
					Energy += std::abs(Degree[v[i].idx()] + 1 - 6);
				}
			}
			if (Energy < preEnergy)
			{
				mesh.flip(eh);
				Degree[v[0].idx()] -= 1;
				Degree[v[1].idx()] -= 1;
				Degree[v[2].idx()] += 1;
				Degree[v[3].idx()] += 1;
				mesh.property(EdgeLenth, eh) = mesh.calc_edge_length(eh);
			}
		}
	}
}

void Remeshing::Smoothing(void)
{
	mesh.update_normals();
	std::vector<OpenMesh::Vec3d> smoothPos(mesh.n_vertices(), OpenMesh::Vec3d(0.0, 0.0, 0.0));
	int deg;
	for (auto& vh : mesh.vertices())
	{
		deg = 0;
		for (auto& vv : vh.vertices())
		{
			smoothPos[vh.idx()] += mesh.point(vv);
			deg++;
		}
		smoothPos[vh.idx()] /= deg;
	}
	OpenMesh::Vec3d desired;
	for (auto& vh : mesh.vertices())
	{
		if (!mesh.is_boundary(vh))
		{
			desired = smoothPos[vh.idx()] + mesh.normal(vh) * (mesh.normal(vh) | (mesh.point(vh) - smoothPos[vh.idx()]));
			mesh.set_point(vh, desired);
		}	
	}
}

void Remeshing::Project(void)
{
	OpenMesh::Vec3d pos, new_pos;
	Vector3f p, nearestP;
	for (auto& vh : mesh.vertices())
	{
		if (!mesh.is_boundary(vh))
		{
			pos = mesh.point(vh);
			p[0] = pos[0];
			p[1] = pos[1];
			p[2] = pos[2];
			ABtree->findNearstPoint(p, nearestP);
			new_pos[0] = double(nearestP[0]);
			new_pos[1] = double(nearestP[1]);
			new_pos[2] = double(nearestP[2]);
			mesh.set_point(vh, new_pos);
		}		
	}
}