#include "HW/DelaunayTriangulation.h"

DelaunayTriangulation::DelaunayTriangulation(Mesh& m) :mesh(m)
{
	mesh.add_property(Circumcenter);
	mesh.add_property(FaceArea);
	mesh.add_property(OneRingArea);
}

DelaunayTriangulation::~DelaunayTriangulation()
{
	mesh.remove_property(Circumcenter);
	mesh.remove_property(FaceArea);
	mesh.remove_property(OneRingArea);
}

void DelaunayTriangulation::DoDelaunay(void)
{
	for (size_t i = 0; i < 3000; i++)
	{		
		UpdateVertexPos();
		UpdateTriangulation();
	}
}

void DelaunayTriangulation::UpdateTriangulation(void)
{
	OpenMesh::Vec3d p0, p1, p2, p3;
	double cos1, cos2, sin1, sin2;
	for (auto& eh : mesh.edges())
	{
		if (!mesh.is_boundary(eh))
		{
			auto h0 = eh.h0();
			auto h1 = eh.h1();
			auto v0 = h0.to();
			auto v1 = h1.to();
			auto v2 = h0.next().to();
			auto v3 = h1.next().to();
			p0 = mesh.point(v0);
			p1 = mesh.point(v1);
			p2 = mesh.point(v2);
			p3 = mesh.point(v3);
			cos1 = cos(p2, p0, p1);
			cos2 = cos(p3, p0, p1);
			sin1 = sin(p2, p0, p1);
			sin2 = sin(p3, p0, p1);
			if ((sin1 * cos2 + cos1 * sin2) < 0 && mesh.is_flip_ok(eh))
			{
				mesh.flip(eh);
			}
		}
	}
}

void DelaunayTriangulation::UpdateVertexPos(void)
{
	UpdateCircumcenter();
	UpdateFaceArea();
	Update1ringArea();
	OpenMesh::Vec3d new_pos;
	for (auto& vh : mesh.vertices())
	{		
		if (!mesh.is_boundary(vh))
		{
			new_pos = OpenMesh::Vec3d(0.0, 0.0, 0.0);
			for (auto& vf : vh.faces())
			{
				new_pos += mesh.property(FaceArea, vf) * mesh.property(Circumcenter, vf);
			}
			new_pos /= mesh.property(OneRingArea, vh);
			mesh.set_point(vh, new_pos);
		}
	}
}

void DelaunayTriangulation::UpdateCircumcenter(void)
{
	mesh.update_normals();
	OpenMesh::Vec3d tmp;
	OpenMesh::Vec3d v0, v1, v2, cross01, cross12;
	for (auto f : mesh.faces())
	{	
		auto heh = mesh.halfedge_handle(f);
		auto p0 = mesh.point(mesh.from_vertex_handle(heh));
		auto p1 = mesh.point(mesh.to_vertex_handle(heh));
		auto p2 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)));	
		if (std::abs(p0[2]) < 1e-5 && std::abs(p1[2]) < 1e-5 && std::abs(p2[2]) < 1e-5) //triangle in xy-plane
		{
			auto n1 = OpenMesh::Vec2d(p0[1] - p1[1], p1[0] - p0[0]);
			auto n2 = OpenMesh::Vec2d(p1[1] - p2[1], p2[0] - p1[0]);
			tmp[0] = (n2[0] * n1[1] * (p0[0] + p1[0]) - n2[1] * n1[0] * (p1[0] + p2[0]) - n1[0] * n2[0] * (p0[1] - p2[1])) / (2 * (n2[0] * n1[1] - n2[1] * n1[0]));
			tmp[1] = (n1[0] * n2[1] * (p0[1] + p1[1]) - n2[0] * n1[1] * (p1[1] + p2[1]) - n1[1] * n2[1] * (p0[0] - p2[0])) / (2 * (n2[1] * n1[0] - n2[0] * n1[1]));
			tmp[2] = 0.0;
				
		}
		else if (std::abs(p0[1]) < 1e-5 && std::abs(p1[1]) < 1e-5 && std::abs(p2[1]) < 1e-5) //triangle in xz-plane
		{
			auto n1 = OpenMesh::Vec2d(p0[2] - p1[2], p1[0] - p0[0]);
			auto n2 = OpenMesh::Vec2d(p1[2] - p2[2], p2[0] - p1[0]);
			tmp[0] = (n2[0] * n1[1] * (p0[0] + p1[0]) - n2[1] * n1[0] * (p1[0] + p2[0]) - n1[0] * n2[0] * (p0[2] - p2[2])) / (2 * (n2[0] * n1[1] - n2[1] * n1[0]));
			tmp[1] = 0.0;
			tmp[2] = (n2[1] * n1[0] * (p0[2] + p1[2]) - n2[0] * n1[1] * (p1[2] + p2[2]) - n1[1] * n2[1] * (p0[0] - p2[0])) / (2 * (n2[1] * n1[0] - n2[0] * n1[1]));
		}
		else if (std::abs(p0[0]) < 1e-5 && std::abs(p1[0]) < 1e-5 && std::abs(p2[0]) < 1e-5) //triangle in yz-plane
		{
			auto n1 = OpenMesh::Vec2d(p0[2] - p1[2], p1[1] - p0[1]);
			auto n2 = OpenMesh::Vec2d(p1[2] - p2[2], p2[1] - p1[1]);
			tmp[0] = 0.0;
			tmp[1] = (n2[0] * n1[1] * (p0[1] + p1[1]) - n2[1] * n1[0] * (p1[1] + p2[1]) - n1[0] * n2[0] * (p0[2] - p2[2])) / (2 * (n2[0] * n1[1] - n2[1] * n1[0]));
			tmp[2] = (n1[0] * n2[1] * (p0[2] + p1[2]) - n2[0] * n1[1] * (p1[2] + p2[2]) - n1[1] * n2[1] * (p0[2] - p2[2])) / (2 * (n2[1] * n1[0] - n2[0] * n1[1]));
		}
		else
		{
			auto n = mesh.normal(f);
			auto n1 = (p1 - p0) % n;
			auto n2 = (p2 - p1) % n;
			tmp[0] = (n2[0] * n1[1] * (p0[0] + p1[0]) - n2[1] * n1[0] * (p1[0] + p2[0]) - n1[0] * n2[0] * (p0[1] - p2[1])) / (2 * (n2[0] * n1[1] - n2[1] * n1[0]));
			tmp[1] = (n1[0] * n2[1] * (p0[1] + p1[1]) - n2[0] * n1[1] * (p1[1] + p2[1]) - n1[1] * n2[1] * (p0[0] - p2[0])) / (2 * (n2[1] * n1[0] - n2[0] * n1[1]));
			tmp[2] = (n2[2] * n1[1] * (p0[2] + p1[2]) - n2[1] * n1[2] * (p1[2] + p2[2]) - n1[2] * n2[2] * (p0[1] - p2[1])) / (2 * (n2[2] * n1[1] - n2[1] * n1[2]));
		}
		if (!mesh.is_boundary(f))
		{
			mesh.property(Circumcenter, f) = tmp;
		}
		else //judge whether the circumcenter is outside the mesh. If so, move the circumcenter at the projection of the boundary edge
		{		
			for (auto& fe : f.edges())
			{
				if (mesh.is_boundary(fe))
				{
					v0 = tmp - mesh.point(fe.v0());
					v1 = tmp - mesh.point(fe.v1());
					v2 = tmp - mesh.point(fe.h0().next().to());
					cross01 = v0 % v1;
					cross12 = v1 % v2;
					if (cross01[0] * cross12[0] <= 0) //Circumcenter is outside the mesh!
					{
						tmp = mesh.point(fe.v0()) + ((tmp - mesh.point(fe.v0())) | ((mesh.point(fe.v1()) - mesh.point(fe.v0())).normalize())) * (mesh.point(fe.v1()) - mesh.point(fe.v0())).normalize(); //projection
						break;
					}
				}
			}
			mesh.property(Circumcenter, f) = tmp;
		}
	}
}

void DelaunayTriangulation::UpdateFaceArea(void)
{
	for (auto& fh : mesh.faces())
	{
		mesh.property(FaceArea, fh) = mesh.calc_face_area(fh);
	}
}

void DelaunayTriangulation::Update1ringArea(void)
{
	double tmp;
	for (auto& vh : mesh.vertices())
	{
		if (!mesh.is_boundary(vh))
		{
			tmp = 0.0;
			for (auto& vf : vh.faces())
			{
				tmp += mesh.property(FaceArea, vf);
			}
			mesh.property(OneRingArea, vh) = tmp;
		}	
	}
}