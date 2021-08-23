#include "MeshCurvature.h"
#include <Eigen\Dense>
#include <iostream>
//#include <OpenMesh/Core/Mesh/DefaultTriMesh.hh>


Curvature::Curvature(Mesh m)
{
	mesh = m;
	
	InitCircumcenter();
	InitMixedVoronoiCellPoint();
	InitArea();	
	InitAngleWeightNormal();
}

Curvature::~Curvature()
{
	if (mesh.get_property_handle(Area,"Area"))
	{
		mesh.remove_property(Area);
	}
	if (mesh.get_property_handle(Circumcenter, "Circumcenter"))
	{
		mesh.remove_property(Circumcenter);
	}
	if (mesh.get_property_handle(MixedVoronoiCellPoint, "MixedVoronoiCellPoint"))
	{
		mesh.remove_property(MixedVoronoiCellPoint);
	}
	if (mesh.get_property_handle(AngleWeightNormal, "AngleWeightNormal"))
	{
		mesh.remove_property(AngleWeightNormal);
	}
}

float Curvature::getMeanCurvature(int id)
{
	if (!mesh.get_property_handle(AngleWeightNormal, "AngleWeightNormal"))
	{
		std::cout << "Unable to retrieve AngleWeightNormal property in MeshCurvature.cpp" << std::endl;		
		return -1.f;
	}
	OpenMesh::Vec3d lap = Laplacian(id);
	auto vh = mesh.vertex_handle(id);
	OpenMesh::Vec3d n = mesh.property(AngleWeightNormal, vh);
	/*float n1 = (lap[0] / (-2 * n[0]));
	float n2 = (lap[1] / (-2 * n[1]));
	float n3 = (lap[2] / (-2 * n[2]));*/
	return (lap | n) / (-2);	
}

float  Curvature::getAbsoluteMeanCurvature(int id)
{
	/*float temp = Laplacian(id).norm();
	return 0.5 * temp;*/
	float mean_cur = getMeanCurvature(id);
	return abs(mean_cur);
}

float  Curvature::getGaussianCurvature(int id)
{
	if (!mesh.get_property_handle(Area, "Area"))
	{
		std::cout << "Unable to retrieve Area property in MeshCurvature.cpp" << std::endl;
		return -1.f;
	}
	double area = mesh.property(Area, mesh.vertex_handle(id));
	auto vh = mesh.vertex_handle(id);
	OpenMesh::Vec3d p0 = mesh.point(vh);
	float theta_sum = 0.0f;
	for (auto vhf : mesh.vf_range(vh))
	{		
		std::vector<OpenMesh::Vec3d> p;
		for (auto fv : mesh.fv_range(vhf))
		{
			if (fv != vh)
			{
				p.push_back(mesh.point(fv));
			}
		}
		theta_sum += acos(cos(p0, p[0], p[1]));		
	}
	return (2 * 3.14159265358979323846 - theta_sum) / area;
}

OpenMesh::Vec3d Curvature::getCircumcenter(int id)
{
	if (!mesh.get_property_handle(Circumcenter, "Circumcenter"))
	{
		std::cout << "Unable to retrieve Circumcenter property in MeshCurvature.cpp" << std::endl;
		OpenMesh::Vec3d vec(0.f, 0.f, 0.f);
		return vec;
	}
	return mesh.property(Circumcenter, mesh.face_handle(id));
}

OpenMesh::Vec3d Curvature::getMixedVoronoiCellPoint(int id)
{
	if (!mesh.get_property_handle(MixedVoronoiCellPoint, "MixedVoronoiCellPoint"))
	{
		std::cout << "Unable to retrieve MixedVoronoiCellPoint property in MeshCurvature.cpp" << std::endl;
		OpenMesh::Vec3d vec(0.f, 0.f, 0.f);
		return vec;
	}
	return mesh.property(MixedVoronoiCellPoint, mesh.face_handle(id));
}

void Curvature::InitMixedVoronoiCellPoint(void)
{
	if (!mesh.get_property_handle(Circumcenter, "Circumcenter"))
	{
		std::cout << "Unable to retrieve Circumcenter property in MeshCurvature.cpp" << std::endl;
		return;
	}
	mesh.add_property(MixedVoronoiCellPoint, "MixedVoronoiCellPoint");
	for (auto fh : mesh.faces())
	{
		auto c = mesh.property(Circumcenter, fh);
		std::vector<OpenMesh::Vec3d> p;
		for (auto fhv:fh.vertices())
		{
			p.push_back(mesh.point(fhv));
		}
		OpenMesh::Vec3d sign0 = (c - p[0]) % (p[1] - p[0]);
		OpenMesh::Vec3d sign1 = (c - p[1]) % (p[2] - p[1]);
		OpenMesh::Vec3d sign2 = (c - p[2]) % (p[0] - p[2]);
		if (sign0[0] * sign1[0] + sign0[1] * sign1[1] + sign0[2] * sign1[2] >= 0)
		{
			if (sign1[0] * sign2[0] + sign1[1] * sign2[1] + sign1[2] * sign2[2] >= 0) 
			{
				mesh.property(MixedVoronoiCellPoint, fh) = c;  //c contain in triangle
			}
			else
			{
				mesh.property(MixedVoronoiCellPoint, fh) = (p[0] + p[2]) / 2; //c is out of side p02
			}
		}
		else
		{
			if (sign1[0] * sign2[0] + sign1[1] * sign2[1] + sign1[2] * sign2[2] >= 0) 
			{
				mesh.property(MixedVoronoiCellPoint, fh) = (p[0] + p[1]) / 2; //c is out of side p01
			}
			else
			{
				mesh.property(MixedVoronoiCellPoint, fh) = (p[1] + p[2]) / 2; // c is out of side p12
			}
		}
	}

}

void Curvature::InitAngleWeightNormal(void)
{
	mesh.add_property(AngleWeightNormal, "AngleWeightNormal");
	for (auto vh : mesh.vertices())
	{
		OpenMesh::Vec3d p0 = mesh.point(vh);
		OpenMesh::Vec3d vh_AngleWeightNormal(0.0, 0.0, 0.0);
		if (vh.is_boundary())
		{
			// 这段代码也适用非边界点(感觉比较低效) 但是else中提供了另一种思路			
			for (auto vhf : mesh.vf_range(vh))
			{
				OpenMesh::Vec3d f_normal = mesh.normal(vhf);
				std::vector<OpenMesh::Vec3d> p;
				for (auto fv : mesh.fv_range(vhf))
				{
					if (fv != vh)
					{
						p.push_back(mesh.point(fv));
					}
				}
				float theta = acos(cos(p0, p[0], p[1]));
				vh_AngleWeightNormal += theta * f_normal;
			}
			mesh.property(AngleWeightNormal, vh) = vh_AngleWeightNormal.normalize();
		}
		else
		{			
			for (auto voh : mesh.voh_range(vh))
			{
				auto f = mesh.face_handle(voh);
				OpenMesh::Vec3d f_normal = mesh.normal(f);
				OpenMesh::Vec3d p1 = mesh.point(voh.to());
				OpenMesh::Vec3d p2 = mesh.point(voh.next().to());
				float theta = acos(cos(p0, p1, p2));
				vh_AngleWeightNormal += theta * f_normal;
			}
			mesh.property(AngleWeightNormal, vh) = vh_AngleWeightNormal.normalize();
		}		
	}
}

void Curvature::InitArea(void)
{
	if (!mesh.get_property_handle(MixedVoronoiCellPoint, "MixedVoronoiCellPoint"))
	{
		std::cout << "Unable to retrieve MixedVoronoiCellPoint property in MeshCurvature.cpp" << std::endl;
		return;
	}
	mesh.add_property(Area, "Area");
	for (auto v:mesh.vertices())
	{
		mesh.property(Area, v) = 0.0f;
	}
	for (auto f:mesh.faces())
	{
		std::vector<OpenMesh::Vec3d> p(3);
		std::vector<OpenMesh::SmartVertexHandle> p_handle(3);
		int cout = 0;
		for (auto vh:f.vertices())
		{
			p[cout] = mesh.point(vh);
			p_handle[cout] = vh;
			cout++;
		}
		auto p01 = (p[0] + p[1]) / 2;
		auto p12 = (p[1] + p[2]) / 2;
		auto p02 = (p[0] + p[2]) / 2;
		auto c = mesh.property(MixedVoronoiCellPoint, f);
		float area0 = 0.5 * (p[0] - p01).norm() * (p[0] - c).norm() * sin(p[0], p01, c) + 0.5 * (p[0] - p02).norm() * (p[0] - c).norm() * sin(p[0], p02, c);
		float area1 = 0.5 * (p[1] - p01).norm() * (p[1] - c).norm() * sin(p[1], p01, c) + 0.5 * (p[1] - p12).norm() * (p[1] - c).norm() * sin(p[1], p12, c);
		float area2 = 0.5 * (p[2] - p02).norm() * (p[2] - c).norm() * sin(p[2], p02, c) + 0.5 * (p[2] - p12).norm() * (p[2] - c).norm() * sin(p[2], p12, c);
		mesh.property(Area, p_handle[0]) += area0;
		mesh.property(Area, p_handle[1]) += area1;
		mesh.property(Area, p_handle[2]) += area2;
	}
}

void Curvature::InitCircumcenter(void)
{		
	mesh.add_property(Circumcenter, "Circumcenter");
	for (auto f:mesh.faces())
	{
		auto n = mesh.normal(f);
		auto heh = mesh.halfedge_handle(f);
		auto p0 = mesh.point(mesh.from_vertex_handle(heh));
		auto p1 = mesh.point(mesh.to_vertex_handle(heh));
		auto p2 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)));
		auto n1 = (p1 - p0) % n;
		auto n2 = (p2 - p1) % n;
		OpenMesh::Vec3d tmp((n2[0] * n1[1] * (p0[0] + p1[0]) - n2[1] * n1[0] * (p1[0] + p2[0]) - n1[0] * n2[0] * (p0[1] - p2[1])) / (2 * (n2[0] * n1[1] - n2[1] * n1[0])),
			(n1[0] * n2[1] * (p0[1] + p1[1]) - n2[0] * n1[1] * (p1[1] + p2[1]) - n1[1] * n2[1] * (p0[0] - p2[0])) / (2 * (n2[1] * n1[0] - n2[0] * n1[1])),
			(n2[2] * n1[1] * (p0[2] + p1[2]) - n2[1] * n1[2] * (p1[2] + p2[2]) - n1[2] * n2[2] * (p0[1] - p2[1])) / (2 * (n2[2] * n1[1] - n2[1] * n1[2])));
		mesh.property(Circumcenter, f) = tmp;
	}
}


OpenMesh::Vec3d Curvature::Laplacian(int id)
{

	if (!mesh.get_property_handle(Area, "Area"))
	{
		std::cout << "Unable to retrieve Area property in MeshCurvature.cpp" << std::endl;
		OpenMesh::Vec3d lap(-1.0, -1.0, -1.0);
		return lap;
	}
	double area = mesh.property(Area, mesh.vertex_handle(id));
	auto xi_handle = mesh.vertex_handle(id);
	OpenMesh::Vec3d xi = mesh.point(xi_handle);
	OpenMesh::Vec3d lap(0.0f, 0.0f, 0.0f);
	for (auto voh:mesh.voh_range(xi_handle))
	{
		auto xj_handle = voh.to();
		auto alpha_handle = voh.next().to();
		auto beta_handle = voh.next().opp().next().to();
		auto xj = mesh.point(xj_handle);
		auto alpha = mesh.point(alpha_handle);
		auto beta = mesh.point(beta_handle);
		double cot_alpha = cot(alpha, xj, xi);
		double cot_beta = cot(beta, xj, xi);
		lap += (cot_alpha + cot_beta) * (xj - xi);
	}
	lap /= 1 / (2 * area);
	return lap;
}

float Curvature::sin(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2)
{
	float c = cos(AnglePoint, Point1, Point2);
	return (sqrt(1 - c * c));
}

float Curvature::cos(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2)
{
	float oppoE = (Point1 - Point2).norm();
	float adjE1 = (Point1 - AnglePoint).norm();
	float adjE2 = (Point2 - AnglePoint).norm();
	return ((adjE1 * adjE1 + adjE2 * adjE2 - oppoE * oppoE) / (2 * adjE1 * adjE2));
}

float Curvature::cot(OpenMesh::Vec3d AnglePoint, OpenMesh::Vec3d Point1, OpenMesh::Vec3d Point2)
{
	float c = cos(AnglePoint, Point1, Point2);;
	float s = sqrt(1 - c * c);
	return (c / s);
}