#include "HW/MeshSmoothing.h"

MeshSmoothing::MeshSmoothing(Mesh& m) :mesh(m)
{
	InitFaceProp();
	K1 = 15;
	K2 = 20;
	original_volume = CalculateVolume();
}

MeshSmoothing::~MeshSmoothing()
{
	mesh.remove_property(FArea);
	mesh.remove_property(FCentroid);
}

void MeshSmoothing::SmoothMesh(void)
{
	UpdateNormal();
	UpdateVertex();
	UpdateVolume();
}

void MeshSmoothing::UpdateNormal(void)
{
	auto Guassian_s = [](double x, double sigma_s) {return exp(-x * x / (2 * sigma_s * sigma_s)); };
	auto Guassian_r = [](double x, double sigma_r) {return exp(-x * x / (2 * sigma_r * sigma_r)); };
	double coef;
	double average_dist;
	for (size_t i = 0; i < K1; i++)
	{
	
		for (auto fh : mesh.faces())
		{
			average_dist = 0.0;
			for (auto ff : fh.faces())
			{
				average_dist += (mesh.property(FCentroid, fh) - mesh.property(FCentroid, ff)).norm();
			}
			average_dist /= fh.valence();
			OpenMesh::Vec3d new_normal = OpenMesh::Vec3d(0.0, 0.0, 0.0);
			double Kp = 0.0;
			for (auto ff : fh.faces())
			{
				coef = mesh.property(FArea, ff) * Guassian_s((mesh.property(FCentroid, ff) - mesh.property(FCentroid, fh)).norm(), average_dist)
					* Guassian_r((mesh.normal(ff) - mesh.normal(fh)).norm(), 0.2);
				Kp += coef;
				new_normal += coef * mesh.normal(ff);
			}
			new_normal /= Kp;
			mesh.set_normal(fh, new_normal.normalize());
		}
	}
}

void MeshSmoothing::UpdateVertex(void)
{
	for (size_t i = 0; i < K2; i++)
	{
		for (auto vh : mesh.vertices())
		{
			int adjFacesNum = 0;
			OpenMesh::Vec3d pos = mesh.point(vh);
			OpenMesh::Vec3d offset = OpenMesh::Vec3d(0.0, 0.0, 0.0);
			for (auto vf : vh.faces())
			{
				offset += mesh.normal(vf) * (mesh.normal(vf) | (mesh.property(FCentroid, vf) - pos));
				adjFacesNum++;
			}
			mesh.set_point(vh, mesh.point(vh) + offset / adjFacesNum);
		}
	}
}

void MeshSmoothing::InitFaceProp()
{
	mesh.add_property(FArea);	
	mesh.add_property(FCentroid);
	std::vector<OpenMesh::Vec3d> facepoints(3);
	for (auto fh : mesh.faces())
	{
		mesh.property(FArea, fh) = mesh.calc_face_area(fh);
		mesh.property(FCentroid, fh) = mesh.calc_centroid(fh);
	}
}

double MeshSmoothing::CalculateVolume(void)
{
	double volume = 0.0;
	for (auto fh : mesh.faces())
	{
		std::vector<OpenMesh::Vec3d> facepoints;
		for (auto fv : fh.vertices())
		{
			facepoints.push_back(mesh.point(fv));
		}
		volume += ((facepoints[0] % facepoints[1]) | facepoints[2]) / 6.0;
	}
	return volume;
}

void MeshSmoothing::UpdateVolume(void)
{
	double current_volume = CalculateVolume();
	double scale = std::pow(original_volume / current_volume, 1.0 / 3.0);
	for (auto vh : mesh.vertices())
	{
		mesh.set_point(vh, scale * mesh.point(vh));
	}
}