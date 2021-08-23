#include "HW/MeshInterpolation.h"
#include <iostream>
# define PI 3.14159265358979323846

MeshInterpolation::MeshInterpolation(Mesh& m0, Mesh& m1) :mesh0(m0), mesh1(m1)
{
	Init();
}

MeshInterpolation::~MeshInterpolation()
{
	mesh0.remove_property(Jacobi0To1);
	mesh0.remove_property(DesiredJacobi);
	mesh0.remove_property(Rotation);
	mesh0.remove_property(Scale);
}

void MeshInterpolation::DoMorphing(double t)
{
	mesh0.add_property(DesiredJacobi);
	
	LocalStep(t);
	GlobalStep();
}

void MeshInterpolation::LocalStep(double t)
{	
	Eigen::Matrix2d rotation_interp, scale_interp;
	double theta, theta_interp;
	for (auto fh : mesh0.faces())
	{	
		if (mesh0.property(Rotation, fh)(1, 0) >= 0)
		{
			theta = acos(mesh0.property(Rotation, fh)(0, 0));
			theta_interp = t * theta;
		}
		else
		{
			theta_interp = 2 * PI - t * acos(mesh0.property(Rotation, fh)(0, 0));
		}
		
		
		rotation_interp << cos(theta_interp), -sin(theta_interp), sin(theta_interp), cos(theta_interp);
		scale_interp = (1 - t) * Eigen::Matrix2d::Identity() + t * mesh0.property(Scale, fh);
		mesh0.property(DesiredJacobi, fh) = rotation_interp * scale_interp;
	}
}

void MeshInterpolation::GlobalStep(void)
{
	int nf = mesh0.n_faces();
	int nv = mesh0.n_vertices();
	Eigen::VectorXd b(4 * nf + 2);
	Eigen::Matrix2d Desired;
	std::vector<int> face_points;
	for (auto fh : mesh0.faces())
	{
		Desired = mesh0.property(DesiredJacobi, fh);
		for (auto fv : fh.vertices())
		{
			face_points.push_back(fv.idx());
		}
		b(4 * fh.idx()) = Desired(0, 0);
		b(4 * fh.idx() + 1) = Desired(0, 1);
		b(4 * fh.idx() + 2) = Desired(1, 0);
		b(4 * fh.idx() + 3) = Desired(1, 1);
		face_points.clear();
	}
	b(4 * nf) = 0.0;
	b(4 * nf + 1) = 0.0;
	Eigen::VectorXd righthand = CoefMatA.transpose() * b;
	Eigen::VectorXd new_pos = lu.solve(righthand);
	int id;
	for (auto vh : mesh0.vertices())
	{
		id = vh.idx();
		mesh0.set_point(vh, OpenMesh::Vec3d(new_pos[id], new_pos[nv + id], 0.0));
	}
}

void MeshInterpolation::Init(void)
{
	
	mesh0.add_property(Jacobi0To1);
	mesh0.add_property(Rotation);
	mesh0.add_property(Scale);

	int nv = mesh0.n_vertices();
	int nf = mesh0.n_faces();
	typedef Eigen::Triplet<double> triple;
	std::vector<triple> TripletLists;
	CoefMatA.resize(4 * nf + 2, 2 * nv);

	Eigen::Matrix2d jacobi;
	Eigen::Matrix<double, 2, 3> coef;
	std::vector<int> face_points;
	double area;
	OpenMesh::Vec3d p0, p1, p2;
	for (auto fh : mesh0.faces())
	{
		//Init Face Property
		area = mesh0.calc_face_area(fh);
		for (auto fv : fh.vertices())
		{
			face_points.push_back(fv.idx());
		}
		coef(0, 0) = mesh0.point(mesh0.vertex_handle(face_points[1]))[1] - mesh0.point(mesh0.vertex_handle(face_points[2]))[1];
		coef(1, 0) = mesh0.point(mesh0.vertex_handle(face_points[2]))[0] - mesh0.point(mesh0.vertex_handle(face_points[1]))[0];
		coef(0, 1) = mesh0.point(mesh0.vertex_handle(face_points[2]))[1] - mesh0.point(mesh0.vertex_handle(face_points[0]))[1];
		coef(1, 1) = mesh0.point(mesh0.vertex_handle(face_points[0]))[0] - mesh0.point(mesh0.vertex_handle(face_points[2]))[0];
		coef(0, 2) = mesh0.point(mesh0.vertex_handle(face_points[0]))[1] - mesh0.point(mesh0.vertex_handle(face_points[1]))[1];
		coef(1, 2) = mesh0.point(mesh0.vertex_handle(face_points[1]))[0] - mesh0.point(mesh0.vertex_handle(face_points[0]))[0];
  		p0 = mesh1.point(mesh1.vertex_handle(face_points[0]));
		p1 = mesh1.point(mesh1.vertex_handle(face_points[1]));
		p2 = mesh1.point(mesh1.vertex_handle(face_points[2]));
		Eigen::Vector3d u(p0[0], p1[0], p2[0]);
		Eigen::Vector3d v(p0[1], p1[1], p2[1]);
		jacobi.row(0) = ((1 / (2 * area)) * coef * u).transpose();
		jacobi.row(1) = ((1 / (2 * area)) * coef * v).transpose();
		mesh0.property(Jacobi0To1, fh) = jacobi;

		//Init Rotation Scale Matrix between 0 to 1
		Eigen::Matrix2d U, V, Sigma;
		Eigen::JacobiSVD<Eigen::Matrix2d> svd(jacobi, Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		V = svd.matrixV();
		auto sigma_values = svd.singularValues();
		Sigma << sigma_values[0], 0, 0, sigma_values[1];
		if ((U * V.transpose()).determinant() < 0)
		{
			V.col(1) = -V.col(1);
		}
		mesh0.property(Rotation, fh) = U * V.transpose();
		mesh0.property(Scale, fh) = V * Sigma * V.transpose();

		//Init CoefMat A			
		TripletLists.emplace_back(4 * fh.idx(), face_points[0], (1 / (2 * area)) * coef(0, 0));
		TripletLists.emplace_back(4 * fh.idx(), face_points[1], (1 / (2 * area)) * coef(0, 1));
		TripletLists.emplace_back(4 * fh.idx(), face_points[2], (1 / (2 * area)) * coef(0, 2));

		TripletLists.emplace_back(4 * fh.idx() + 1, face_points[0], (1 / (2 * area)) * coef(1, 0));
		TripletLists.emplace_back(4 * fh.idx() + 1, face_points[1], (1 / (2 * area)) * coef(1, 1));
		TripletLists.emplace_back(4 * fh.idx() + 1, face_points[2], (1 / (2 * area)) * coef(1, 2));

		TripletLists.emplace_back(4 * fh.idx() + 2, nv + face_points[0], (1 / (2 * area)) * coef(0, 0));
		TripletLists.emplace_back(4 * fh.idx() + 2, nv + face_points[1], (1 / (2 * area)) * coef(0, 1));
		TripletLists.emplace_back(4 * fh.idx() + 2, nv + face_points[2], (1 / (2 * area)) * coef(0, 2));

		TripletLists.emplace_back(4 * fh.idx() + 3, nv + face_points[0], (1 / (2 * area)) * coef(1, 0));
		TripletLists.emplace_back(4 * fh.idx() + 3, nv + face_points[1], (1 / (2 * area)) * coef(1, 1));
		TripletLists.emplace_back(4 * fh.idx() + 3, nv + face_points[2], (1 / (2 * area)) * coef(1, 2));
		
		face_points.clear();
	}
	//fix a point
	TripletLists.emplace_back(4 * nf, 0, 1.0);
	TripletLists.emplace_back(4 * nf + 1, nv, 1.0);
	CoefMatA.setFromTriplets(TripletLists.begin(), TripletLists.end());
	//std::cout << CoefMatA << std::endl << std::endl;
	CoefMat = CoefMatA.transpose() * CoefMatA;
	//std::cout << CoefMat << std::endl;
	lu.compute(CoefMat);
}


