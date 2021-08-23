#include "HW/MeanValueCoordinatePara.h"
#include <iostream>
#define PI 3.14159265358979323846

MeanValuePara::~MeanValuePara()
{

}

void MeanValuePara::DoPara(void)
{
	CalcBoundary();
	CalcInteriorPos();
}

void MeanValuePara::CalcBoundary(void)
{
	auto he = mesh.halfedges_begin();
	while (!mesh.is_boundary(*he)) he++;
	boundary_start = *he;
	auto boundary_he = *he;
	boundaryNum = 1;
	while (boundary_he.next() != boundary_start)
	{
		boundaryNum++;		
		boundary_he = boundary_he.next();
	}
	b.resize(mesh.n_vertices(), 3);
	b.setZero();
	for (size_t i = 0; i < boundaryNum; i++)
	{
		auto boundary_v = boundary_he.to();
		b.row(boundary_v.idx()) = Eigen::RowVector3d(std::cos(2 * i * PI / boundaryNum), std::sin(2 * PI * i / boundaryNum), 0.0);
		boundary_he = boundary_he.next();
	}
}

void MeanValuePara::CalcInteriorPos(void)
{
	//Coefficient Matrix
	typedef Eigen::Triplet<double> triple;
	int n = mesh.n_vertices();
	Eigen::SparseMatrix<double> coefMat(n, n);
	std::vector<triple> TripletLists;

	std::vector<OpenMesh::Vec3d> face_points;
	std::vector<int> points_idx;
	std::vector<int> boundary_tags;
	double t;
	for (auto fh : mesh.faces()) // i think traverse halfedges also work
	{
		for (auto fv : fh.vertices())
		{
			face_points.push_back(mesh.point(fv));
			boundary_tags.push_back(mesh.is_boundary(fv));
			points_idx.push_back(fv.idx());
		}
		if (boundary_tags[0] == 0)
		{
			t = tan_half(face_points[0], face_points[1], face_points[2]);
			TripletLists.emplace_back(points_idx[0], points_idx[1], t / (face_points[0] - face_points[1]).norm());
			TripletLists.emplace_back(points_idx[0], points_idx[0], -t / (face_points[0] - face_points[1]).norm());
			TripletLists.emplace_back(points_idx[0], points_idx[2], t / (face_points[0] - face_points[2]).norm());
			TripletLists.emplace_back(points_idx[0], points_idx[0], -t / (face_points[0] - face_points[2]).norm());
		}
		if (boundary_tags[1] == 0)
		{
			t = tan_half(face_points[1], face_points[0], face_points[2]);
			TripletLists.emplace_back(points_idx[1], points_idx[0], t / (face_points[1] - face_points[0]).norm());
			TripletLists.emplace_back(points_idx[1], points_idx[1], -t / (face_points[1] - face_points[0]).norm());
			TripletLists.emplace_back(points_idx[1], points_idx[2], t / (face_points[1] - face_points[2]).norm());
			TripletLists.emplace_back(points_idx[1], points_idx[1], -t / (face_points[1] - face_points[2]).norm());
		}
		if (boundary_tags[2] == 0)
		{
			t = tan_half(face_points[2], face_points[0], face_points[1]);
			TripletLists.emplace_back(points_idx[2], points_idx[0], t / (face_points[2] - face_points[0]).norm());
			TripletLists.emplace_back(points_idx[2], points_idx[2], -t / (face_points[2] - face_points[0]).norm());
			TripletLists.emplace_back(points_idx[2], points_idx[1], t / (face_points[2] - face_points[1]).norm());
			TripletLists.emplace_back(points_idx[2], points_idx[2], -t / (face_points[2] - face_points[1]).norm());
		}
		face_points.clear();
		points_idx.clear();
		boundary_tags.clear();
	}
	
	auto he = boundary_start;
	for (size_t i = 0; i < boundaryNum; i++)
	{
		auto v = he.to();
		TripletLists.emplace_back(v.idx(), v.idx(), 1.0);
		he = he.next();
	}
	coefMat.setFromTriplets(TripletLists.begin(), TripletLists.end());

	//set new position
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(coefMat);
	Eigen::MatrixXd new_pos(n, 3);
	new_pos.col(0) = solver.solve(b.col(0));
	new_pos.col(1) = solver.solve(b.col(1));
	int id;
	for (auto vh : mesh.vertices())
	{
		id = vh.idx();
		mesh.set_point(vh, OpenMesh::Vec3d(new_pos.row(id)[0], new_pos.row(id)[1], 0.0));
	}
}