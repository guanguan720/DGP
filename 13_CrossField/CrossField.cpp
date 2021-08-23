#include <Eigen/Sparse>
#include <iostream>
#include "CrossField.h"

CrossField::CrossField(Mesh& m) :mesh(m)
{
	Init();
}

CrossField::~CrossField()
{
	mesh.remove_property(centroid);
	mesh.remove_property(local_x_axis);
	mesh.remove_property(local_y_axis);
	mesh.remove_property(local_halfedge_4);
}

std::vector<OpenMesh::Vec3d> CrossField::GetCrossField(void)
{
	Eigen::VectorXd sol = Calc_uf4();
	std::vector<OpenMesh::Vec3d> field(4 * mesh.n_faces());
	std::complex<double> uf_1, uf_2, uf_3, uf_4;
	double len;
	for (auto& fh : mesh.faces())
	{
		std::complex<double> uf4(sol[2 * fh.idx()], sol[2 * fh.idx() + 1]);
		uf_1 = std::pow(uf4, 0.25);
		uf_2 = -uf_1;
		uf_3 = std::complex<double>(-uf_1.imag(), uf_1.real());
		uf_4 = -uf_3;
		len = std::numeric_limits<double>::max();		
		for (auto& eh : fh.edges())
		{
			if (mesh.calc_edge_length(eh) < len)
			{
				len = mesh.calc_edge_length(eh);
			}
		}
		field[4 * fh.idx()] = mesh.property(centroid, fh) + (uf_1.real() * mesh.property(local_x_axis, fh) + uf_1.imag() * mesh.property(local_y_axis, fh)).normalize() * 0.2 * len;
		field[4 * fh.idx() + 1] = mesh.property(centroid, fh) + (uf_2.real() * mesh.property(local_x_axis, fh) + uf_2.imag() * mesh.property(local_y_axis, fh)).normalize() * 0.2 * len;
		field[4 * fh.idx() + 2] = mesh.property(centroid, fh) + (uf_3.real() * mesh.property(local_x_axis, fh) + uf_3.imag() * mesh.property(local_y_axis, fh)).normalize() * 0.2 * len;
		field[4 * fh.idx() + 3] = mesh.property(centroid, fh) + (uf_4.real() * mesh.property(local_x_axis, fh) + uf_4.imag() * mesh.property(local_y_axis, fh)).normalize() * 0.2 * len;
	}
	return field;
}

void CrossField::Init(void)
{
	//initialize face property
	mesh.add_property(centroid);
	mesh.add_property(local_x_axis);
	mesh.add_property(local_y_axis);
	std::vector<OpenMesh::Vec3d> p;
	for (auto& fh : mesh.faces())
	{
		for (auto& fv : fh.vertices())
		{
			p.push_back(mesh.point(fv));
		}	
		mesh.property(centroid, fh) = (p[0] + p[1] + p[2]) / 3;
		mesh.property(local_x_axis, fh) = (p[1] - p[0]).normalize();
		mesh.property(local_y_axis, fh) = ((p[2] - p[0]) - ((p[2] - p[0]) | mesh.property(local_x_axis, fh)) * mesh.property(local_x_axis, fh)).normalize();
		p.clear();
	}
	//initialize halfedge property
	mesh.add_property(local_halfedge_4);
	OpenMesh::Vec3d to, from;
	for (auto& he : mesh.halfedges())
	{	
		if (!mesh.is_boundary(he))
		{
			auto f = he.face();
			to = mesh.point(he.to());
			from = mesh.point(he.from());
			auto x = mesh.property(local_x_axis, f);
			auto y = mesh.property(local_y_axis, f);
			std::complex<double> he_complex((to - from) | mesh.property(local_x_axis, f), (to - from) | mesh.property(local_y_axis, f));
			mesh.property(local_halfedge_4, he) = std::pow(conj(he_complex), 4);
		}
	}
}  

Eigen::VectorXd CrossField::Calc_uf4(void)
{
	//assembly coefficient matrix
	uf_4.resize(mesh.n_faces());
	typedef Eigen::Triplet<double> triple;
	std::vector<triple> TripletList;
	int interior_edge = 0;
	std::complex<double> ef4, eg4;
	for (auto& eh : mesh.edges())
	{
		if (!mesh.is_boundary(eh))
		{
			auto h0 = eh.h0();
			auto h1 = eh.h1();
			auto f = h0.face();
			auto g = h1.face();
			ef4 = mesh.property(local_halfedge_4, h0);
			eg4 = mesh.property(local_halfedge_4, h1);
			TripletList.emplace_back(2 * interior_edge, 2 * f.idx(), ef4.real());
			TripletList.emplace_back(2 * interior_edge, 2 * f.idx() + 1, -ef4.imag());
			TripletList.emplace_back(2 * interior_edge, 2 * g.idx(), -eg4.real());
			TripletList.emplace_back(2 * interior_edge, 2 * g.idx() + 1, eg4.imag());

			TripletList.emplace_back(2 * interior_edge + 1, 2 * f.idx(), ef4.imag());
			TripletList.emplace_back(2 * interior_edge + 1, 2 * f.idx() + 1, ef4.real());
			TripletList.emplace_back(2 * interior_edge + 1, 2 * g.idx(), -eg4.imag());
			TripletList.emplace_back(2 * interior_edge + 1, 2 * g.idx() + 1, -eg4.real());
			interior_edge++;
		}
	}
	//fix a cross field in 0-th face
	TripletList.emplace_back(2 * interior_edge, 0, 1.0);
	TripletList.emplace_back(2 * interior_edge + 1, 1, 1.0);
	Eigen::SparseMatrix<double> A(2 * (interior_edge + 1), 2 * mesh.n_faces());
	A.setFromTriplets(TripletList.begin(), TripletList.end());
	Eigen::SparseMatrix<double> coef = A.transpose() * A;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> lu(coef);
	//assembly righthand term
	Eigen::VectorXd b(2 * (interior_edge + 1));
	b.setZero();	
	b[2 * interior_edge] = 1;
	Eigen::VectorXd righthand = A.transpose() * b;
	Eigen::VectorXd sol = lu.solve(righthand);
	return sol;
}