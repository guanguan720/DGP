#include "HW\TutteEmbedding.h"

#define PI 3.14159265358979323846

TutteEmbedding::TutteEmbedding(Mesh& m):mesh(m)
{

}

TutteEmbedding::~TutteEmbedding()
{

}

void TutteEmbedding::DoTutte(void)
{
	BindBoundToCircle();
	InnerVEmbedding();
}

void TutteEmbedding::BindBoundToCircle(void)
{
	auto he = mesh.halfedges_begin();
	while (! mesh.is_boundary(*he)) he++;
	auto boundary_he = *he;
	auto boundary_heStart = *he;
	int boundaryNum = 1;
	while (boundary_he.next() != boundary_heStart)
	{
		boundaryNum++;
		boundary_he = boundary_he.next();
	}
	for (size_t i = 0; i < boundaryNum; i++)
	{
		auto boundary_v = boundary_he.to();
		mesh.set_point(boundary_v, OpenMesh::Vec3d(cos((2 * PI * i) / (double)boundaryNum), sin((2 * PI * i) / (double)boundaryNum), 0.0));
		boundary_he = boundary_he.next();
	}
}

void TutteEmbedding::InnerVEmbedding(void)
{
	typedef Eigen::Triplet<double> triple;
	int n = mesh.n_vertices();
	Eigen::SparseMatrix<double> A(n, n);
	Eigen::VectorXd bu=Eigen::VectorXd::Zero(n);
	Eigen::VectorXd bv = Eigen::VectorXd::Zero(n);
	std::vector<triple> tripletList;
	int idx;

	//assemble A and b
	for (auto vh : mesh.vertices())
	{
		idx = vh.idx();
		if (mesh.is_boundary(vh))
		{
			tripletList.push_back(triple(idx, idx, 1));
			bu(idx) = mesh.point(vh)[0];
			bv(idx) = mesh.point(vh)[1];
		}
		else
		{
			tripletList.push_back(triple(idx, idx, vh.valence()));
			for (auto vv : vh.vertices())
			{
				if (mesh.is_boundary(vv))
				{
					bu(idx) += mesh.point(vv)[0];
					bv(idx) += mesh.point(vv)[1];
				}
				else
				{
					tripletList.push_back(triple(idx, vv.idx(), -1.0));
				}
			}
		}
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end());

	//solve uv
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(A);
	Eigen::VectorXd u = chol.solve(bu);
	Eigen::VectorXd v = chol.solve(bv);

	//set uv point
	for (auto vh : mesh.vertices())
	{
		int idx = vh.idx();
		mesh.set_point(vh, OpenMesh::Vec3d(u(idx), v(idx), 0));
	}
}