#include"src/HW/Deformation.h"
#include <iostream>

Deformation::Deformation(Mesh& m):mesh(m)
{
	
	Calc_Cot();	
	InitProp();
}

Deformation::~Deformation()
{
	mesh.remove_property(cot_value);
	mesh.remove_property(R);
}

void Deformation::DoDeformation(int iterNum, std::vector<int> Fixed, OpenMesh::Vec3d Active)
{	
	FixedPoints = Fixed;
	ActivePos = Active;
	Calc_CoefMat();
	for (size_t i = 0; i < iterNum; i++)
	{
		UpdateRotation();
		UpadatePos();		
	}
}

void Deformation::UpdateRotation(void)
{
	Eigen::Matrix3d S = Eigen::Matrix3d::Zero();
	for (auto vh : mesh.vertices())
	{
		double wij, min;
		Eigen::Vector3d eij, eij_prime;
		int i = vh.idx();
		int j;
		for (auto voh : vh.outgoing_halfedges())
		{
			auto vv = voh.to();
			j = vv.idx();
			wij = 0.5 * (mesh.property(cot_value, voh) + mesh.property(cot_value, voh.opp()));
			eij[0] = (LastFramePos[i] - LastFramePos[j])[0];
			eij[1] = (LastFramePos[i] - LastFramePos[j])[1];
			eij[2] = (LastFramePos[i] - LastFramePos[j])[2];
			eij_prime[0] = (mesh.point(vh) - mesh.point(vv))[0];
			eij_prime[1] = (mesh.point(vh) - mesh.point(vv))[1];
			eij_prime[2] = (mesh.point(vh) - mesh.point(vv))[2];
			S += wij * eij * eij_prime.transpose();
		}
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Matrix3d Rv = svd.matrixV() * svd.matrixU().transpose();
		if (Rv.determinant() < 0)
		{
			/*for (size_t i = 0; i < 3; i++)
			{
				if (std::abs(singular_values[i] - min) < 1e-5)
				{
					minele_col = i;
					break;
				}
			}*/
			Rv.col(2) = -Rv.col(2);
		}
		mesh.property(R, vh) = Rv;
	}
}


void Deformation::UpadatePos(void)
{
	//assembly right-hand term
	Eigen::MatrixXd b = Eigen::MatrixXd::Zero(mesh.n_vertices(), 3);
	double wij;
	Eigen::Matrix3d Ri, Rj;
	Eigen::Vector3d pi, pj;
	int i, j;
	for (auto he : mesh.halfedges())
	{
		if (!mesh.is_boundary(he))
		{
			auto vi = he.from();
			auto vj = he.to();
			i = vi.idx();
			j = vj.idx();
			Ri = mesh.property(R, vi);
			Rj = mesh.property(R, vj);
			wij = 0.5 * mesh.property(cot_value, he);
			pi[0] = LastFramePos[i][0]; pi[1] = LastFramePos[i][1]; pi[2] = LastFramePos[i][2];
			pj[0] = LastFramePos[j][0]; pj[1] = LastFramePos[j][1]; pj[2] = LastFramePos[j][2];
			if (std::find(FixedPoints.begin(), FixedPoints.end(), i) == FixedPoints.end())
			{
				b.row(i) += 0.5 * wij * (Ri + Rj) * (pi - pj);
			}
			if (std::find(FixedPoints.begin(), FixedPoints.end(), j) == FixedPoints.end())
			{
				b.row(j) += 0.5 * wij * (Ri + Rj) * (pj - pi);
			}
		}
	}
	Eigen::Vector3d tmp;
	for (size_t i = 0; i < FixedPoints.size() - 1; i++)
	{
		tmp[0] = LastFramePos[FixedPoints[i]][0];
		tmp[1] = LastFramePos[FixedPoints[i]][1];
		tmp[2] = LastFramePos[FixedPoints[i]][2];
		b.row(FixedPoints[i]) = tmp;
	}
	b.row(FixedPoints[FixedPoints.size() - 1]) = Eigen::Vector3d(ActivePos[0], ActivePos[1], ActivePos[2]);
	//solve new position
	Eigen::MatrixXd new_pos(mesh.n_vertices(), 3);
	new_pos.col(0) = solver.solve(b.col(0));
	new_pos.col(1) = solver.solve(b.col(1));
	new_pos.col(2) = solver.solve(b.col(2));
	//std::cout << "new pos:" << new_pos << std::endl << std::endl;
	for (auto vh : mesh.vertices())
	{
		OpenMesh::Vec3d pos(new_pos.row(vh.idx())[0], new_pos.row(vh.idx())[1], new_pos.row(vh.idx())[2]);
		mesh.set_point(vh, pos);
	}
}

void Deformation::InitProp(void)
{
	mesh.add_property(R);
	LastFramePos.resize(mesh.n_vertices());
	for (auto vh : mesh.vertices())
	{
		LastFramePos[vh.idx()] = mesh.point(vh);
	}
}

void Deformation::Calc_Cot(void)
{
	mesh.add_property(cot_value);
	for (auto he : mesh.halfedges())
	{
		OpenMesh::Vec3d AnglePoint, Point1, Point2;
		if (! mesh.is_boundary(he))
		{
			Point1 = mesh.point(he.from());
			Point2 = mesh.point(he.to());
			AnglePoint = mesh.point(he.next().to());
			mesh.property(cot_value, he) = cot(AnglePoint, Point1, Point2);
		}
		else
		{
			mesh.property(cot_value, he) = 0.0;
		}
	}
}

void Deformation::Calc_CoefMat(void)
{
	typedef Eigen::Triplet<double> triple;
	std::vector<triple> TripletList;
	int n = mesh.n_vertices();
	CoefMat.resize(n, n);
	double cot;
	int i, j;
	for (auto he : mesh.halfedges())
	{
		if (!mesh.is_boundary(he))
		{
			i = he.from().idx();
			j = he.to().idx();
			cot = mesh.property(cot_value, he);
			if (std::find(FixedPoints.begin(), FixedPoints.end(), i) == FixedPoints.end())
			{
				TripletList.emplace_back(i, i, 0.5 * cot);
				TripletList.emplace_back(i, j, -0.5 * cot);
			}	
			if (std::find(FixedPoints.begin(), FixedPoints.end(), j) == FixedPoints.end())
			{
				TripletList.emplace_back(j, j, 0.5 * cot);
				TripletList.emplace_back(j, i, -0.5 * cot);
			}		
		}		
	}
	for (size_t i = 0; i < FixedPoints.size(); i++)
	{
		TripletList.emplace_back(FixedPoints[i], FixedPoints[i], 1.0);
	}
	CoefMat.setFromTriplets(TripletList.begin(), TripletList.end());
	//std::cout << "Matrix:" << CoefMat << std::endl << std::endl;
	solver.compute(CoefMat);
}

