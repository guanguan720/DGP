#include"ARAP.h"
#include <iostream>

ARAP::ARAP(Mesh& m):mesh(m)
{
	CalcCot();
	InitFaceProp();
	InitGlobalMat();
	InitPara();
}

ARAP::~ARAP()
{
	mesh.remove_property(face_area);
	mesh.remove_property(jacobi);
	mesh.remove_property(rotationMat);
	mesh.remove_property(ori_face_uv);
	mesh.remove_property(gradientMat);
}

void ARAP::DoARAP(void)
{
	int n = mesh.n_faces();	
	double energy = std::numeric_limits<double>::max();
	double pre_energy = std::numeric_limits<double>::max();
	do
	{
		pre_energy = energy;
		LocalStep();
		GlobalStep();
		energy = getEnergy();
		double diff = energy - pre_energy;
	} while ((pre_energy- energy) > 1e-3);
}

void ARAP::InitPara(void)
{
	std::unique_ptr<TutteEmbedding> tutte(new TutteEmbedding(mesh));
	tutte->DoTutte();
	OpenMesh::Vec3d b0_3d = mesh.point(mesh.vertex_handle(0));
	b0 = Eigen::Vector2d(b0_3d[0], b0_3d[1]);
}

void ARAP::LocalStep(void)
{
	Eigen::RowVector3d u;
	Eigen::RowVector3d v;
	Eigen::Matrix2d U;
	Eigen::Matrix2d V;
	Eigen::Matrix<double, 3, 2> coefMat;
	int cout;
	for (auto fh : mesh.faces())
	{
		cout = 0;
		for (auto fv : fh.vertices())
		{
			u(cout) = mesh.point(fv)[0];
			v(cout) = mesh.point(fv)[1];
			cout++;
		}
		coefMat = mesh.property(gradientMat, fh);
		mesh.property(jacobi, fh).row(0) = u * coefMat;
		mesh.property(jacobi, fh).row(1) = v * coefMat;
		Eigen::JacobiSVD<Eigen::Matrix2d> svd(mesh.property(jacobi,fh), Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		V = svd.matrixV();
		if ((U * V.transpose()).determinant() < 0)
		{
			V.col(1) = -V.col(1);
		}
		mesh.property(rotationMat, fh) = U * V.transpose();
	}
}

// Jacobi form given in paper
//void ARAP::LocalStep(void)
//{
//	Eigen::Matrix<double, 3, 2> original_uv;
//	Eigen::Matrix<double, 3, 2> uv;
//	std::vector<OpenMesh::Vec3d> v_pts(3);
//	Eigen::Matrix2d St = Eigen::Matrix2d::Zero();
//	Eigen::Matrix2d U;
//	Eigen::Matrix2d V;
//	for (auto fh : mesh.faces())
//	{
//		original_uv = mesh.property(ori_face_uv, fh);
//		int cout = 0;
//		for (auto fv : fh.vertices())
//		{
//			uv(cout, 0) = mesh.point(fv)[0];
//			uv(cout, 1) = mesh.point(fv)[1];
//			v_pts[cout] = mesh.point(fv);
//			cout++;
//		}
//		for (size_t i = 0; i < 3; i++)
//		{
//			double c = cot(v_pts[(i + 2) % 3], v_pts[i % 3], v_pts[(i + 1) % 3]);
//			St += c * (uv.row(i) - uv.row((i + 1) % 3)).transpose() * (original_uv.row(i) - original_uv.row((i + 1) % 3));
//		}
//		mesh.property(jacobi, fh) = St;
//		Eigen::JacobiSVD<Eigen::Matrix2d> svd(mesh.property(jacobi, fh), Eigen::ComputeFullU | Eigen::ComputeFullV);
//		U = svd.matrixU();
//		V = svd.matrixV();
//		if ((U * V.transpose()).determinant() < 0)
//		{
//			V.col(1) = -V.col(1);
//		}
//		mesh.property(rotationMat, fh) = U * V.transpose();
//	}
//}

void ARAP::GlobalStep(void)
{
	int n = mesh.n_vertices();
	GenRighthandTerm();
 	Eigen::MatrixXd uv = solver.solve(b);
	//std::cout << "uv:" << uv << std::endl << std::endl;
	for (size_t i = 0; i < n; i++)
	{
		auto v = mesh.vertex_handle(i);
		mesh.set_point(v, OpenMesh::Vec3d(uv(i, 0), uv(i, 1), 0.0));
	}
}

void ARAP::InitFaceProp(void)
{
	mesh.add_property(face_area);
	mesh.add_property(jacobi);
	mesh.add_property(rotationMat);
	mesh.add_property(gradientMat);
	mesh.add_property(ori_face_uv);
	double area;
	std::vector<OpenMesh::Vec3d> face_points;
	for (auto fh : mesh.faces())
	{
		//face_area
		area = mesh.calc_face_area(fh);
		mesh.property(face_area, fh) = area;
		//ori_face_uv
		for (auto fv : fh.vertices())
		{
			face_points.push_back(mesh.point(fv));
		}
		mesh.property(ori_face_uv, fh).row(0) = Eigen::RowVector2d(0.0, 0.0);
		mesh.property(ori_face_uv, fh).row(1) = Eigen::RowVector2d((face_points[1] - face_points[0]).norm(), 0.0);

		Eigen::Vector3d e1((face_points[1] - face_points[0]).norm(), 0.0, 0.0);
		Eigen::Vector3d e2((face_points[2] - face_points[0]).norm() * cos(face_points[0], face_points[1], face_points[2]), 
							(face_points[2] - face_points[0]).norm() * sin(face_points[0], face_points[1], face_points[2]), 0.0);
		/*Eigen::Vector3d cross = e1.cross(e2);
		if (cross[2] > 0)
		{*/
			mesh.property(ori_face_uv, fh).row(2) = Eigen::RowVector2d((face_points[2] - face_points[0]).norm() * cos(face_points[0], face_points[1], face_points[2]),
																		(face_points[2] - face_points[0]).norm() * sin(face_points[0], face_points[1], face_points[2]));
		/*}
		else
		{
			mesh.property(ori_face_uv, fh).row(2) = Eigen::RowVector2d((face_points[2] - face_points[0]).norm() * cos(face_points[0], face_points[1], face_points[2]),
																		-(face_points[2] - face_points[0]).norm() * sin(face_points[0], face_points[1], face_points[2]));
		}*/
		face_points.clear();
		//gradientMat
		Eigen::Matrix<double, 3, 2> ori_uv = mesh.property(ori_face_uv, fh);
		Eigen::Matrix<double, 3, 2> tmp;
		tmp(0, 0) = ori_uv(1, 1) - ori_uv(2, 1);
		tmp(1, 0) = ori_uv(2, 1);
		tmp(2, 0) = -ori_uv(1, 1);
		tmp(0, 1) = ori_uv(2, 0) - ori_uv(1, 0);
		tmp(1, 1) = -ori_uv(2, 0);
		tmp(2, 1) = ori_uv(1, 0);
		mesh.property(gradientMat, fh) = (1 / (2 * area)) * tmp;
	}
}

void ARAP::InitGlobalMat(void)
{
	int n = mesh.n_vertices();
	globalMat.resize(n, n);
	typedef Eigen::Triplet<double> triple;
	std::vector<triple> tripletList;
	OpenMesh::Vec3d vi, vj, vij, vji;
	double cij, cji;
	int idxi, idxj;
	double coeff_i, coeff_j;
	for (auto vh : mesh.vertices())
	{
		idxi = vh.idx();
		if (idxi == 0)
		{
			tripletList.emplace_back(0, 0, 1.0);
			continue;
		}
		vi = mesh.point(vh);
		coeff_i = 0.0;
		for (auto voh : vh.outgoing_halfedges())
		{
			coeff_j = 0.0;
			auto vv = voh.to();		
			vj = mesh.point(vv);		
			cij = cot_value[voh.idx()];
			coeff_i += cij;
			coeff_j -= cij;
			
			cji = cot_value[voh.opp().idx()];
			coeff_i += cji;
			coeff_j -= cji;			
			idxj = vv.idx();  
			tripletList.emplace_back(idxi, idxj, coeff_j);		
		}
		tripletList.emplace_back(idxi, idxi, coeff_i);
	}
	globalMat.setFromTriplets(tripletList.begin(), tripletList.end());
	//std::cout << "A:" << globalMat << std::endl << std::endl;
	solver.compute(globalMat);
}

void ARAP::CalcCot(void)
{
	std::vector<double> tmp(mesh.n_halfedges(), 0);
	cot_value = tmp;
	for (auto he : mesh.halfedges())
	{
		if (!mesh.is_boundary(he))
		{
			auto vi = he.from();
			auto vj = he.to();
			auto vij = he.next().to();
			cot_value[he.idx()] = cot(mesh.point(vij), mesh.point(vi), mesh.point(vj));
		}
	}
}

void ARAP::GenRighthandTerm(void)
{
	int n = mesh.n_vertices();
	b = Eigen::MatrixXd::Zero(n, 2);
	int idxi, idxj, cout;
	double cij, cji;
	OpenMesh::Vec3d vi, vj, vij, vji;
	Eigen::Vector2d xi, xj;
	for (auto vh : mesh.vertices())
	{
		idxi = vh.idx();
		if (idxi == 0)
		{
			b.row(0) = b0.transpose();
			continue;
		}
		for (auto voh : vh.outgoing_halfedges())
		{
			vi = mesh.point(vh);
				
			auto vv = voh.to();
			idxj = vv.idx();
			vj = mesh.point(vv);
			if (!mesh.is_boundary(voh))
			{				
				auto vvf = voh.face();				
				cij = cot_value[voh.idx()];
				cout = 0;
				for (auto fv : vvf.vertices())
				{
					if (idxi == fv.idx())
					{
						xi = mesh.property(ori_face_uv, vvf).row(cout).transpose();
					}
					if (idxj == fv.idx())
					{
						xj = mesh.property(ori_face_uv, vvf).row(cout).transpose();
					}
					cout++;
				}
				b.row(idxi) += cij * mesh.property(rotationMat, vvf) * (xi - xj);
			}			
			if (!mesh.is_boundary(voh.opp()))
			{
				auto voppf = voh.opp().face();
				cji = cot_value[voh.opp().idx()];
				cout = 0;
				for (auto fv : voppf.vertices())
				{
					if (idxi == fv.idx())
					{
						xi = mesh.property(ori_face_uv, voppf).row(cout).transpose();
					}
					if (idxj == fv.idx())
					{
						xj = mesh.property(ori_face_uv, voppf).row(cout).transpose();
					}
					cout++;
				}
				b.row(idxi) += cji * mesh.property(rotationMat, voppf) * (xi - xj);
			}
		}
	}
	
	//std::cout << "b:" << b << std::endl << std::endl;
}

double ARAP::getEnergy(void)
{
	double energy = 0.0;
	for (auto fh : mesh.faces())
	{
		energy += mesh.property(face_area, fh) * (mesh.property(rotationMat, fh) - mesh.property(jacobi, fh)).norm();
	}
	return energy;
}

