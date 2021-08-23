#include "HW/LloydIteration.h"



LloydIteration::LloydIteration(Mesh& m, int k, std::vector<int> InitSeed) :mesh(m)
{
	kNum = k;
	CurProxyNormal.resize(k);
	faceNum = mesh.n_faces();
	CurEnergy = 0.0;
	CurSeedIdx = InitSeed;
}

LloydIteration::~LloydIteration()
{
	mesh.remove_property(FaceCluster);
	mesh.remove_property(FaceArea);
}

void LloydIteration::GetApproximation(void)
{
	Init();
	double LastEnergy = std::numeric_limits<double>::max();
	int cout = 0;
	while ((cout < 100) && (abs((LastEnergy - CurEnergy) / CurEnergy) > 1e-3))
	{
		LastEnergy = CurEnergy;
		UpdatePartition();
		UpdateProxyFit();
		UpdateSeedTri();
		cout++;
	}

}

void LloydIteration::Init(void)
{
	//Init Proxy
	for (size_t i = 0; i < kNum; i++)
	{
		CurProxyNormal[i] = mesh.normal(mesh.face_handle(CurSeedIdx[i]));
	}
	//Compute Face Area
	mesh.add_property(FaceArea);
	for (auto fh : mesh.faces())
	{
		mesh.property(FaceArea, fh) = mesh.calc_face_area(fh);
	}
	//add property
	mesh.add_property(FaceCluster);
	/*for (auto& fh : mesh.faces())
	{
		mesh.property(FaceCluster, fh) = -1;
	}*/
}

void LloydIteration::UpdatePartition(void)
{
	//data for priority queue
	struct triangle
	{
		int FaceIdx;
		int ClusterIdx;
	};
	typedef std::pair<double, triangle> pair;
	struct cmp
	{
		bool operator()(pair a, pair b)
		{
			return a.first > b.first;
		}
	};
	std::priority_queue<pair, std::vector<pair>, cmp> queue;
	//set all faces' tags to be false  ** do NOT forget!!! **
	for (auto fh : mesh.faces())
	{
		mesh.status(fh).set_tagged(false);
	}
	//do partition
	int seedIdx;
	double FaceError;
	for (int i = 0; i < kNum; i++)
	{
		seedIdx = CurSeedIdx[i];
		auto fh = mesh.face_handle(seedIdx);
		mesh.property(FaceCluster, fh) = i;
		mesh.status(fh).set_tagged(true);		
		for (auto& ff : mesh.ff_range(fh))
		{
			if (!mesh.status(ff).tagged())			
			{
				FaceError = (mesh.normal(ff) - CurProxyNormal[i]).sqrnorm() * mesh.property(FaceArea, ff);
				queue.push(std::make_pair(FaceError, triangle{ ff.idx(), i }));
			}
		}
	}
	//grow the region
	pair best;
	while (!queue.empty())
	{
		best = queue.top();
		queue.pop(); //NOTICE the position!!!
		auto fh = mesh.face_handle(best.second.FaceIdx);
		if (!mesh.status(fh).tagged())	
		{
			mesh.property(FaceCluster, fh) = best.second.ClusterIdx;
			mesh.status(fh).set_tagged(true);			
			for (auto& ff : mesh.ff_range(fh))
			{
				if (!mesh.status(ff).tagged())				
				{
					FaceError = (mesh.normal(ff) - CurProxyNormal[best.second.ClusterIdx]).sqrnorm() * mesh.property(FaceArea, ff);
					queue.push(std::make_pair(FaceError, triangle{ ff.idx(), best.second.ClusterIdx }));
				}
			}
		}	
	}
}

void LloydIteration::UpdateProxyFit(void)
{
	CurProxyNormal.clear();
	for (size_t i = 0; i < kNum; i++)
	{
		CurProxyNormal.push_back(OpenMesh::Vec3d(0.0, 0.0, 0.0));
	}
	int clusterId;
	for (auto& fh : mesh.faces())
	{
		clusterId = mesh.property(FaceCluster, fh);
		CurProxyNormal[clusterId] += mesh.property(FaceArea, fh) * mesh.normal(fh);	
	}
	for (size_t i = 0; i < kNum; i++)
	{
		CurProxyNormal[i] = CurProxyNormal[i].normalize();
	}
}

void LloydIteration::UpdateSeedTri(void)
{
	typedef std::pair<double, int> pair;
	struct compare
	{
		bool operator()(pair a, pair b)
		{
			return a.first > b.first;
		}
	};
	std::vector<std::priority_queue<pair, std::vector<pair>, compare>> k_queue(kNum);
	int clusterId;
	double FaceError;
	CurEnergy = 0.0;
	for (auto& fh : mesh.faces())
	{	
		clusterId = mesh.property(FaceCluster, fh);
		FaceError = (mesh.normal(fh) - CurProxyNormal[clusterId]).sqrnorm() * mesh.property(FaceArea, fh);
		k_queue[clusterId].push(std::make_pair(FaceError, fh.idx()));
		//update current energy
		CurEnergy += FaceError;		
	}
	pair best;
	for (size_t i = 0; i < kNum; i++)
	{
		best = k_queue[i].top();
		CurSeedIdx[i] = best.second;
	}
}