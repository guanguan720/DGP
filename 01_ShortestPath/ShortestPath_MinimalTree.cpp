#include "HW\ShorestPath_MinimalTree.h"
#include <queue>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
using namespace std;

ShortestPath_MinimalTree::ShortestPath_MinimalTree()
{
	IsChoosing = true;
}

ShortestPath_MinimalTree::~ShortestPath_MinimalTree()
{
}

vector<int> ShortestPath_MinimalTree::findPath(Mesh& m, int a, int b, float& dist)
{
	typedef pair<float, int> pair;
	typedef tuple<float, int, int> triple; //shortest distance, idx, pre_footstep
	priority_queue<triple, vector<triple>, greater<triple>> U;
	vector<pair> S;
	vector<int> S_idx;
	int vNum = m.n_vertices();
	vector<int> FootPrint(vNum); //record which points the ShoretstPath pass through
	
	auto start_point = m.vertex_handle(a);
	S.push_back(pair(0.0, a));
	S_idx.push_back(a);
	FootPrint[a] = a;
	for (auto vv_ite = m.vv_iter(start_point); vv_ite.is_valid(); ++vv_ite)
	{
		float dist = (m.point(*vv_ite) - m.point(start_point)).norm();
		int id = (*vv_ite).idx();
		U.push(triple(dist, id, a));
	}

	triple nearestPoint;
	do{
		nearestPoint = U.top();
		if (find(S_idx.begin(),S_idx.end(),get<1>(nearestPoint)) == S_idx.end())
		{
			S.push_back(pair(get<0>(nearestPoint), get<1>(nearestPoint)));			
			FootPrint[get<1>(nearestPoint)] = get<2>(nearestPoint);
			S_idx.push_back(get<1>(nearestPoint));
			for (auto vv : m.vv_range(m.vertex_handle(get<1>(nearestPoint))))
			{
				// update the set U
				if (find(S_idx.begin(), S_idx.end(), vv.idx()) == S_idx.end())
				{
					float dist = (m.point(m.vertex_handle(get<1>(nearestPoint))) - m.point(vv)).norm();
					U.push(triple(dist + get<0>(nearestPoint), vv.idx(), get<1>(nearestPoint)));
				}
			}
		}
		U.pop();  
	} while (get<1>(nearestPoint) != b);

	dist = get<0>(nearestPoint);
	return FootPrint;
}

vector<pair<int, int>> ShortestPath_MinimalTree::findMinimalTree(Mesh& m)
{	
	typedef pair<int, int> pair;
	vector<pair> connect;
	if (IsChoosing)
	{
		//get input
		int meshNum = m.n_vertices();
		cout << "Please input id numbers of the vertices that spanning the tree(" << 0 << " to " << meshNum - 1 << "):";
		int inputnum;
		input.clear();
		while (cin >> inputnum)
		{
			input.push_back(inputnum);
			if (cin.get() == '\n')
			{
				break;
			}
		}		
	}
	//compute minimal tree	
	size_t n = input.size();
	vector<vector<float>> completeG_weight(n, vector<float>(n));
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			if (i != j)
			{
				float dist = 0.0f;
				vector<int> tmp = findPath(m, input[i], input[j], dist);
				completeG_weight[i][j] = completeG_weight[j][i] = dist;
			}
		}

	}
	//prim algorithm
	typedef tuple<float, int, int> triple;

	vector<int> setA;
	priority_queue<triple, vector<triple>, greater<triple>> setB; //distance,id,connect point id
	vector<int> setB_idx = input; //record id number in setB


	setA.push_back(input[n - 1]);
	setB_idx.pop_back();
	for (size_t i = 0; i < n - 1; i++)
	{
		float dist = completeG_weight[i][n - 1];
		triple tmp(dist, input[i], input[n - 1]);
		setB.push(tmp);
	}
	triple shortest;
	while (!setB.empty())
	{
		shortest = setB.top();
		if (find(setA.begin(), setA.end(), get<1>(shortest)) == setA.end())
		{
			setB.pop();
			setA.push_back(get<1>(shortest));
			connect.push_back(pair(get<1>(shortest), get<2>(shortest)));
			auto ite = find(setB_idx.begin(), setB_idx.end(), get<1>(shortest));
			setB_idx.erase(ite);
			//update set B
			int id1 = distance(input.begin(), find(input.begin(), input.end(), get<1>(shortest)));
			for (size_t i = 0; i < setB_idx.size(); i++)
			{
				int id2 = distance(input.begin(), find(input.begin(), input.end(), setB_idx[i]));

				float dist = completeG_weight[id1][id2];
				triple tmp(dist, setB_idx[i], get<1>(shortest));
				setB.push(tmp);
			}
		}
		else
		{
			setB.pop();
		}
		
	}
	IsChoosing = false;
	return connect;
	
}