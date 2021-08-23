#pragma once
#include "MeshViewer\MeshDefinition.h"

class ShortestPath_MinimalTree 
{
public:
	ShortestPath_MinimalTree();
	~ShortestPath_MinimalTree();
	std::vector<int> findPath(Mesh& m, int a, int b, float& dist);
	std::vector<std::pair<int, int>> findMinimalTree(Mesh& m);

public:
	bool IsChoosing;

private:
	std::vector<int> input;

};

