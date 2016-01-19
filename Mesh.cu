#include "Mesh.cuh"
#include "book.cuh"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <map>

using namespace std;

Mesh::Mesh()
{
	edges = NULL; angles = NULL;
	edgeNum = 0; vertNum = 0;
}

Mesh::~Mesh()
{

}

bool Mesh::LoadFromFile(const char *fileName)
{
	ifstream input(fileName);
	if (!input)
	{
		cout << "Cannot open file " << fileName << endl;
		return false;
	}

	string curLine;
	char t;
	Vector3D point;
	unsigned v0, v1, v2;
	stringstream sin;
	vector<Vector3D> vertsVec;
	vector<int> facesVec;
	vector<Edge> edgesVec;

	while (getline(input, curLine))
	{
		sin.clear();
		if (curLine[0] == 'v')
		{
			if (curLine[1] == ' ' || curLine[1] == '\t')
			{
				sin << curLine;
				sin >> t >> point.x >> point.y >> point.z;
				vertsVec.push_back(point);
			}
		}
		else if (curLine[0] == 'f')
		{
			sin << curLine;
			sin >> t >> v0 >> v1 >> v2;
			facesVec.push_back(v0 - 1); facesVec.push_back(v1 - 1); facesVec.push_back(v2 - 1);
		}
	}
	
	edgesVec.resize(facesVec.size());

	int edgeIdx = 0;
	map< pair<int, int>, int > verts2edge;

	for (int i = 0; i < facesVec.size(); i+=3)
	{
		for (int j = 0; j < 3; ++j)
		{
			int v0 = facesVec[i + j], v1 = facesVec[i + (j + 1) % 3];
			edgesVec[edgeIdx].verts[0] = v0; edgesVec[edgeIdx].verts[1] = v1;
			edgesVec[edgeIdx].edgeLen = (vertsVec[v0] - vertsVec[v1]).length();
			edgesVec[edgeIdx].nextEdge = i + (j + 1) % 3;
			edgesVec[edgeIdx].prevEdge = i + (j + 2) % 3;

			if (v0 > v1) std::swap(v0, v1);
			auto iter = verts2edge.find(make_pair(v0, v1));
			if (iter != verts2edge.end())
			{
				edgesVec[edgeIdx].twinEdge = iter->second;
				edgesVec[iter->second].twinEdge = edgeIdx;
			}
			else
			{
				verts2edge[make_pair(v0, v1)] = edgeIdx;
			}
			++edgeIdx;
		}
	}

	edgeNum = edgesVec.size();
	edges = new Edge[edgeNum];
	copy(edgesVec.begin(), edgesVec.end(), edges);

	vertNum = vertsVec.size();
	angles = new double[vertNum];
	memset(angles, 0, vertNum * sizeof(double));

	edgeAdjToVert = new int[vertNum];

	faceNum = facesVec.size() / 3;

	for (int i = 0; i < edgeNum; ++i)
	{
		double l0 = edges[i].edgeLen;
		double l1 = edges[edges[i].nextEdge].edgeLen;
		double l2 = edges[edges[i].prevEdge].edgeLen;

		double curAngle = (l0*l0 + l2*l2 - l1*l1) / (2.0 * l0*l2);
		if (curAngle > 1.0) curAngle = 1.0; else if (curAngle < -1.0) curAngle = -1.0;
		curAngle = acos(curAngle);
		angles[edges[i].verts[0]] += curAngle;

		curAngle = (l0*l0 + l1*l1 - l2*l2) / (2.0 * l0*l1);
		if (curAngle > 1.0) curAngle = 1.0; else if (curAngle < -1.0) curAngle = -1.0;
		curAngle = acos(curAngle);
		angles[edges[i].verts[1]] += curAngle;

		edgeAdjToVert[edges[i].verts[0]] = i;
	}

	return true;
}

bool Mesh::copyToGPU(Mesh *d_mesh)
{
	cudaError_t cudaStatus;
	d_mesh->edgeNum = edgeNum;
	d_mesh->vertNum = vertNum;
	d_mesh->faceNum = faceNum;
	HANDLE_ERROR(cudaMalloc((void**)&(d_mesh->edges), edgeNum * sizeof(Edge)));
	HANDLE_ERROR(cudaMalloc((void**)(&d_mesh->angles), vertNum * sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void**)&(d_mesh->edgeAdjToVert), vertNum * sizeof(int)));

	HANDLE_ERROR(cudaMemcpy(d_mesh->edges, edges, edgeNum * sizeof(Edge), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_mesh->angles, angles, vertNum * sizeof(double), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_mesh->edgeAdjToVert, edgeAdjToVert, vertNum * sizeof(int), cudaMemcpyHostToDevice));
	return true;
}

void Mesh::clear()
{
	if (edges) delete[] edges;
	if (angles) delete[] angles;
	if (edgeAdjToVert) delete[] edgeAdjToVert;

	edgeNum = 0; vertNum = 0; faceNum = 0;
}

void Mesh::clearGPU()
{
	HANDLE_ERROR(cudaFree(edges));
	HANDLE_ERROR(cudaFree(angles));
	HANDLE_ERROR(cudaFree(edgeAdjToVert));
	edgeNum = 0; vertNum = 0; faceNum = 0;
}