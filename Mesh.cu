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
	edges = NULL; verts = NULL; faces = NULL;
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
	Vertex point; Face face;
	stringstream sin;
	vector<Vertex> vertsVec;
	vector<Face> facesVec;
	vector<Edge> edgesVec;

	while (getline(input, curLine))
	{
		sin.clear();
		if (curLine[0] == 'v')
		{
			if (curLine[1] == ' ' || curLine[1] == '\t')
			{
				sin << curLine;
				sin >> t >> point.pos.x >> point.pos.y >> point.pos.z;
				vertsVec.push_back(point);
			}
		}
		else if (curLine[0] == 'f')
		{
			sin << curLine;
			sin >> t >> face.verts[0] >> face.verts[1] >> face.verts[2];
			--face.verts[0]; --face.verts[1]; --face.verts[2];
			facesVec.push_back(face);
		}
	}

	edgesVec.resize(facesVec.size() * 3);

	int edgeIdx = 0;
	map< pair<int, int>, int > verts2edge;

	for (int i = 0; i < facesVec.size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			int v0 = facesVec[i].verts[j], v1 = facesVec[i].verts[(j + 1) % 3];
			edgesVec[edgeIdx].verts[0] = v0; edgesVec[edgeIdx].verts[1] = v1;
			edgesVec[edgeIdx].edgeLen = (vertsVec[v0].pos - vertsVec[v1].pos).length();
			edgesVec[edgeIdx].nextEdge = i * 3 + (j + 1) % 3;
			edgesVec[edgeIdx].prevEdge = i * 3 + (j + 2) % 3;
			edgesVec[edgeIdx].faceId = i;
			facesVec[i].edges[j] = edgeIdx;

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
	verts = new Vertex[vertNum];
	copy(vertsVec.begin(), vertsVec.end(), verts);

	faceNum = facesVec.size();
	faces = new Face[faceNum];
	copy(facesVec.begin(), facesVec.end(), faces);

	for (int i = 0; i < edgeNum; ++i)
	{
		double l0 = edges[i].edgeLen;
		double l1 = edges[edges[i].nextEdge].edgeLen;
		double l2 = edges[edges[i].prevEdge].edgeLen;

		double curAngle = (l0*l0 + l2*l2 - l1*l1) / (2.0 * l0*l2);
		if (curAngle > 1.0) curAngle = 1.0; else if (curAngle < -1.0) curAngle = -1.0;
		curAngle = acos(curAngle);
		verts[edges[i].verts[0]].angle += curAngle;
		/*
		curAngle = (l0*l0 + l1*l1 - l2*l2) / (2.0 * l0*l1);
		if (curAngle > 1.0) curAngle = 1.0; else if (curAngle < -1.0) curAngle = -1.0;
		curAngle = acos(curAngle);
		angles[edges[i].verts[1]] += curAngle;
		*/
		verts[edges[i].verts[0]].firstEdge = i;
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
	HANDLE_ERROR(cudaMalloc((void**)&(d_mesh->verts), vertNum * sizeof(Vertex)));
	HANDLE_ERROR(cudaMalloc((void**)&(d_mesh->faces), faceNum * sizeof(Face)));

	HANDLE_ERROR(cudaMemcpy(d_mesh->edges, edges, edgeNum * sizeof(Edge), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_mesh->verts, verts, vertNum * sizeof(Vertex), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_mesh->faces, faces, faceNum * sizeof(Face), cudaMemcpyHostToDevice));
	return true;
}

void Mesh::clear()
{
	if (edges) delete[] edges;
	if (verts) delete[] verts;
	if (faces) delete[] faces;

	edgeNum = 0; vertNum = 0; faceNum = 0;
}

void Mesh::clearGPU()
{
	HANDLE_ERROR(cudaFree(edges));
	HANDLE_ERROR(cudaFree(verts));
	HANDLE_ERROR(cudaFree(faces));
	edgeNum = 0; vertNum = 0; faceNum = 0;
}