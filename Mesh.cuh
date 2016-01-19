#ifndef MESH_CUH
#define MESH_CUH

#include "Geometry.cuh"

class Edge
{
public:
	int verts[2];
	int twinEdge, nextEdge, prevEdge;
	double edgeLen;
};

class Mesh
{
public:
	Mesh();
	~Mesh();

	bool LoadFromFile(const char *fileName);
	bool copyToGPU(Mesh *d_mesh);		// note that d_mesh should still be a host variable
	void clear();
	void clearGPU();

public:
	Edge *edges;
	int edgeNum;
	double *angles;
	int vertNum;
	int *edgeAdjToVert;
	int faceNum;
};
#endif