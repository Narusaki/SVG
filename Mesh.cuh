#ifndef MESH_CUH
#define MESH_CUH

#include "Geometry.cuh"

class Edge	//
{
public:
	int verts[2];
	int twinEdge, nextEdge, prevEdge;
	int faceId;
	double edgeLen;
};

class Vertex
{
public:
	Vector3D pos;
	double angle;
	int firstEdge;

	Vertex() { angle = 0.0; firstEdge = -1; }
};

class Face
{
public:
	unsigned verts[3];
	unsigned edges[3];
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

	Vertex *verts;
	int vertNum;

	Face *faces;
	int faceNum;
};
#endif