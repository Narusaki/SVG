#ifndef SVG_CUH
#define SVG_CUH

#include "Mesh.cuh"
#include "PriorityQueue.cuh"
#include "ICH.cuh"
#include "InitialValueGeodesic.cuh"

#define BLOCK_NUM 1
#define THREAD_NUM 1

#define WINPQ_SIZE 2048
#define PSEUDOWINPQ_SIZE 2048
#define STORED_WIN_BUF_SIZE 100
#define KEPT_FACE_SIZE 10

class SVG
{
public:

	typedef PriorityQueues<ICH::Window>::PQItem PQWinItem;
	typedef PriorityQueues<ICH::PseudoWindow>::PQItem PQPseudoWinItem;

	SVG();
	~SVG();

	void AssignMesh(Mesh *mesh_, Mesh *d_mesh_);
	bool Allocation();

	void ConstructSVG();

private:
	Mesh *mesh, *d_mesh;

	// priority queues for Windows
	PQWinItem *d_winPQs;
	// priority queues for Pseudo Windows
	PQPseudoWinItem *d_pseudoWinPQs;

	// split info buffer for SplitInfo
	ICH::SplitInfo *d_splitInfoBuf;
	ICH::VertInfo *d_vertInfoBuf;

	// other buffers for ICH
	ICH::Window *d_storedWindowsBuf;
	unsigned *d_keptFacesBuf;
};

__global__ void constructSVG(Mesh mesh, 
	SVG::PQWinItem *d_winPQs, SVG::PQPseudoWinItem *d_pseudoWinPQs, 
	ICH::SplitInfo *d_splitInfoBuf, ICH::VertInfo *d_vertInfoBuf, 
	ICH::Window *d_storedWindowsBuf, unsigned *d_keptFacesBuf, 
	InitialValueGeodesic::GeodesicKeyPoint *d_dstPoints);
#endif