#ifndef SVG_CUH
#define SVG_CUH

#include "Mesh.cuh"
#include "PriorityQueue.cuh"
#include "PriorityQueueWithHandle.cuh"
#include "ICH.cuh"
#include "InitialValueGeodesic.cuh"

#define BLOCK_NUM 8
#define THREAD_NUM 128

#define WINPQ_SIZE 1024
#define PSEUDOWINPQ_SIZE 1024
#define STORED_WIN_BUF_SIZE 100
#define KEPT_FACE_SIZE 10

class SVG
{
public:
	struct SVGNode	// 28 bytes
	{
		int adjNode;
		double geodDist;
		unsigned nextToSrcEdge, nextToDstEdge;
		double nextToSrcX, nextToDstX;
	};

	struct GraphDistInfo
	{
		double dist;
		int pathParentIndex;
		int indexInPQ;

		__host__ __device__ GraphDistInfo() { dist = DBL_MAX; pathParentIndex = -1; indexInPQ = -1; }
	};

	enum SEARCHTYPE
	{
		ASTAR, DIJKSTRA, 
	};

public:

	typedef PriorityQueues<ICH::Window>::PQItem PQWinItem;
	typedef PriorityQueues<ICH::PseudoWindow>::PQItem PQPseudoWinItem;

	SVG();
	~SVG();

	void AssignMesh(Mesh *mesh_, Mesh *d_mesh_);
	void SetParameters(int K_);
	bool Allocation();
	void Free();
	void FreeSVGStructure();

	void ConstructSVG();

	// must invoke this method before searching SVG on host
	void CopySVGToHost();

	__host__ __device__ void SolveSSSD(int s, int t, Mesh *mesh, 
		GraphDistInfo * graphDistInfos, PriorityQueuesWithHandle<int> pq);

	__host__ __device__ void SolveMSMD(int *sources, int Ns, int *destinations, int Nd, Mesh *mesh,
		GraphDistInfo * graphDistInfos, PriorityQueuesWithHandle<int> pq);

private:
	// Astar algorithm, degenerate to Dijkstra if choose the so-far distance only as priority
	__host__ __device__ void Astar(Mesh *mesh, int t, GraphDistInfo * graphDistInfos, PriorityQueuesWithHandle<int> pq);

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

	// search parameter
	SEARCHTYPE searchType;

private:
	SVGNode *d_svg, *svg;
	int *d_svg_tails, *svg_tails;
	int K;
};

__global__ void constructSVG(Mesh mesh, int K, 
	SVG::PQWinItem *d_winPQs, SVG::PQPseudoWinItem *d_pseudoWinPQs, 
	ICH::SplitInfo *d_splitInfoBuf, ICH::VertInfo *d_vertInfoBuf, 
	ICH::Window *d_storedWindowsBuf, unsigned *d_keptFacesBuf, 
	SVG::SVGNode *d_svg, int *d_svg_tails);
#endif