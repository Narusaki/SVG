#ifndef SVG_CUH
#define SVG_CUH

#include "Mesh.cuh"
#include "PriorityQueue.cuh"
#include "PriorityQueueWithHandle.cuh"
#include "ICHDevice.cuh"
#include "InitialValueGeodesic.cuh"

const int SVG_BLOCK_NUM = 64;
const int SVG_THREAD_NUM = 128;

#define WINPQ_SIZE 1024
#define PSEUDOWINPQ_SIZE 1024
#define STORED_WIN_BUF_SIZE 100
#define KEPT_FACE_SIZE 10

class SVG
{
public:
	struct SVGNode	// 36 bytes
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
		int srcId;
		int indexInPQ;

		__host__ __device__ GraphDistInfo() { dist = DBL_MAX; pathParentIndex = -1; srcId = -1; indexInPQ = -1; }
	};

	enum SEARCHTYPE
	{
		ASTAR, DIJKSTRA,
	};

public:

	typedef PriorityQueues<ICHDevice::Window>::PQItem PQWinItem;
	typedef PriorityQueues<ICHDevice::PseudoWindow>::PQItem PQPseudoWinItem;

	__host__ __device__ SVG();
	__host__ __device__ ~SVG();

	void AssignMesh(Mesh *mesh_, Mesh *d_mesh_);
	// K represents the number of vertices in the geodesic-disk covered by each vertex
	// splitInfoCoef represents the (average) valance of each vertex (6 for instance)
	// vertInfoCoef has the same meaning as K (1 for instance)
	void SetParameters(int K_, unsigned splitInfoCoef_, unsigned vertInfoCoef_);
	bool Allocation();
	void Free();
	void FreeSVGStructure();

	void ConstructSVG();

	// must invoke this method before searching SVG on host
	void CopySVGToHost();

	void SaveSVGToFile(const char *fileName);
	void LoadSVGFromFile(const char *fileName);

	__host__ __device__ void SolveSSSD(int s, int t, Mesh mesh,
		GraphDistInfo * graphDistInfos, PriorityQueuesWithHandle<int> pq);

	__host__ __device__ void SolveSSSD(int f0, Vector3D p0, int f1, Vector3D p1, Mesh mesh,
		ICHDevice::SplitItem *d_splitInfos, unsigned splitInfoSize,
		ICHDevice::VertItem *d_vertInfos, unsigned vertInfoSize,
		PQWinItem *winPQBuf, PQPseudoWinItem *pseudoWinPQBuf,
		ICHDevice::Window *storedWindows, unsigned int *keptFaces,
		GraphDistInfo * graphDistInfos, PriorityQueuesWithHandle<int> pq,
		SVGNode *res, int *lastVert);

	__host__ __device__ void SolveMSAD(int *sources, double *sourceWeights, int Ns,
		Mesh mesh,
		GraphDistInfo * graphDistInfos, PriorityQueuesWithHandle<int> pq);

private:
	// Astar algorithm, degenerate to Dijkstra if choose the so-far distance only as priority
	__host__ __device__ void Astar(Mesh *mesh, int t, GraphDistInfo * graphDistInfos, PriorityQueuesWithHandle<int> pq);

public:
	Mesh *mesh, *d_mesh;

public:
	// priority queues for Windows
	PQWinItem *d_winPQs;
	// priority queues for Pseudo Windows
	PQPseudoWinItem *d_pseudoWinPQs;

	// split info buffer for SplitInfo
	ICHDevice::SplitItem *d_splitInfoBuf; unsigned splitInfoCoef;
	ICHDevice::VertItem *d_vertInfoBuf; unsigned vertInfoCoef;

	// other buffers for ICH
	ICHDevice::Window *d_storedWindowsBuf;
	unsigned *d_keptFacesBuf;

	// search parameter
	SEARCHTYPE searchType;

private:
	SVGNode *d_svg, *svg;
	int *d_svg_tails, *svg_tails;
public:
	int K;
};

__global__ void constructSVG(Mesh mesh, int K,
	SVG::PQWinItem *d_winPQs, SVG::PQPseudoWinItem *d_pseudoWinPQs,
	ICHDevice::SplitItem *d_splitInfoBuf, unsigned splitInfoCoef,
	ICHDevice::VertItem *d_vertInfoBuf, unsigned vertInfoCoef,
	ICHDevice::Window *d_storedWindowsBuf, unsigned *d_keptFacesBuf,
	SVG::SVGNode *d_svg, int *d_svg_tails);
#endif