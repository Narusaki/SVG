#ifndef SVG_CUH
#define SVG_CUH

#include "Mesh.cuh"
#include "PriorityQueue.cuh"
#include "ICH.cuh"

#define BLOCK_NUM 2
#define THREAD_NUM 256
#define WIN_PQ_SIZE 512
#define PSEUDOWIN_PQ_SIZE 512

class SVG
{
public:

	typedef PriorityQueues<ICH::Window>::PQItem PQWinItem;
	typedef PriorityQueues<ICH::PseudoWindow>::PQItem PQPseudoWinItem;

	SVG();
	~SVG();

	bool Allocation();

	void ConstructSVG();

private:
	Mesh *mesh;

	PQWinItem *d_winPQs;
	PQPseudoWinItem *d_pseudoWinPQs;

	ICH::SplitInfo *d_splitInfoBuf;
	ICH::VertInfo *d_vertInfoBuf;
};

__global__ void constructSVG(Mesh *mesh, 
	SVG::PQWinItem *d_winPQs, SVG::PQPseudoWinItem *d_pseudoWinPQs, 
	ICH::SplitInfo *d_splitInfoBuf, ICH::VertInfo *d_vertInfoBuf);
#endif