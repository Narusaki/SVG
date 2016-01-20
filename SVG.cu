#include "SVG.cuh"
#include "book.cuh"
#include <iostream>
#include <fstream>

using namespace std;

__global__ void constructSVG(Mesh mesh,
	SVG::PQWinItem *d_winPQs, SVG::PQPseudoWinItem *d_pseudoWinPQs,
	ICH::SplitInfo *d_splitInfoBuf, ICH::VertInfo *d_vertInfoBuf)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int totalThreadNum = blockDim.x * gridDim.x;
	int vertPerThread = (mesh.vertNum + totalThreadNum - 1) / totalThreadNum;

	PriorityQueues<ICH::Window> winPQ;
	PriorityQueues<ICH::PseudoWindow> pseudoWinPQ;
	winPQ.AssignMemory(d_winPQs + idx * WIN_PQ_SIZE); 
	pseudoWinPQ.AssignMemory(d_pseudoWinPQs + idx * PSEUDOWIN_PQ_SIZE);

	ICH ich;
	ich.AssignMesh(&mesh); 
	ich.AssignBuffers(d_splitInfoBuf + idx * mesh.edgeNum, 
		d_vertInfoBuf + idx * mesh.vertNum, 
		winPQ, pseudoWinPQ);

	for (int i = idx; i < mesh.vertNum; i += vertPerThread)
	{
		// TODO: run ICH
		ich.Clear();
		ich.AddSource(i);
		ich.Execute();
		break;
	}
}

SVG::SVG()
{

}

SVG::~SVG()
{

}

void SVG::AssignMesh(Mesh *mesh_, Mesh *d_mesh_)
{
	mesh = mesh_; d_mesh = d_mesh_;
}

bool SVG::Allocation()
{
	// allocation memories for PriorityQueues

	int totalThreadNum = THREAD_NUM * BLOCK_NUM; 

	cudaError_t cudaStatus;

	HANDLE_ERROR(cudaMalloc((void**)&d_winPQs, totalThreadNum * WIN_PQ_SIZE * sizeof(PQWinItem)));
	HANDLE_ERROR(cudaMalloc((void**)&d_pseudoWinPQs, totalThreadNum * PSEUDOWIN_PQ_SIZE * sizeof(PQPseudoWinItem)));

	// allocation info buffers for ICH
	HANDLE_ERROR(cudaMalloc((void**)&d_splitInfoBuf, totalThreadNum * mesh->edgeNum * sizeof(ICH::SplitInfo)));
	HANDLE_ERROR(cudaMalloc((void**)&d_vertInfoBuf, totalThreadNum * mesh->vertNum * sizeof(ICH::VertInfo)));

	return true;
}

void SVG::ConstructSVG()
{
	constructSVG <<<BLOCK_NUM, THREAD_NUM >>>(*d_mesh, d_winPQs, d_pseudoWinPQs, d_splitInfoBuf, d_vertInfoBuf);
	// TODO: organize the constructed SVG
}