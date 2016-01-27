#include "SVG.cuh"
#include "book.cuh"
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

__global__ void constructSVG(Mesh mesh, int K, 
	SVG::PQWinItem *d_winPQs, SVG::PQPseudoWinItem *d_pseudoWinPQs,
	ICH::SplitInfo *d_splitInfoBuf, ICH::VertInfo *d_vertInfoBuf, 
	ICH::Window *d_storedWindowsBuf, unsigned *d_keptFacesBuf)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int totalThreadNum = blockDim.x * gridDim.x;

	PriorityQueues<ICH::Window> winPQ(WINPQ_SIZE-1);
	PriorityQueues<ICH::PseudoWindow> pseudoWinPQ(PSEUDOWINPQ_SIZE-1);
	winPQ.AssignMemory(d_winPQs + idx * WINPQ_SIZE); 
	pseudoWinPQ.AssignMemory(d_pseudoWinPQs + idx * PSEUDOWINPQ_SIZE);

	ICH ich;
	ich.AssignMesh(&mesh); 
	ich.AssignBuffers(d_splitInfoBuf + idx * mesh.edgeNum, 
		d_vertInfoBuf + idx * mesh.vertNum, 
		winPQ, pseudoWinPQ, 
		d_storedWindowsBuf + idx * STORED_WIN_BUF_SIZE, d_keptFacesBuf + idx * KEPT_FACE_SIZE);

	InitialValueGeodesic initGeodesic;
	initGeodesic.AssignMesh(&mesh);

	for (int i = idx; i < mesh.vertNum; i += totalThreadNum)
	{
		// TODO: run ICH
		ich.Clear();
		ich.AddSource(i);
 		ich.Execute(K);

		if (ich.GetDistanceTo(100) == DBL_MAX)
		{
			continue;
		}

		unsigned srcId;
		unsigned nextToSrcEdge, nextToDstEdge;
		double nextSrcX, nextDstX;
		ich.BuildGeodesicPathTo(100, srcId,
			nextToSrcEdge, nextSrcX, nextToDstEdge, nextDstX);

		// organize nextToSrcEdge & nextSrcX; and if they don't exist, it means src&dst are on the same face
		if (nextToSrcEdge == -1)
		{
			nextToSrcEdge = nextToDstEdge;
			nextSrcX = nextDstX;
		}

		InitialValueGeodesic::GeodesicKeyPoint dstPoint;
		if (nextToSrcEdge != -1 && !ich.pathPassVert)
		{
			initGeodesic.AssignLength(ich.GetDistanceTo(100));
			initGeodesic.AssignStartPoint(i);
			initGeodesic.AssignFirstKeyPoint(
				mesh.edges[nextToSrcEdge].twinEdge,
				mesh.edges[nextToSrcEdge].edgeLen - nextSrcX);
			dstPoint = initGeodesic.BuildGeodesicPath();
		}
		else
		{
			dstPoint.isInterior = true;
			dstPoint.faceIndex = -1;
			dstPoint.facePos3D = Vector3D();
		}
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

void SVG::SetParameters(int K_)
{
	K = K_;
}

bool SVG::Allocation()
{
	// allocation memories for PriorityQueues

	int totalThreadNum = THREAD_NUM * BLOCK_NUM; 

	HANDLE_ERROR(cudaMalloc((void**)&d_winPQs, totalThreadNum * WINPQ_SIZE * sizeof(PQWinItem)));
	HANDLE_ERROR(cudaMalloc((void**)&d_pseudoWinPQs, totalThreadNum * PSEUDOWINPQ_SIZE * sizeof(PQPseudoWinItem)));

	// allocation info buffers for ICH
	HANDLE_ERROR(cudaMalloc((void**)&d_splitInfoBuf, totalThreadNum * mesh->edgeNum * sizeof(ICH::SplitInfo)));
	HANDLE_ERROR(cudaMalloc((void**)&d_vertInfoBuf, totalThreadNum * mesh->vertNum * sizeof(ICH::VertInfo)));
	HANDLE_ERROR(cudaMalloc((void**)&d_storedWindowsBuf, totalThreadNum * STORED_WIN_BUF_SIZE * sizeof(ICH::Window)));
	HANDLE_ERROR(cudaMalloc((void**)&d_keptFacesBuf, totalThreadNum * KEPT_FACE_SIZE * sizeof(unsigned)));

	// allocation for SVG structure
	HANDLE_ERROR(cudaMalloc((void**)&d_svg, mesh->vertNum * K * sizeof(SVGNode)));
	return true;
}

void SVG::ConstructSVG()
{
	clock_t start = clock();
	constructSVG <<<BLOCK_NUM, THREAD_NUM >>>(*d_mesh, K, 
		d_winPQs, d_pseudoWinPQs, d_splitInfoBuf, d_vertInfoBuf, 
		d_storedWindowsBuf, d_keptFacesBuf);
	// TODO: organize the constructed SVG
	HANDLE_ERROR(cudaGetLastError());
	HANDLE_ERROR(cudaDeviceSynchronize());
	clock_t end = clock();
	cout << "Time consumed: " << (double)(end - start) / (double)CLOCKS_PER_SEC << endl;
	system("pause");
}