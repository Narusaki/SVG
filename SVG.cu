#include "SVG.cuh"
#include "book.cuh"
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

__global__ void constructSVG(Mesh mesh, int K, 
	SVG::PQWinItem *d_winPQs, SVG::PQPseudoWinItem *d_pseudoWinPQs,
	ICH::SplitInfo *d_splitInfoBuf, ICH::VertInfo *d_vertInfoBuf, 
	ICH::Window *d_storedWindowsBuf, unsigned *d_keptFacesBuf, 
	SVG::SVGNode *d_svg, int *d_svg_tails)
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

		d_svg_tails[i] = 0;

		for (int j = 0; j < mesh.vertNum; ++j)
		{
			if (ich.GetDistanceTo(j) == DBL_MAX) continue;

			unsigned srcId;
			unsigned nextToSrcEdge, nextToDstEdge;
			double nextSrcX, nextDstX;
			ich.BuildGeodesicPathTo(j, srcId,
				nextToSrcEdge, nextSrcX, nextToDstEdge, nextDstX);

			// if path passes a (saddle) vertex, then continue
			if (ich.pathPassVert) continue;

			// organize nextToSrcEdge & nextSrcX; and if they don't exist, it means src&dst are the same
			// if nextToDstEdge & nextDstX does not exist either (this should not happen when constructing svg (vert-vert-pair))
			if (nextToSrcEdge == -1)
			{
				nextToSrcEdge = nextToDstEdge;
				nextSrcX = nextDstX;
			}
			d_svg[d_svg_tails[i]].adjNode = j;
			d_svg[d_svg_tails[i]].geodDist = ich.GetDistanceTo(j);
			d_svg[d_svg_tails[i]].nextToSrcX = nextSrcX;
			d_svg[d_svg_tails[i]].nextToSrcEdge = nextToSrcEdge;
			d_svg[d_svg_tails[i]].nextToDstX = nextDstX;
			d_svg[d_svg_tails[i]].nextToDstEdge = nextToDstEdge;

			++d_svg_tails[i];

			// in case of overflow
			if (d_svg_tails[i] == K) break;
		}

// 		InitialValueGeodesic::GeodesicKeyPoint dstPoint;
// 		if (nextToSrcEdge != -1 && !ich.pathPassVert)
// 		{
// 			initGeodesic.AssignLength(ich.GetDistanceTo(100));
// 			initGeodesic.AssignStartPoint(i);
// 			initGeodesic.AssignFirstKeyPoint(
// 				mesh.edges[nextToSrcEdge].twinEdge,
// 				mesh.edges[nextToSrcEdge].edgeLen - nextSrcX);
// 			dstPoint = initGeodesic.BuildGeodesicPath();
// 		}
// 		else
// 		{
// 			dstPoint.isInterior = true;
// 			dstPoint.faceIndex = -1;
// 			dstPoint.facePos3D = Vector3D();
// 		}
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
	HANDLE_ERROR(cudaMalloc((void**)&d_svg_tails, mesh->vertNum * sizeof(int)));

	return true;
}

void SVG::ConstructSVG()
{
	clock_t start = clock();
	constructSVG <<<BLOCK_NUM, THREAD_NUM >>>(*d_mesh, K, 
		d_winPQs, d_pseudoWinPQs, d_splitInfoBuf, d_vertInfoBuf, 
		d_storedWindowsBuf, d_keptFacesBuf, 
		d_svg, d_svg_tails);
	// TODO: organize the constructed SVG
	HANDLE_ERROR(cudaGetLastError());
	HANDLE_ERROR(cudaDeviceSynchronize());
	clock_t end = clock();
	cout << "Time consumed: " << (double)(end - start) / (double)CLOCKS_PER_SEC << endl;
	system("pause");
}