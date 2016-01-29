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
	winPQ.AssignMemory(d_winPQs + idx * WINPQ_SIZE, WINPQ_SIZE); 
	pseudoWinPQ.AssignMemory(d_pseudoWinPQs + idx * PSEUDOWINPQ_SIZE, PSEUDOWINPQ_SIZE);

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
		SVG::SVGNode *d_cur_svg = d_svg + K * i;

		for (int j = 0; j < mesh.vertNum; ++j)
		{
			if (ich.GetDistanceTo(j) == DBL_MAX) continue;
			if (j == i) continue;

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
			d_cur_svg[d_svg_tails[i]].adjNode = j;
			d_cur_svg[d_svg_tails[i]].geodDist = ich.GetDistanceTo(j);
			d_cur_svg[d_svg_tails[i]].nextToSrcX = nextSrcX;
			d_cur_svg[d_svg_tails[i]].nextToSrcEdge = nextToSrcEdge;
			d_cur_svg[d_svg_tails[i]].nextToDstX = nextDstX;
			d_cur_svg[d_svg_tails[i]].nextToDstEdge = nextToDstEdge;

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
	d_svg = NULL; svg = NULL;
	d_svg_tails = NULL; svg_tails = NULL;
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

void SVG::Free()
{
	HANDLE_ERROR(cudaFree(d_winPQs));
	HANDLE_ERROR(cudaFree(d_pseudoWinPQs));
	HANDLE_ERROR(cudaFree(d_splitInfoBuf));
	HANDLE_ERROR(cudaFree(d_vertInfoBuf));
	HANDLE_ERROR(cudaFree(d_storedWindowsBuf));
	HANDLE_ERROR(cudaFree(d_keptFacesBuf));
}

void SVG::FreeSVGStructure()
{
	HANDLE_ERROR(cudaFree(d_svg));
	HANDLE_ERROR(cudaFree(d_svg_tails));
	if (svg) delete[] svg;
	if (svg_tails) delete[] svg_tails;
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

	int *svg_tails = new int[mesh->vertNum];
	HANDLE_ERROR(cudaMemcpy(svg_tails, d_svg_tails, mesh->vertNum * sizeof(int), cudaMemcpyDeviceToHost));
	double degree = 0;
	for (int i = 0; i < mesh->vertNum; ++i)
		degree += svg_tails[i];
	degree /= mesh->vertNum;
	cout << "Average degree of node in SVG: " << degree << endl;
}

void SVG::CopySVGToHost()
{
	svg = new SVGNode[mesh->vertNum * K];
	svg_tails = new int[mesh->vertNum];

	HANDLE_ERROR(cudaMemcpy(svg, d_svg, mesh->vertNum * K * sizeof(SVGNode), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(svg_tails, d_svg_tails, mesh->vertNum * sizeof(int), cudaMemcpyDeviceToHost));
// 	cout << mesh->vertNum << endl;
// 	for (int j = 0; j < mesh->vertNum; ++j)
// 	{
// 		cout << "Node " << j << " degree: " << svg_tails[j] << endl;
// 		SVGNode *curSVGList = svg + K * j;
// 		for (int i = 0; i < svg_tails[j]; ++i)
// 		{
// 			printf("%d %f %d %f %d %f\n", curSVGList[i].adjNode, curSVGList[i].geodDist, 
// 				curSVGList[i].nextToSrcEdge, curSVGList[i].nextToSrcX, 
// 				curSVGList[i].nextToDstEdge, curSVGList[i].nextToDstX);
// 			system("pause");
// 		}
// 	}
}

__host__ __device__ void SVG::SolveSSSD(int s, int t, Mesh *mesh, GraphDistInfo * graphDistInfos, PriorityQueuesWithHandle<int> pq)
{
	// TODO: use Astar algorithm to search dist&path to t; use Euclidean dist as heuristic prediction
	// TODO: initialize graphDistInfos & pq
	searchType = ASTAR;
	graphDistInfos[s].dist = 0.0;
	graphDistInfos[s].pathParentIndex = -1;
	pq.push(s, &graphDistInfos[s].indexInPQ, 0.0 + (mesh->verts[s].pos - mesh->verts[t].pos).length());
	Astar(mesh, t, graphDistInfos, pq);
}

__host__ __device__ void SVG::SolveMSMD(int *sources, int Ns, int *destinations, int Nd, Mesh *mesh, GraphDistInfo * graphDistInfos, PriorityQueuesWithHandle<int> pq)
{
	// TODO: use Dijkstra algorithm to search min-dist&path to destinations
	// TODO: initialize graphDistInfos & pq
	searchType = DIJKSTRA;
	Astar(mesh, -1, graphDistInfos, pq);
}

__host__ __device__ void SVG::Astar(Mesh *mesh, int t, GraphDistInfo * graphDistInfos, PriorityQueuesWithHandle<int> pq)
{
	SVGNode *local_svg = NULL;
	int *local_svg_tails = NULL;

#ifdef __CUDA_ARCH__
	local_svg = d_svg; local_svg_tails = d_svg_tails;
#else
	local_svg = svg; local_svg_tails = svg_tails;
#endif

	while (!pq.empty())
	{
		int curNodeIndex = pq.pop();
		if (searchType == ASTAR && curNodeIndex == t)
		{
			pq.clear();
			break;
		}
		SVGNode *curNodeList = local_svg + K * curNodeIndex;
		for (int i = 0; i < local_svg_tails[curNodeIndex]; ++i)
		{
			int adjNodeIndex = curNodeList[i].adjNode;
			double newDist = graphDistInfos[curNodeIndex].dist + curNodeList[i].geodDist;
			if (newDist >= graphDistInfos[adjNodeIndex].dist) continue;

			graphDistInfos[adjNodeIndex].dist = newDist;
			graphDistInfos[adjNodeIndex].pathParentIndex = curNodeIndex;

			double priority = 0.0;
			switch (searchType)
			{
			case SVG::ASTAR: priority = newDist + (mesh->verts[adjNodeIndex].pos - mesh->verts[t].pos).length(); break;
			case SVG::DIJKSTRA: priority = newDist; break;
			default: priority = newDist; break;
			}

			if (graphDistInfos[adjNodeIndex].indexInPQ == -1)
				pq.push(adjNodeIndex, &graphDistInfos[adjNodeIndex].indexInPQ, priority);
			else
				pq.decrease(graphDistInfos[adjNodeIndex].indexInPQ, priority);
		}
	}
}