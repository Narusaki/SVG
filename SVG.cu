#include "SVG.cuh"
#include "book.cuh"
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

__global__ void constructSVG(Mesh mesh,
	SVG::PQWinItem *d_winPQs, SVG::PQPseudoWinItem *d_pseudoWinPQs,
	ICH::SplitInfo *d_splitInfoBuf, ICH::VertInfo *d_vertInfoBuf, 
	ICH::Window *d_storedWindowsBuf, unsigned *d_keptFacesBuf,
	InitialValueGeodesic::GeodesicKeyPoint *d_dstPoints)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int totalThreadNum = blockDim.x * gridDim.x;

	PriorityQueues<ICH::Window> winPQ(WINPQ_SIZE);
	PriorityQueues<ICH::PseudoWindow> pseudoWinPQ(PSEUDOWINPQ_SIZE);
	winPQ.AssignMemory(d_winPQs + idx * mesh.vertNum); 
	pseudoWinPQ.AssignMemory(d_pseudoWinPQs + idx * mesh.vertNum);

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
		ich.Execute();
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
		d_dstPoints[i] = dstPoint;
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

	HANDLE_ERROR(cudaMalloc((void**)&d_winPQs, totalThreadNum * WINPQ_SIZE * sizeof(PQWinItem)));
	HANDLE_ERROR(cudaMalloc((void**)&d_pseudoWinPQs, totalThreadNum * PSEUDOWINPQ_SIZE * sizeof(PQPseudoWinItem)));

	// allocation info buffers for ICH
	HANDLE_ERROR(cudaMalloc((void**)&d_splitInfoBuf, totalThreadNum * mesh->edgeNum * sizeof(ICH::SplitInfo)));
	HANDLE_ERROR(cudaMalloc((void**)&d_vertInfoBuf, totalThreadNum * mesh->vertNum * sizeof(ICH::VertInfo)));
	HANDLE_ERROR(cudaMalloc((void**)&d_storedWindowsBuf, totalThreadNum * STORED_WIN_BUF_SIZE * sizeof(ICH::Window)));
	HANDLE_ERROR(cudaMalloc((void**)&d_keptFacesBuf, totalThreadNum * KEPT_FACE_SIZE * sizeof(unsigned)));

	return true;
}

void SVG::ConstructSVG()
{
	InitialValueGeodesic::GeodesicKeyPoint *dstPoints = new InitialValueGeodesic::GeodesicKeyPoint[mesh->vertNum];
	InitialValueGeodesic::GeodesicKeyPoint *d_dstPoints;
	HANDLE_ERROR(cudaMalloc((void**)&d_dstPoints, mesh->vertNum * sizeof(InitialValueGeodesic::GeodesicKeyPoint)));
// 	dstPoints[0].facePos3D.x = 1.414; dstPoints[0].facePos3D.y = -1.414; dstPoints[0].facePos3D.z = 1.414;
// 	HANDLE_ERROR(cudaMemcpy(d_dstPoints, dstPoints, mesh->vertNum * sizeof(InitialValueGeodesic::GeodesicKeyPoint), cudaMemcpyHostToDevice));

	clock_t start = clock();
	constructSVG <<<BLOCK_NUM, THREAD_NUM >>>(*d_mesh, 
		d_winPQs, d_pseudoWinPQs, d_splitInfoBuf, d_vertInfoBuf, 
		d_storedWindowsBuf, d_keptFacesBuf, 
		d_dstPoints);
	// TODO: organize the constructed SVG
	HANDLE_ERROR(cudaGetLastError());
	HANDLE_ERROR(cudaDeviceSynchronize());
	clock_t end = clock();
	cout << "Time consumed: " << (double)(end - start) / (double)CLOCKS_PER_SEC << endl;
	system("pause");
	HANDLE_ERROR(cudaMemcpy(dstPoints, d_dstPoints, mesh->vertNum * sizeof(InitialValueGeodesic::GeodesicKeyPoint), cudaMemcpyDeviceToHost));

 	ofstream output("outputDstPoints.obj");
	for (int i = 0; i < mesh->vertNum; ++i)
	{
		if (dstPoints[i].isInterior)
			output << "v " << dstPoints[i].facePos3D << endl;
		else
		{
			Vector3D p0 = mesh->verts[mesh->edges[dstPoints[i].edgeIndex].verts[0]].pos;
			Vector3D p1 = mesh->verts[mesh->edges[dstPoints[i].edgeIndex].verts[1]].pos;
			Vector3D univ = p1 - p0; univ.normalize();

			output << "v " << p0 + dstPoints[i].pos * univ << endl;
		}
	}

	delete[] dstPoints;
	HANDLE_ERROR(cudaFree(d_dstPoints));
}