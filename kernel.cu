#include <iostream>
#include <fstream>
#include "SVG.cuh"
#include "book.cuh"
#include <ctime>

using namespace std;

/*#define TEST_SVG*/
/*#define TEST_ICHHOST*/
#define TEST_SVG_SSSD2
/*#define TEST_ICH*/

#ifdef TEST_SVG
int main(int argc, char **argv)
{
	if (argc < 2)
	{
		cout << "[.exe] [in.obj]" << endl;
		return -1;
	}
	HANDLE_ERROR(cudaSetDevice(0));
	Mesh mesh, d_mesh;
	mesh.LoadFromFile(argv[1]);
	cout << "Mesh is loaded." << endl;
	mesh.copyToGPU(&d_mesh);

	cout << "Mesh is copied into GPU." << endl;

	SVG svg;
	svg.AssignMesh(&mesh, &d_mesh);
	svg.SetParameters(500, 10, 2);
	svg.Allocation();

	cout << "(Please check initial memory.)" << endl;
	system("pause");

	svg.ConstructSVG();

	cout << "SVG is built." << endl;
	cout << "(Please check final memory.)" << endl;
	system("pause");

	svg.Free();

	cout << "(Please check SVG-structure + Mesh memory.)" << endl;
	system("pause");

	svg.CopySVGToHost();
	cout << "SVG-structure is copied to host." << endl;
	system("pause");

	SVG::GraphDistInfo *graphDistInfos = new SVG::GraphDistInfo[mesh.vertNum];
	PriorityQueuesWithHandle<int> pq(mesh.vertNum);
	PriorityQueuesWithHandle<int>::PQItem *pqBuf = new PriorityQueuesWithHandle<int>::PQItem[mesh.vertNum + 1];
	pq.AssignMemory(pqBuf, mesh.vertNum);

	int s = -1, t = -1;

	while (true)
	{
		cout << "Input src & dst (input -1 -1 to exit): ";
		cin >> s >> t;
		if (s == -1 && t == -1) break;
		cout << "Solving ... " << endl;
		clock_t start = clock();
		svg.SolveSSSD(s, t, mesh, graphDistInfos, pq);
		clock_t end = clock();
		cout << "Solved." << endl;
		cout << "Time elapsed: " << (double)(end - start) / (double)CLOCKS_PER_SEC << endl;
		cout << "Distance: " << graphDistInfos[t].dist << endl;
		cout << "Passed vertices: ";
		while (t != s)
		{
			cout << t << " ";
			t = graphDistInfos[t].pathParentIndex;
		}
		cout << s << endl;

		// clear graphDistInfos
		for (int i = 0; i < mesh.vertNum; ++i)
		{
			graphDistInfos[i].dist = DBL_MAX;
			graphDistInfos[i].indexInPQ = -1;
			graphDistInfos[i].pathParentIndex = -1;
		}
		// clear pq
		pq.clear();
	}


	delete[] graphDistInfos;
	delete[] pqBuf;

	svg.FreeSVGStructure();
	mesh.clear();
	d_mesh.clearGPU();

	HANDLE_ERROR(cudaDeviceReset());
	return 0;
}
#endif

#ifdef TEST_SVG_SSSD2
int main(int argc, char **argv)
{
	if (argc < 2)
	{
		cout << "[.exe] [in.obj]" << endl;
		return -1;
	}
	HANDLE_ERROR(cudaSetDevice(0));
	Mesh mesh, d_mesh;
	mesh.LoadFromFile(argv[1]);
	cout << "Mesh is loaded." << endl;
	mesh.copyToGPU(&d_mesh);

	cout << "Mesh is copied into GPU." << endl;

	SVG svg;
	svg.AssignMesh(&mesh, &d_mesh);
	svg.SetParameters(500, 10, 2);
	svg.Allocation();

	cout << "(Please check initial memory.)" << endl;
	system("pause");

	svg.ConstructSVG();

	cout << "SVG is built." << endl;
	cout << "(Please check final memory.)" << endl;
	system("pause");

	svg.Free();

	cout << "(Please check SVG-structure + Mesh memory.)" << endl;
	system("pause");

	svg.CopySVGToHost();
	cout << "SVG-structure is copied to host." << endl;
	system("pause");

	SVG::GraphDistInfo *graphDistInfos = new SVG::GraphDistInfo[mesh.vertNum];
	PriorityQueuesWithHandle<int> pq(mesh.vertNum);
	PriorityQueuesWithHandle<int>::PQItem *pqBuf = new PriorityQueuesWithHandle<int>::PQItem[mesh.vertNum + 1];
	pq.AssignMemory(pqBuf, mesh.vertNum);

	ICH::SplitItem *splitInfos = new ICH::SplitItem[500 * 10 + 1];
	ICH::VertItem *vertInfos = new ICH::VertItem[500 * 2 + 1];
	SVG::PQWinItem *winPQ = new SVG::PQWinItem[WINPQ_SIZE];
	SVG::PQPseudoWinItem *pseudoWinPQs = new SVG::PQPseudoWinItem[PSEUDOWINPQ_SIZE];
	ICH::Window *storedWindows = new ICH::Window[STORED_WIN_BUF_SIZE];
	unsigned int *keptFaces = new unsigned int[KEPT_FACE_SIZE];

	int srcFace = -1, dstFace = -1;
	Vector3D srcPoint, dstPoint;

	while (true)
	{
		cout << "Input srcFace & dstFace (input -1 -1 to exit): ";
		cin >> srcFace >> dstFace;
		if (srcFace == -1 && dstFace == -1) break;
		cout << "Input srcPoint: ";
		cin >> srcPoint.x >> srcPoint.y >> srcPoint.z;
		cout << "Input dstPoint: ";
		cin >> dstPoint.x >> dstPoint.y >> dstPoint.z;
		// 		srcFace = 2382; dstFace = 1847;
		// 		srcPoint = Vector3D(0.077627, 74.551063, 22.270255);
		// 		dstPoint = Vector3D(-17.694569, 7.870463, 12.209889);
		cout << "Solving ... " << endl;
		SVG::SVGNode res;
		int lastVert = -1;
		clock_t start = clock();
		svg.SolveSSSD(srcFace, srcPoint, dstFace, dstPoint, mesh, 
			splitInfos, 500 * 10, 
			vertInfos, 500 * 2, 
			winPQ, pseudoWinPQs, storedWindows, keptFaces, graphDistInfos, pq, &res, &lastVert);
		clock_t end = clock();
		cout << "Solved." << endl;
		cout << "Time elapsed: " << (double)(end - start) / (double)CLOCKS_PER_SEC << endl;
		cout << "Distance: " << res.geodDist << endl;
		cout << "Next to src edge: " << res.nextToSrcEdge << endl;
		cout << "Next to src pos: " << res.nextToSrcX << endl;
		cout << "Next to dst edge: " << res.nextToDstEdge << endl;
		cout << "Next to dst pos: " << res.nextToDstX << endl;

		cout << "Passed vertices: ";
		while (lastVert != -1)
		{
			cout << lastVert << " ";
			lastVert = graphDistInfos[lastVert].pathParentIndex;
		}
		cout << endl;

		// clear graphDistInfos
		for (int i = 0; i < mesh.vertNum; ++i)
		{
			graphDistInfos[i].dist = DBL_MAX;
			graphDistInfos[i].indexInPQ = -1;
			graphDistInfos[i].pathParentIndex = -1;
		}
		// clear pq
		pq.clear();
	}


	delete[] graphDistInfos;
	delete[] pqBuf;

	svg.FreeSVGStructure();
	mesh.clear();
	d_mesh.clearGPU();

	HANDLE_ERROR(cudaDeviceReset());
	return 0;
}
#endif

#ifdef TEST_ICH
int main(int argc, char **argv)
{
	if (argc < 2)
	{
		cout << "[.exe] [in.obj]" << endl;
		return -1;
	}
	HANDLE_ERROR(cudaSetDevice(0));
	Mesh mesh;
	mesh.LoadFromFile(argv[1]);
	cout << "Mesh is loaded." << endl;

	unsigned K = 5000;
	unsigned splitInfoCoef = 10, vertInfoCoef = 2;
	ICH::SplitItem *splitInfoBuf = new ICH::SplitItem[K * splitInfoCoef + 1];
	ICH::VertItem *vertInfoBuf = new ICH::VertItem[K * vertInfoCoef + 1];

	PriorityQueues<ICH::Window>::PQItem *winPQs = new PriorityQueues<ICH::Window>::PQItem[WINPQ_SIZE];
	PriorityQueues<ICH::PseudoWindow>::PQItem *pseudoWinPQs = new PriorityQueues<ICH::PseudoWindow>::PQItem[PSEUDOWINPQ_SIZE];

	PriorityQueues<ICH::Window> winPQ(WINPQ_SIZE - 1);
	PriorityQueues<ICH::PseudoWindow> pseudoWinPQ(PSEUDOWINPQ_SIZE - 1);
	winPQ.AssignMemory(winPQs, WINPQ_SIZE - 1);
	pseudoWinPQ.AssignMemory(pseudoWinPQs, PSEUDOWINPQ_SIZE - 1);

	ICH::Window *storedWindowsBuf = new ICH::Window[STORED_WIN_BUF_SIZE];
	unsigned *keptFacesBuf = new unsigned[KEPT_FACE_SIZE];
	
	ICH ich;
	ich.AssignMesh(&mesh);
	ich.AssignBuffers(splitInfoBuf, K * splitInfoCoef, 
		vertInfoBuf, K * vertInfoCoef,
		winPQ, pseudoWinPQ,
		storedWindowsBuf, keptFacesBuf);

	ich.Clear();
	ich.AddSource(100);

	clock_t start = clock();
	ich.Execute(K);
	clock_t end = clock();
	cout << "Time: " << static_cast<double>(end - start) / static_cast<double>(CLOCKS_PER_SEC) << endl;

	ofstream output("result.obj");
	output << "mtllib texture.mtl" << endl;
	for (int i = 0; i < mesh.vertNum; ++i)
		output << "v " << mesh.verts[i].pos << endl;
	double maxDist = 0.0;
	for (int i = 0; i < mesh.vertNum; ++i)
		if (ich.GetDistanceTo(i) != DBL_MAX) maxDist = max(maxDist, ich.GetDistanceTo(i));
	for (int i = 0; i < mesh.vertNum; ++i)
	{
		double dist = ich.GetDistanceTo(i);
		if (dist == DBL_MAX) dist = 0.0;
		output << "vt " << dist / maxDist << " " << dist / maxDist << endl;
	}
	for (int i = 0; i < mesh.faceNum; ++i)
		output << "f" << " " << mesh.faces[i].verts[0] + 1 << "/" << mesh.faces[i].verts[0] + 1
		<< " " << mesh.faces[i].verts[1] + 1 << "/" << mesh.faces[i].verts[1] + 1
		<< " " << mesh.faces[i].verts[2] + 1 << "/" << mesh.faces[i].verts[2] + 1 << endl;
	return 0;
}
#endif