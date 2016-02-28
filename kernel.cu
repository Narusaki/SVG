#include <iostream>
#include "SVG.cuh"
#include "book.cuh"
#include <ctime>

using namespace std;

#define TEST_SVG
/*#define TEST_ICHHOST*/
/*#define TEST_SVG_SSSD2*/
/*#define TEST_SVG_SSSD2_MULTIPLE_RUNNING*/

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
	svg.SetParameters(500);
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
	svg.SetParameters(500);
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

	ICH::SplitInfo *splitInfos = new ICH::SplitInfo[mesh.edgeNum];
	ICH::VertInfo *vertInfos = new ICH::VertInfo[mesh.vertNum];
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
		svg.SolveSSSD(srcFace, srcPoint, dstFace, dstPoint, mesh, splitInfos, vertInfos, winPQ, pseudoWinPQs, storedWindows, keptFaces, graphDistInfos, pq, &res, &lastVert);
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

#ifdef TEST_SVG_SSSD2_MULTIPLE_RUNNING
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
	svg.SetParameters(500);
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

	ICH::SplitInfo *splitInfos = new ICH::SplitInfo[mesh.edgeNum];
	ICH::VertInfo *vertInfos = new ICH::VertInfo[mesh.vertNum];
	SVG::PQWinItem *winPQ = new SVG::PQWinItem[WINPQ_SIZE];
	SVG::PQPseudoWinItem *pseudoWinPQs = new SVG::PQPseudoWinItem[PSEUDOWINPQ_SIZE];
	ICH::Window *storedWindows = new ICH::Window[STORED_WIN_BUF_SIZE];
	unsigned int *keptFaces = new unsigned int[KEPT_FACE_SIZE];

	ifstream input("r0r1r2.txt");
	string curLine;
	int r[3] = { -1, -1, -1 };
	int srcFace = -1, dstFace = -1;
	Vector3D srcPoint, dstPoint;
	while (getline(input, curLine))
	{
		stringstream sin;
		sin << curLine;
		if (r[0] == -1) {
			sin >> r[0]; continue;
		}
		else if (r[1] == -1) {
			sin >> r[1]; continue;
		}
		else sin >> r[2];
		cout << "Geodesic difference between " << r[0] << " and " << r[1] << " ..." << endl;
		ifstream inputIndividual0(("individual_" + to_string(r[0])));
		ifstream inputIndividual1(("individual_" + to_string(r[1])));

		string curLine1, curLine2;
		while (getline(inputIndividual0, curLine1))
		{
			getline(inputIndividual1, curLine2);
			sin.clear();
			sin << curLine1; sin >> srcFace >> srcPoint.x >> srcPoint.y >> srcPoint.z;
			sin.clear();
			sin << curLine2; sin >> dstFace >> dstPoint.x >> dstPoint.y >> dstPoint.z;

			cout << "Solving (" << srcFace << " " << srcPoint << ") and (" << dstFace << " " << dstPoint << ") ... " << endl;
			SVG::SVGNode res;
			int lastVert = -1;
			clock_t start = clock();
			svg.SolveSSSD(srcFace, srcPoint, dstFace, dstPoint, mesh, splitInfos, vertInfos, winPQ, pseudoWinPQs, storedWindows, keptFaces, graphDistInfos, pq, &res, &lastVert);
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

			// clear
			for (int i = 0; i < mesh.edgeNum; ++i)
			{
				splitInfos[i].dist = DBL_MAX;
				splitInfos[i].pseudoSrcId = -1;
				splitInfos[i].x = DBL_MAX;
			}
			for (int i = 0; i < mesh.vertNum; ++i)
			{
				vertInfos[i].birthTime = -1;
				vertInfos[i].dist = DBL_MAX;
				vertInfos[i].enterEdge = -1;
			}
			for (int i = 0; i < mesh.vertNum; ++i)
			{
				graphDistInfos[i].dist = DBL_MAX;
				graphDistInfos[i].indexInPQ = -1;
				graphDistInfos[i].pathParentIndex = -1;
			}
			pq.clear();
			/*system("pause");*/
		}
		r[0] = r[1] = r[2] = -1;
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