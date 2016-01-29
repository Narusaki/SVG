#include <iostream>
#include "SVG.cuh"
#include "book.cuh"
#include <ctime>

using namespace std;

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
		svg.SolveSSSD(s, t, &mesh, graphDistInfos, pq);
		clock_t end = clock();
		cout << "Solved." << endl;
		cout << "Time elapsed: " << (double)(end - start) / (double)CLOCKS_PER_SEC << endl;
		cout << "Distance: " << graphDistInfos[t].dist << endl;

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