// ICH.cpp : Defines the exported functions for the DLL application.
//

#include "ICH.cuh"
#include <iostream>
#include <list>

using namespace std;

// This is the constructor of a class that has been exported.
// see ICH.h for the class definition
__host__ __device__ ICH::ICH()
{
	sourceVert = -1; sourcePointFace = -1;

	mesh = NULL;
	storedWindows = NULL; keptFaces = NULL;
	storedWindowsIdx = 0; keptFacesIdx = 0;

	numOfWinGen = 0;
	maxWinQSize = 0;
	maxPseudoQSize = 0;
	totalCalcVertNum = 0;
	return;
}

__host__ __device__ ICH::~ICH()
{
	return;
}

__host__ __device__ void ICH::AssignMesh(Mesh *mesh_)
{
	mesh = mesh_;
}

__host__ __device__ void ICH::AssignBuffers(SplitItem *splitInfos_, unsigned splitInfoSize,
	VertItem *vertInfos_, unsigned vertInfoSize,
	PriorityQueues< Window > winQ_, PriorityQueues< PseudoWindow > pseudoSrcQ_,
	Window* storedWindows_, unsigned *keptFaces_)
{
	splitInfos = MyHashTable<SplitInfo>(splitInfoSize, splitInfos_);
	vertInfos = MyHashTable<VertInfo>(vertInfoSize, vertInfos_);
	winQ = winQ_;
	pseudoSrcQ = pseudoSrcQ_;
	storedWindows = storedWindows_;
	keptFaces = keptFaces_;
}

__host__ __device__ void ICH::AddSource(unsigned vertId)
{
	sourceVert = vertId;
}

__host__ __device__ void ICH::AddSource(unsigned faceId, Vector3D pos)
{
	sourcePointFace = faceId;
	sourcePointPos = pos;
}

__host__ __device__ void ICH::AddFacesKeptWindow(unsigned faceId)
{
	if (keptFacesIdx == keptFacesSize) return;
	keptFaces[keptFacesIdx++] = faceId;
}

__host__ __device__ void ICH::Execute(int totalCalcVertNum_)
{
	// Initialize
	Initialize();

	while (!winQ.empty() || !pseudoSrcQ.empty())
	{
		// Get valid window (for window whose pseudoSrcBirthTime is not equal (which means smaller/older) than
		// the current one, it must be an old window, which can be safely skipped)
		/*cout << "\r" << winQ.size() << " " << pseudoSrcQ.size();*/
		maxWinQSize = max(maxWinQSize, winQ.size());
		maxPseudoQSize = max(maxPseudoQSize, pseudoSrcQ.size());

		while (!winQ.empty() && winQ.top().pseudoSrcId < mesh->vertNum &&
			winQ.top().pseudoSrcBirthTime != vertInfos.get(winQ.top().pseudoSrcId).birthTime)
			winQ.pop();

		while (!pseudoSrcQ.empty() && winQ.top().pseudoSrcId < mesh->vertNum &&
			pseudoSrcQ.top().pseudoBirthTime != vertInfos.get(pseudoSrcQ.top().vertID).birthTime)
			pseudoSrcQ.pop();

		if (!winQ.empty() && (pseudoSrcQ.empty() || winQ.top().minDist < pseudoSrcQ.top().dist))
		{
			Window win = winQ.top(); winQ.pop();
			if (win.level > mesh->faceNum) continue;
			// save windows for arbitrary dst geodesic construction
			unsigned twinEdge = mesh->edges[win.edgeID].twinEdge;
			if (twinEdge != -1)
			{
				bool found = false;
				for (int i = 0; i < keptFacesIdx; ++i)
				{
					if (keptFaces[i] == mesh->edges[twinEdge].faceId)
					{
						found = true;
						break;
					}
				}
				if (found && storedWindowsIdx < storedWindowsSize)
					storedWindows[storedWindowsIdx++] = win;
			}
			PropagateWindow(win);
		}
		else if (!pseudoSrcQ.empty() && (winQ.empty() || winQ.top().minDist >= pseudoSrcQ.top().dist))
		{
			PseudoWindow pseudoWin = pseudoSrcQ.top(); pseudoSrcQ.pop();
			if (pseudoWin.level >= mesh->faceNum) continue;
			GenSubWinsForPseudoSrc(pseudoWin);
		}

		if (totalCalcVertNum_ != -1 && totalCalcVertNum >= totalCalcVertNum_)
			break;
	}

	auto vertData = vertInfos.Data();
	unsigned cnt = 0;
	for (int i = 0; i < vertInfos.Size(); ++i)
		if (vertData[i].index != -1) ++cnt;
	cout << "Vertex info hash-table size: " << cnt << "/" << vertInfos.Size() << endl;

	auto splitData = splitInfos.Data();
	cnt = 0;
	for (int i = 0; i < splitInfos.Size(); ++i)
		if (splitData[i].index != -1) ++cnt;
	cout << "Split info hash-table: " << cnt << "/" << splitInfos.Size() << endl;
}

__host__ __device__ void ICH::OutputStatisticInfo()
{
	/*
	cout << "Total generated window number: " << numOfWinGen << endl;
	cout << "Max windows queue size: " << maxWinQSize << endl;
	cout << "Max pseudo-source queue size: " << maxPseudoQSize << endl;
	*/
}

__host__ __device__ double ICH::BuildGeodesicPathTo(unsigned faceId, Vector3D pos, unsigned &srcId,
	unsigned &nextToSrcEdge, double &nextToSrcX, unsigned &nextToDstEdge, double &nextToDstX)
{
	pathPassVert = false;
	// find the window provide the nearest distance
	nextToSrcEdge = -1; nextToDstEdge = -1;
	double minDist = DBL_MAX;
	Window minWin; double xInter; Vector2D pos2D;
	unsigned dstVert = -1;
	bool throughAWindow = true;

	// traverse the surrounded windows
	for (int i = 0; i < storedWindowsIdx; ++i)
	{
		unsigned twinEdge = mesh->edges[storedWindows[i].edgeID].twinEdge;
		if (twinEdge == -1) continue;
		if (mesh->edges[twinEdge].faceId != faceId) continue;

		unsigned e0 = twinEdge;
		unsigned e1 = mesh->edges[e0].nextEdge;
		unsigned e2 = mesh->edges[e1].nextEdge;

		double l0 = mesh->edges[e0].edgeLen;
		double l1 = mesh->edges[e1].edgeLen;
		double l2 = mesh->edges[e2].edgeLen;

		unsigned v0 = mesh->edges[e0].verts[1];
		unsigned v1 = mesh->edges[e0].verts[0];
		unsigned v2 = mesh->edges[e1].verts[1];

		Vector2D p0(0.0, 0.0), p1(l0, 0.0), p2;
		p2.x = (l1*l1 + l0*l0 - l2*l2) / (2.0*l0);
		p2.y = -sqrt(fabs(l1*l1 - p2.x*p2.x));

		// window's pseudo source's 2D planar coordinate
		Vector2D src2D = storedWindows[i].FlatenedSrc();

		// dst point's centroid coordinates
		double a = (pos - mesh->verts[v0].pos).length();
		double b = (pos - mesh->verts[v1].pos).length();
		double c = (pos - mesh->verts[v2].pos).length();

		double s0 = (b + c + l2) / 2.0;
		double s1 = (a + c + l1) / 2.0;
		double s2 = (a + b + l0) / 2.0;

		s0 = sqrt(fabs(s0 * (s0 - b) * (s0 - c) * (s0 - l2)));
		s1 = sqrt(fabs(s1 * (s1 - a) * (s1 - c) * (s1 - l1)));
		s2 = sqrt(fabs(s2 * (s2 - a) * (s2 - b) * (s2 - l0)));

		double w0 = s0 / (s0 + s1 + s2);
		double w1 = s1 / (s0 + s1 + s2);
		double w2 = s2 / (s0 + s1 + s2);

		Vector2D curPos2D = w0 * p0 + w1 * p1 + w2 * p2;

		// calculate the shortest distance
		double curXInter = src2D.x - (curPos2D.x - src2D.x) / (curPos2D.y - src2D.y) * src2D.y;
		double curMinDist = DBL_MAX;
		if (curXInter > storedWindows[i].b0 && curXInter < storedWindows[i].b1)
			curMinDist = (curPos2D - src2D).length() + storedWindows[i].pseudoSrcDist;
		else if (curXInter <= storedWindows[i].b0)
			curMinDist = (curPos2D - Vector2D(storedWindows[i].b0, 0.0)).length() + storedWindows[i].d0 + storedWindows[i].pseudoSrcDist;
		else
			curMinDist = (curPos2D - Vector2D(storedWindows[i].b1, 0.0)).length() + storedWindows[i].d1 + storedWindows[i].pseudoSrcDist;

		if (curMinDist < minDist)
		{
			minDist = curMinDist;
			minWin = storedWindows[i];
			xInter = curXInter;
			pos2D = curPos2D;
		}
	}

	// traverse the surrounded vertices
	for (int i = 0; i < 3; ++i)
	{
		unsigned opVert = mesh->edges[mesh->faces[faceId].edges[i]].verts[0];
		if (mesh->verts[opVert].angle < 2.0 * PI) continue;

		double curDist = (pos - mesh->verts[opVert].pos).length() + vertInfos.get(opVert).dist;
		if (curDist < minDist)
		{
			throughAWindow = false;
			dstVert = opVert;
			minDist = curDist;
		}
	}

	if (minDist == DBL_MAX) {
		return minDist;
	}

	if (!throughAWindow)
	{
		BuildGeodesicPathTo(dstVert, srcId, nextToSrcEdge, nextToSrcX, nextToDstEdge, nextToDstX);
		GeodesicKeyPoint gkp;
		gkp.isVertex = true; gkp.id = dstVert;
		pathPassVert = true;
		nextToDstEdge = mesh->verts[gkp.id].firstEdge; nextToDstX = 0.0;
	}
	else
	{
		// next key point is on an edge
		GeodesicKeyPoint gkp;
		gkp.isVertex = false;
		gkp.id = mesh->edges[minWin.edgeID].twinEdge;
		gkp.pos = mesh->edges[gkp.id].edgeLen - xInter;

		nextToDstEdge = gkp.id; nextToDstX = gkp.pos;

		unsigned enterEdge = gkp.id;
		unsigned opVert = mesh->edges[gkp.id].twinEdge;
		opVert = mesh->edges[mesh->edges[opVert].nextEdge].verts[1];
		double l0 = mesh->edges[gkp.id].edgeLen;
		double l1 = mesh->edges[mesh->edges[gkp.id].nextEdge].edgeLen;
		double l2 = mesh->edges[mesh->edges[mesh->edges[gkp.id].nextEdge].nextEdge].edgeLen;

		Vector2D lastPoint = pos2D, curPoint;
		curPoint.x = l0 - gkp.pos; curPoint.y = 0.0;

		while (minWin.pseudoSrcId < mesh->vertNum && opVert != minWin.pseudoSrcId ||
			minWin.pseudoSrcId >= mesh->vertNum &&
			mesh->edges[mesh->edges[gkp.id].twinEdge].faceId != sourcePointFace)
		{
			// trace back
			unsigned e0 = mesh->edges[gkp.id].twinEdge;
			unsigned e1 = mesh->edges[e0].nextEdge;
			unsigned e2 = mesh->edges[e1].nextEdge;
			double l0 = mesh->edges[e0].edgeLen;
			double l1 = mesh->edges[e1].edgeLen;
			double l2 = mesh->edges[e2].edgeLen;

			Vector2D opVert2D;
			opVert2D.x = (l0*l0 + l2*l2 - l1*l1) / (2.0*l0);
			opVert2D.y = sqrt(fabs(l2*l2 - opVert2D.x*opVert2D.x));

			if (toLeft(opVert2D, lastPoint, curPoint))
			{
				Vector2D p0, p1;
				p0.x = (l2*l2 + l1*l1 - l0*l0) / (2.0*l1);
				p0.y = -sqrt(fabs(l2*l2 - p0.x*p0.x));
				p1.x = l1; p1.y = 0.0;
				Vector2D newlastPoint = gkp.pos / l0 * p0 + (1.0 - gkp.pos / l0) * p1;

				gkp.pos = Intersect(lastPoint, curPoint, Vector2D(l0, 0.0), opVert2D);
				gkp.pos = (1.0 - gkp.pos) * l1;
				gkp.id = e1;
				curPoint.x = l1 - gkp.pos; curPoint.y = 0.0;
				lastPoint = newlastPoint;
			}
			else
			{
				Vector2D p0, p1;
				p0.x = 0.0; p0.y = 0.0;
				p1.x = (l2*l2 + l0*l0 - l1*l1) / (2.0*l2);
				p1.y = -sqrt(fabs(l0*l0 - p1.x*p1.x));
				Vector2D newlastPoint = gkp.pos / l0 * p0 + (1.0 - gkp.pos / l0) * p1;

				gkp.pos = Intersect(lastPoint, curPoint, opVert2D, Vector2D(0.0, 0.0));
				gkp.pos = (1.0 - gkp.pos) * l2;
				gkp.id = e2;
				curPoint.x = l2 - gkp.pos; curPoint.y = 0.0;
				lastPoint = newlastPoint;
			}
			nextToSrcEdge = gkp.id; nextToSrcX = gkp.pos;

			opVert = mesh->edges[gkp.id].twinEdge;
			opVert = mesh->edges[mesh->edges[opVert].nextEdge].verts[1];
		}

		if (minWin.pseudoSrcId >= mesh->vertNum) {
			dstVert = minWin.pseudoSrcId;
			srcId = dstVert;
		}
		else if (vertInfos.get(opVert).dist != 0.0)
		{
			gkp.isVertex = true;
			gkp.id = opVert;
			nextToSrcEdge = mesh->verts[gkp.id].firstEdge; nextToSrcX = 0.0;
			dstVert = opVert;
			pathPassVert = true;
			unsigned nextToDstEdge_; double nextToDstX_;
			BuildGeodesicPathTo(opVert, srcId, nextToSrcEdge, nextToSrcX, nextToDstEdge_, nextToDstX_);
		}
		srcId = minWin.srcID;
	}

	return minDist;
}

__host__ __device__ void ICH::BuildGeodesicPathTo(unsigned vertId, unsigned &srcId,
	unsigned &nextToSrcEdge, double &nextToSrcX, unsigned &nextToDstEdge, double &nextToDstX)
{
	// TODO: build geodesic path from vertex vertId to source
	pathPassVert = false;
	nextToSrcEdge = -1; nextToDstEdge = -1;
	unsigned curVert = vertId;
	GeodesicKeyPoint gkp;
	while (vertInfos.get(curVert).dist != 0.0)
	{
		unsigned enterEdge = vertInfos.get(curVert).enterEdge;
		if (enterEdge == -1)
		{
			// trace back to an arbitrary point
			double curPlanarDist = (sourcePointPos - mesh->verts[curVert].pos).length();
			srcId = mesh->vertNum;
			if (nextToDstEdge == -1)
			{
				nextToDstEdge = mesh->verts[curVert].firstEdge;
				nextToDstX = 0.0;
			}
			else
			{
				nextToSrcEdge = mesh->verts[curVert].firstEdge;
				nextToSrcX = 0.0;
			}
			pathPassVert = true;
			return;
		}
		else if (mesh->edges[enterEdge].verts[0] == curVert)
		{
			// next key point is still a vertex
			unsigned nextVert = mesh->edges[enterEdge].verts[1];
			if (!vertInfos.get(nextVert).isSource)
			{
				gkp.isVertex = true;
				gkp.id = nextVert;

				if (nextToDstEdge == -1)
				{
					nextToDstEdge = mesh->edges[enterEdge].nextEdge;
					nextToDstX = 0.0;
				}
				else
				{
					nextToSrcEdge = mesh->edges[enterEdge].nextEdge;
					nextToSrcX = 0.0;
				}
				pathPassVert = true;
			}
			curVert = nextVert;
		}
		else
		{
			// next key point is on an edge
			gkp.isVertex = false;
			gkp.id = enterEdge; gkp.pos = splitInfos.get(enterEdge).x;

			if (nextToDstEdge == -1)
			{
				nextToDstEdge = gkp.id; nextToDstX = gkp.pos;
			}
			else
			{
				nextToSrcEdge = gkp.id; nextToSrcX = gkp.pos;
			}

			unsigned opVert = mesh->edges[gkp.id].twinEdge;
			opVert = mesh->edges[mesh->edges[opVert].nextEdge].verts[1];
			double l0 = mesh->edges[gkp.id].edgeLen;
			double l1 = mesh->edges[mesh->edges[gkp.id].nextEdge].edgeLen;
			double l2 = mesh->edges[mesh->edges[mesh->edges[gkp.id].nextEdge].nextEdge].edgeLen;

			Vector2D lastPoint, curPoint;
			lastPoint.x = (l1*l1 + l0*l0 - l2*l2) / (2.0*l0);
			lastPoint.y = -sqrt(fabs(l1*l1 - lastPoint.x*lastPoint.x));
			curPoint.x = l0 - gkp.pos; curPoint.y = 0.0;

			while (splitInfos.get(enterEdge).pseudoSrcId < mesh->vertNum &&
				opVert != splitInfos.get(enterEdge).pseudoSrcId ||
				splitInfos.get(enterEdge).pseudoSrcId >= mesh->vertNum &&
				mesh->edges[mesh->edges[gkp.id].twinEdge].faceId != sourcePointFace)
			{
				// trace back
				unsigned e0 = mesh->edges[gkp.id].twinEdge;
				unsigned e1 = mesh->edges[e0].nextEdge;
				unsigned e2 = mesh->edges[e1].nextEdge;
				double l0 = mesh->edges[e0].edgeLen;
				double l1 = mesh->edges[e1].edgeLen;
				double l2 = mesh->edges[e2].edgeLen;

				Vector2D opVert2D;
				opVert2D.x = (l0*l0 + l2*l2 - l1*l1) / (2.0*l0);
				opVert2D.y = sqrt(fabs(l2*l2 - opVert2D.x*opVert2D.x));

				if (toLeft(opVert2D, lastPoint, curPoint))
				{
					Vector2D p0, p1;
					p0.x = (l2*l2 + l1*l1 - l0*l0) / (2.0*l1);
					p0.y = -sqrt(fabs(l2*l2 - p0.x*p0.x));
					p1.x = l1; p1.y = 0.0;
					Vector2D newlastPoint = gkp.pos / l0 * p0 + (1.0 - gkp.pos / l0) * p1;

					gkp.pos = Intersect(lastPoint, curPoint, Vector2D(l0, 0.0), opVert2D);
					gkp.pos = (1.0 - gkp.pos) * l1;
					gkp.id = e1;
					curPoint.x = l1 - gkp.pos; curPoint.y = 0.0;
					lastPoint = newlastPoint;
				}
				else
				{
					Vector2D p0, p1;
					p0.x = 0.0; p0.y = 0.0;
					p1.x = (l2*l2 + l0*l0 - l1*l1) / (2.0*l2);
					p1.y = -sqrt(fabs(l0*l0 - p1.x*p1.x));
					Vector2D newlastPoint = gkp.pos / l0 * p0 + (1.0 - gkp.pos / l0) * p1;

					gkp.pos = Intersect(lastPoint, curPoint, opVert2D, Vector2D(0.0, 0.0));
					gkp.pos = (1.0 - gkp.pos) * l2;
					gkp.id = e2;
					curPoint.x = l2 - gkp.pos; curPoint.y = 0.0;
					lastPoint = newlastPoint;
				}
				nextToSrcEdge = gkp.id; nextToSrcX = gkp.pos;

				opVert = mesh->edges[gkp.id].twinEdge;
				opVert = mesh->edges[mesh->edges[opVert].nextEdge].verts[1];
			}

			if (splitInfos.get(enterEdge).pseudoSrcId >= mesh->vertNum) {
				curVert = splitInfos.get(enterEdge).pseudoSrcId;
				break;
			}
			if (vertInfos.get(opVert).dist != 0.0)
			{
				gkp.isVertex = true;
				gkp.id = opVert;
				nextToSrcEdge = mesh->verts[gkp.id].firstEdge; nextToSrcX = 0.0;
				pathPassVert = true;
			}
			curVert = opVert;
		}
	}
	srcId = curVert;
}

__host__ __device__ double ICH::GetDistanceTo(unsigned vertId)
{
	return vertInfos.get(vertId).dist;
}

__host__ __device__ void ICH::Clear()
{
	winQ.clear(); pseudoSrcQ.clear();

	unsigned size = splitInfos.Size();
	SplitItem *splitData = splitInfos.Data();
	for (int i = 0; i < size + 1; ++i)
	{
		splitData[i].item.dist = DBL_MAX;
		splitData[i].item.x = DBL_MAX;
		splitData[i].item.pseudoSrcId = -1;
		splitData[i].index = -1;
	}

	size = vertInfos.Size();
	VertItem *vertData = vertInfos.Data();
	for (int i = 0; i < size + 1; ++i)
	{
		vertData[i].item.birthTime = -1;
		vertData[i].item.dist = DBL_MAX;
		vertData[i].item.enterEdge = -1;
		vertData[i].item.isSource = false;
		vertData[i].index = -1;
	}
	sourceVert = -1;
	sourcePointFace = -1;

	storedWindowsIdx = 0;
	keptFacesIdx = 0;

	numOfWinGen = 0;
	maxWinQSize = 0;
	maxPseudoQSize = 0;
	totalCalcVertNum = 0;
}

__host__ __device__ void ICH::Initialize()
{
	if (sourceVert != -1)
	{
		int startEdge = mesh->verts[sourceVert].firstEdge;
		int curEdge = startEdge;
		do
		{
			unsigned opEdge = mesh->edges[curEdge].nextEdge;
			Window win;
			win.edgeID = opEdge;
			win.b0 = 0.0; win.b1 = mesh->edges[opEdge].edgeLen;
			win.d0 = mesh->edges[curEdge].edgeLen;
			win.d1 = mesh->edges[mesh->edges[curEdge].prevEdge].edgeLen;
			win.pseudoSrcDist = 0.0; win.calcMinDist();
			win.srcID = sourceVert; win.pseudoSrcId = sourceVert;
			win.pseudoSrcBirthTime = 0;
			win.level = 0;
			winQ.push(win, win.minDist);
			++numOfWinGen;

			unsigned opVert = mesh->edges[curEdge].verts[1];
			vertInfos[opVert].birthTime = 0;
			vertInfos[opVert].dist = mesh->edges[curEdge].edgeLen;
			vertInfos[opVert].enterEdge = mesh->edges[curEdge].twinEdge;

			if (mesh->verts[opVert].angle < 2.0 * PI)
			{
				curEdge = mesh->edges[curEdge].twinEdge;
				if (curEdge != -1) curEdge = mesh->edges[curEdge].nextEdge;
				continue;
			}

			PseudoWindow pseudoWin;
			pseudoWin.vertID = opVert; pseudoWin.dist = mesh->edges[curEdge].edgeLen;
			pseudoWin.srcId = sourceVert; pseudoWin.pseudoSrcId = sourceVert;
			pseudoWin.pseudoBirthTime = vertInfos.get(opVert).birthTime;
			pseudoWin.level = 0;
			pseudoSrcQ.push(pseudoWin, pseudoWin.dist);

			curEdge = mesh->edges[curEdge].twinEdge;
			if (curEdge != -1) curEdge = mesh->edges[curEdge].nextEdge;

		} while (curEdge != startEdge && curEdge != -1);
		vertInfos[sourceVert].birthTime = 0;
		vertInfos[sourceVert].dist = 0.0;
		vertInfos[sourceVert].enterEdge = -1;
		vertInfos[sourceVert].isSource = true;
	}

	if (sourcePointFace != -1)
	{
		for (int j = 0; j < 3; ++j)
		{
			unsigned opEdge = mesh->faces[sourcePointFace].edges[j];
			Window win;
			win.edgeID = opEdge;
			win.b0 = 0.0; win.b1 = mesh->edges[opEdge].edgeLen;
			win.d0 = (sourcePointPos - mesh->verts[mesh->edges[opEdge].verts[0]].pos).length();
			win.d1 = (sourcePointPos - mesh->verts[mesh->edges[opEdge].verts[1]].pos).length();
			win.pseudoSrcDist = 0.0; win.calcMinDist();
			win.srcID = mesh->vertNum; win.pseudoSrcId = win.srcID;
			win.pseudoSrcBirthTime = 0; win.level = 0;
			winQ.push(win, win.minDist);

			unsigned opVert = mesh->edges[opEdge].verts[0];
			vertInfos[opVert].birthTime = 0;
			vertInfos[opVert].dist = (sourcePointPos - mesh->verts[opVert].pos).length();
			vertInfos[opVert].enterEdge = -1;

			if (mesh->verts[opVert].angle < 2.0 * PI) continue;

			PseudoWindow pseudoWin;
			pseudoWin.vertID = opVert;
			pseudoWin.dist = (mesh->verts[opVert].pos - sourcePointPos).length();
			pseudoWin.srcId = win.srcID; pseudoWin.pseudoSrcId = win.srcID;
			pseudoWin.pseudoBirthTime = vertInfos.get(opVert).birthTime;
			pseudoWin.level = 0;
			pseudoSrcQ.push(pseudoWin, pseudoWin.dist);
		}
	}
}

__host__ __device__ void ICH::PropagateWindow(const Window &win)
{
	unsigned e0 = mesh->edges[win.edgeID].twinEdge;
	if (e0 == -1) return;
	unsigned e1 = mesh->edges[e0].nextEdge;
	unsigned e2 = mesh->edges[e1].nextEdge;
	unsigned opVert = mesh->edges[e1].verts[1];

	Vector2D src2D = win.FlatenedSrc();

	Vector2D left(win.b0, 0.0), right(win.b1, 0.0);
	double l0 = mesh->edges[e0].edgeLen;
	double l1 = mesh->edges[e1].edgeLen;
	double l2 = mesh->edges[e2].edgeLen;
	Vector2D v0(0.0, 0.0), v1(l0, 0.0), v2;
	v2.x = (l1*l1 + l0*l0 - l2*l2) / (2.0 * l0);
	v2.y = -sqrt(fabs(l1*l1 - v2.x*v2.x));

	double interX = v2.x - v2.y * (v2.x - src2D.x) / (v2.y - src2D.y);
	Window leftChildWin, rightChildWin;
	bool hasLeftChild = true, hasRightChild = true;
	// only generate right window
	if (interX <= left.x)
	{
		hasLeftChild = false;
		double t0 = Intersect(src2D, left, v2, v1);
		double t1 = Intersect(src2D, right, v2, v1);
		BuildWindow(win, e2, t0, t1, v2, v1, rightChildWin);
		if (!IsValidWindow(rightChildWin, false)) hasRightChild = false;
	}
	// only generate left window
	else if (interX >= right.x)
	{
		hasRightChild = false;
		double t0 = Intersect(src2D, left, v0, v2);
		double t1 = Intersect(src2D, right, v0, v2);
		BuildWindow(win, e1, t0, t1, v0, v2, leftChildWin);
		if (!IsValidWindow(leftChildWin, true)) hasLeftChild = false;
	}
	// generate both left and right window
	else
	{
		double directDist = (v2 - src2D).length();
		// ONE ANGLE, ONE SPLIT
		if (directDist + win.pseudoSrcDist > splitInfos.get(e0).dist &&
			(directDist + win.pseudoSrcDist) / splitInfos.get(e0).dist - 1.0 > RELATIVE_ERROR)
		{
			hasLeftChild = splitInfos.get(e0).x < interX;
			hasRightChild = !hasLeftChild;
			/*cout << "Filter 1 works..." << endl;*/
		}
		else
		{
			splitInfos[e0].dist = directDist + win.pseudoSrcDist;
			splitInfos[e0].pseudoSrcId = win.pseudoSrcId;
			splitInfos[e0].x = l0 - interX;

			if (directDist + win.pseudoSrcDist < vertInfos.get(opVert).dist)
			{
				if (vertInfos.get(opVert).dist == DBL_MAX
					/*&& (win.pseudoSrcId == sourceVert || win.pseudoSrcId == mesh->vertNum)*/)
					++totalCalcVertNum;

				++vertInfos[opVert].birthTime;
				vertInfos[opVert].dist = directDist + win.pseudoSrcDist;
				vertInfos[opVert].enterEdge = e0;
				if (mesh->verts[opVert].angle > 2.0 * PI)
				{
					PseudoWindow pseudoWin;
					pseudoWin.vertID = opVert; pseudoWin.dist = vertInfos.get(opVert).dist;
					pseudoWin.srcId = win.srcID; pseudoWin.pseudoSrcId = win.pseudoSrcId;
					pseudoWin.pseudoBirthTime = vertInfos.get(opVert).birthTime;
					pseudoWin.level = win.level + 1;
					pseudoSrcQ.push(pseudoWin, pseudoWin.dist);
				}
			}
		}
		if (hasLeftChild)
		{
			// left child window
			double t0 = Intersect(src2D, left, v0, v2);
			BuildWindow(win, e1, t0, 0.0, v0, v2, leftChildWin);
			if (!IsValidWindow(leftChildWin, true)) hasLeftChild = false;
		}
		if (hasRightChild)
		{
			// right child window
			double t1 = Intersect(src2D, right, v2, v1);
			BuildWindow(win, e2, 1.0, t1, v2, v1, rightChildWin);
			if (!IsValidWindow(rightChildWin, false)) hasRightChild = false;
		}
	}

	if (hasLeftChild)
	{
		++numOfWinGen;
		winQ.push(leftChildWin, leftChildWin.minDist);
	}
	if (hasRightChild)
	{
		++numOfWinGen;
		winQ.push(rightChildWin, rightChildWin.minDist);
	}

}

__host__ __device__ void ICH::GenSubWinsForPseudoSrc(const PseudoWindow &pseudoWin)
{
	unsigned startEdge, endEdge;
	if (vertInfos.get(pseudoWin.vertID).enterEdge == -1 && vertInfos.get(pseudoWin.vertID).birthTime != -1)
	{
		startEdge = mesh->verts[pseudoWin.vertID].firstEdge;
		endEdge = startEdge;
	}
	else if (mesh->edges[vertInfos.get(pseudoWin.vertID).enterEdge].verts[0] == pseudoWin.vertID)
		GenSubWinsForPseudoSrcFromPseudoSrc(pseudoWin, startEdge, endEdge);
	else if (mesh->edges[mesh->edges[vertInfos.get(pseudoWin.vertID).enterEdge].nextEdge].verts[1] == pseudoWin.vertID)
		GenSubWinsForPseudoSrcFromWindow(pseudoWin, startEdge, endEdge);
	else assert(false);

	// generate windows
	do
	{
		Window win;
		win.edgeID = mesh->edges[startEdge].nextEdge;
		win.b0 = 0.0; win.b1 = mesh->edges[win.edgeID].edgeLen;
		win.d0 = mesh->edges[startEdge].edgeLen;
		win.d1 = mesh->edges[mesh->edges[win.edgeID].nextEdge].edgeLen;
		win.pseudoSrcDist = pseudoWin.dist; win.calcMinDist();
		win.srcID = pseudoWin.srcId; win.pseudoSrcId = pseudoWin.vertID;
		win.pseudoSrcBirthTime = pseudoWin.pseudoBirthTime;
		win.level = pseudoWin.level + 1;
		winQ.push(win, win.minDist);
		++numOfWinGen;

		startEdge = mesh->edges[mesh->edges[mesh->edges[startEdge].nextEdge].nextEdge].twinEdge;
	} while (startEdge != endEdge);

	// generate adjacent pseudo sources
	startEdge = mesh->verts[pseudoWin.vertID].firstEdge;
	int curEdge = startEdge;
	do
	{
		unsigned opVert = mesh->edges[curEdge].verts[1];
		if (mesh->verts[opVert].angle < 2.0 * PI ||
			vertInfos.get(opVert).dist < pseudoWin.dist + mesh->edges[curEdge].edgeLen)
		{
			curEdge = mesh->edges[curEdge].twinEdge;
			if (curEdge != -1) curEdge = mesh->edges[curEdge].nextEdge;
			continue;
		}

		if (vertInfos.get(opVert).dist == DBL_MAX
			/*&& (pseudoWin.pseudoSrcId == sourceVert || pseudoWin.pseudoSrcId == mesh->vertNum)*/)
			++totalCalcVertNum;

		vertInfos[opVert].dist = pseudoWin.dist + mesh->edges[curEdge].edgeLen;
		++vertInfos[opVert].birthTime;
		vertInfos[opVert].enterEdge = mesh->edges[curEdge].twinEdge;

		PseudoWindow childPseudoWin;
		childPseudoWin.vertID = opVert; childPseudoWin.dist = vertInfos.get(opVert).dist;
		childPseudoWin.srcId = pseudoWin.srcId; childPseudoWin.pseudoSrcId = pseudoWin.vertID;
		childPseudoWin.pseudoBirthTime = vertInfos.get(opVert).birthTime;
		childPseudoWin.level = pseudoWin.level;
		pseudoSrcQ.push(childPseudoWin, childPseudoWin.dist);

		curEdge = mesh->edges[curEdge].twinEdge;
		if (curEdge != -1) curEdge = mesh->edges[curEdge].nextEdge;
	} while (curEdge != startEdge && curEdge != -1);
}

__host__ __device__ void ICH::GenSubWinsForPseudoSrcFromWindow(const PseudoWindow &pseudoWin, unsigned &startEdge, unsigned &endEdge)
{
	unsigned e0 = vertInfos.get(pseudoWin.vertID).enterEdge;
	unsigned e1 = mesh->edges[e0].nextEdge;
	unsigned e2 = mesh->edges[e1].nextEdge;

	double l0 = mesh->edges[e0].edgeLen;
	double l1 = mesh->edges[e1].edgeLen;
	double l2 = mesh->edges[e2].edgeLen;

	unsigned pseudoSrc = pseudoWin.vertID;
	Vector2D enterPoint;
	enterPoint.x = l0 - splitInfos.get(e0).x;
	enterPoint.y = 0.0;

	Vector2D v0(0.0, 0.0), v1(l0, 0.0), v2;
	v2.x = (l1*l1 + l0*l0 - l2*l2) / (2.0*l0);
	v2.y = -sqrt(fabs(l1*l1 - v2.x*v2.x));

	// TODO: generate windows using opVert as pseudo sources
	double angle0 = (enterPoint - v2) * (v0 - v2) / (enterPoint - v2).length() / l1;
	double angle1 = (enterPoint - v2) * (v1 - v2) / (enterPoint - v2).length() / l2;
	if (angle0 > 1.0) angle0 = 1.0; else if (angle0 < -1.0) angle0 = -1.0;
	if (angle1 > 1.0) angle1 = 1.0; else if (angle1 < -1.0) angle1 = -1.0;
	angle0 = acos(angle0); angle1 = acos(angle1);

	startEdge = -1, endEdge = -1;
	// traverse from left
	unsigned curEdge = mesh->edges[e1].twinEdge;
	while (angle0 < PI && curEdge != -1)
	{
		unsigned opEdge = mesh->edges[curEdge].nextEdge;
		unsigned nextEdge = mesh->edges[opEdge].nextEdge;
		double L0 = mesh->edges[curEdge].edgeLen, L1 = mesh->edges[nextEdge].edgeLen;
		double L2 = mesh->edges[opEdge].edgeLen;
		double curAngle = (L0*L0 + L1*L1 - L2*L2) / (2.0 * L0 * L1);
		if (curAngle > 1.0) curAngle = 1.0; else if (curAngle < -1.0) curAngle = -1.0;
		curAngle = acos(curAngle);
		angle0 += curAngle;
		curEdge = mesh->edges[nextEdge].twinEdge;
	}
	if (curEdge != -1)
		startEdge = mesh->edges[mesh->edges[curEdge].twinEdge].nextEdge;

	// traverse from right
	curEdge = mesh->edges[e2].twinEdge;
	while (angle1 < PI && curEdge != -1)
	{
		unsigned nextEdge = mesh->edges[curEdge].nextEdge;
		unsigned opEdge = mesh->edges[nextEdge].nextEdge;
		double L0 = mesh->edges[curEdge].edgeLen, L1 = mesh->edges[nextEdge].edgeLen;
		double L2 = mesh->edges[opEdge].edgeLen;
		double curAngle = (L0*L0 + L1*L1 - L2*L2) / (2.0 * L0 * L1);
		if (curAngle > 1.0) curAngle = 1.0; else if (curAngle < -1.0) curAngle = -1.0;
		curAngle = acos(curAngle);
		angle1 += curAngle;
		curEdge = mesh->edges[nextEdge].twinEdge;
	}
	if (curEdge != -1)
	{
		endEdge = mesh->edges[mesh->edges[mesh->edges[curEdge].twinEdge].nextEdge].nextEdge;
		endEdge = mesh->edges[endEdge].twinEdge;
	}
}

__host__ __device__ void ICH::GenSubWinsForPseudoSrcFromPseudoSrc(const PseudoWindow &pseudoWin, unsigned &startEdge, unsigned &endEdge)
{
	unsigned pseudoSrc = pseudoWin.vertID;

	// TODO: generate windows using opVert as pseudo sources
	double angle0 = 0.0, angle1 = 0.0;

	startEdge = -1, endEdge = -1;
	// traverse from left
	unsigned curEdge = vertInfos.get(pseudoWin.vertID).enterEdge;
	while (angle0 < PI && curEdge != -1)
	{
		unsigned opEdge = mesh->edges[curEdge].nextEdge;
		unsigned nextEdge = mesh->edges[opEdge].nextEdge;
		double L0 = mesh->edges[curEdge].edgeLen, L1 = mesh->edges[nextEdge].edgeLen;
		double L2 = mesh->edges[opEdge].edgeLen;
		double curAngle = (L0*L0 + L1*L1 - L2*L2) / (2.0 * L0 * L1);
		if (curAngle > 1.0) curAngle = 1.0; else if (curAngle < -1.0) curAngle = -1.0;
		curAngle = acos(curAngle);
		angle0 += curAngle;
		curEdge = mesh->edges[nextEdge].twinEdge;
	}
	if (curEdge != -1)
		startEdge = mesh->edges[mesh->edges[curEdge].twinEdge].nextEdge;

	// traverse from right
	curEdge = mesh->edges[vertInfos.get(pseudoWin.vertID).enterEdge].twinEdge;
	while (angle1 < PI && curEdge != -1)
	{
		unsigned nextEdge = mesh->edges[curEdge].nextEdge;
		unsigned opEdge = mesh->edges[nextEdge].nextEdge;
		double L0 = mesh->edges[curEdge].edgeLen, L1 = mesh->edges[nextEdge].edgeLen;
		double L2 = mesh->edges[opEdge].edgeLen;
		double curAngle = (L0*L0 + L1*L1 - L2*L2) / (2.0 * L0 * L1);
		if (curAngle > 1.0) curAngle = 1.0; else if (curAngle < -1.0) curAngle = -1.0;
		curAngle = acos(curAngle);
		angle1 += curAngle;
		curEdge = mesh->edges[nextEdge].twinEdge;
	}
	if (curEdge != -1)
	{
		endEdge = mesh->edges[mesh->edges[mesh->edges[curEdge].twinEdge].nextEdge].nextEdge;
		endEdge = mesh->edges[endEdge].twinEdge;
	}
}

__host__ __device__ bool ICH::IsValidWindow(const Window &win, bool isLeftChild)
{
	if (win.b1 <= win.b0) return false;
	// apply ICH's filter
	unsigned v1 = mesh->edges[win.edgeID].verts[0];
	unsigned v2 = mesh->edges[win.edgeID].verts[1];
	unsigned v3 = mesh->edges[mesh->edges[win.edgeID].nextEdge].verts[1];
	double l0 = mesh->edges[win.edgeID].edgeLen;
	double l1 = mesh->edges[mesh->edges[win.edgeID].nextEdge].edgeLen;
	double l2 = mesh->edges[mesh->edges[mesh->edges[win.edgeID].nextEdge].nextEdge].edgeLen;
	Vector2D p1(0.0, 0.0), p2(l0, 0.0), p3;
	p3.x = (l2*l2 + l0*l0 - l1*l1) / (2.0 * l0);
	p3.y = sqrt(fabs(l2*l2 - p3.x*p3.x));

	Vector2D A(win.b0, 0.0), B(win.b1, 0.0);
	Vector2D src2D = win.FlatenedSrc();


	if (win.pseudoSrcDist + (src2D - B).length() > vertInfos.get(v1).dist + win.b1 &&
		(win.pseudoSrcDist + (src2D - B).length()) / (vertInfos.get(v1).dist + win.b1) - 1.0 > 0.0)
	{
		/*cout << "Filter 2 works..." << endl;*/
		return false;
	}
	if (win.pseudoSrcDist + (src2D - A).length() > vertInfos.get(v2).dist + l0 - win.b0 &&
		(win.pseudoSrcDist + (src2D - A).length()) / (vertInfos.get(v2).dist + l0 - win.b0) - 1.0 > 0.0)
	{
		/*cout << "Filter 2 works..." << endl;*/
		return false;
	}
	if (isLeftChild)
	{
		if (win.pseudoSrcDist + (src2D - A).length() > vertInfos.get(v3).dist + (p3 - A).length() &&
			(win.pseudoSrcDist + (src2D - A).length()) / (vertInfos.get(v3).dist + (p3 - A).length()) - 1.0 > 0.0)
		{
			/*cout << "Filter 2 works..." << endl;*/
			return false;
		}
	}
	else
	{
		if (win.pseudoSrcDist + (src2D - B).length() > vertInfos.get(v3).dist + (p3 - B).length() &&
			(win.pseudoSrcDist + (src2D - B).length()) / (vertInfos.get(v3).dist + (p3 - B).length()) - 1.0 > RELATIVE_ERROR)
		{
			/*cout << "Filter 2 works..." << endl;*/
			return false;
		}
	}
	return true;
}

__host__ __device__ void ICH::BuildWindow(const Window &fatherWin,
	unsigned edge,
	double t0, double t1,
	const Vector2D &v0, const Vector2D &v1,
	Window &win)
{
	Vector2D src2D = fatherWin.FlatenedSrc();
	win.edgeID = edge;
	win.b0 = (1 - t0) * mesh->edges[edge].edgeLen; win.b1 = (1 - t1) * mesh->edges[edge].edgeLen;
	win.d0 = (src2D - (t0 * v0 + (1 - t0)*v1)).length();
	win.d1 = (src2D - (t1 * v0 + (1 - t1)*v1)).length();
	win.pseudoSrcDist = fatherWin.pseudoSrcDist;
	win.calcMinDist();
	win.srcID = fatherWin.srcID; win.pseudoSrcId = fatherWin.pseudoSrcId;
	win.pseudoSrcBirthTime = fatherWin.pseudoSrcBirthTime;
	win.level = fatherWin.level + 1;
}

__host__ __device__ double ICH::Intersect(const Vector2D &v0, const Vector2D &v1, const Vector2D &p0, const Vector2D &p1)
{
	double a00 = p0.x - p1.x, a01 = v1.x - v0.x, b0 = v1.x - p1.x;
	double a10 = p0.y - p1.y, a11 = v1.y - v0.y, b1 = v1.y - p1.y;
	return (b0*a11 - b1*a01) / (a00*a11 - a10*a01);
}