// ICH.cpp : Defines the exported functions for the DLL application.
//

#include "ICH.cuh"
#include <iostream>
#include <list>

using namespace std;

// This is the constructor of a class that has been exported.
// see ICH.h for the class definition
__device__ ICH::ICH()
{
	numOfWinGen = 0;
	maxWinQSize = 0;
	maxPseudoQSize = 0;
	return;
}

__device__ ICH::~ICH()
{
	return;
}

__device__ void ICH::AssignMesh(Mesh *mesh_)
{
	mesh = mesh_;
}

__device__ void ICH::AssignBuffers(SplitInfo *splitInfos_, VertInfo *vertInfos_)
{
	splitInfos = splitInfos_;
	vertInfos = vertInfos_;
}

__device__ void ICH::AddSource(unsigned vertId)
{
	sourceVert = vertId;
}

__device__ void ICH::Execute()
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

		while (!winQ.empty() &&
			winQ.top().pseudoSrcBirthTime != vertInfos[winQ.top().pseudoSrcId].birthTime)
			winQ.pop();

		while (!pseudoSrcQ.empty() &&
			pseudoSrcQ.top().pseudoBirthTime != vertInfos[pseudoSrcQ.top().vertID].birthTime)
			pseudoSrcQ.pop();

		if (!winQ.empty() && (pseudoSrcQ.empty() || winQ.top().minDist < pseudoSrcQ.top().dist))
		{
			Window win = winQ.top(); winQ.pop();
			if (win.level > mesh->faceNum) continue;
			PropagateWindow(win);
		}
		else if (!pseudoSrcQ.empty() && (winQ.empty() || winQ.top().minDist >= pseudoSrcQ.top().dist))
		{
			PseudoWindow pseudoWin = pseudoSrcQ.top(); pseudoSrcQ.pop();
			if (pseudoWin.level >= mesh->faceNum) continue;
			GenSubWinsForPseudoSrc(pseudoWin);
		}
	}
}

__device__ void ICH::OutputStatisticInfo()
{
	/*
	cout << "Total generated window number: " << numOfWinGen << endl;
	cout << "Max windows queue size: " << maxWinQSize << endl;
	cout << "Max pseudo-source queue size: " << maxPseudoQSize << endl;
	*/
}

/*
list<ICH::GeodesicKeyPoint> ICH::BuildGeodesicPathTo(unsigned vertId, unsigned &srcId)
{
	// TODO: build geodesic path from vertex vertId to source
	list < GeodesicKeyPoint > path;
	unsigned curVert = vertId;
	GeodesicKeyPoint gkp;
	while (vertInfos[curVert].dist != 0.0)
	{
		unsigned enterEdge = vertInfos[curVert].enterEdge;
		if (mesh->edges[enterEdge].m_iVertex[0] == curVert)
		{
			// next key point is still a vertex
			gkp.isVertex = true;
			gkp.id = mesh->edges[enterEdge].m_iVertex[1];
			path.push_back(gkp);
			curVert = gkp.id;
		}
		else
		{
			// next key point is on an edge
			gkp.isVertex = false;
			gkp.id = enterEdge; gkp.pos = splitInfos[enterEdge].x;
			path.push_back(gkp);

			unsigned opVert = mesh->edges[gkp.id].twinEdge;
			opVert = mesh->edges[mesh->edges[opVert].nextEdge].m_iVertex[1];
			double l0 = mesh->edges[gkp.id].edgeLen;
			double l1 = mesh->edges[mesh->edges[gkp.id].nextEdge].edgeLen;
			double l2 = mesh->edges[mesh->edges[mesh->edges[gkp.id].nextEdge].nextEdge].edgeLen;
			
			Vector2D lastPoint, curPoint;
			lastPoint.x = (l1*l1 + l0*l0 - l2*l2) / (2.0*l0);
			lastPoint.y = -sqrt(fabs(l1*l1 - lastPoint.x*lastPoint.x));
			curPoint.x = l0 - gkp.pos; curPoint.y = 0.0;

			while (opVert != splitInfos[enterEdge].pseudoSrcId)
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
				path.push_back(gkp);

				opVert = mesh->edges[gkp.id].twinEdge;
				opVert = mesh->edges[mesh->edges[opVert].nextEdge].m_iVertex[1];
			}

			if (vertInfos[opVert].dist != 0.0)
			{
				gkp.isVertex = true;
				gkp.id = opVert;
				path.push_back(gkp);
			}
			curVert = opVert;
		}
	}
	srcId = curVert;
	return path;
}
*/

__device__ double ICH::GetDistanceTo(unsigned vertId)
{
	return vertInfos[vertId].dist;
}

__device__ void ICH::Initialize()
{
	int startEdge = mesh->edgeAdjToVert[sourceVert];
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

		if (mesh->angles[opVert] < 2.0 * PI) continue;

		PseudoWindow pseudoWin;
		pseudoWin.vertID = opVert; pseudoWin.dist = mesh->edges[curEdge].edgeLen;
		pseudoWin.srcId = sourceVert; pseudoWin.pseudoSrcId = sourceVert;
		pseudoWin.pseudoBirthTime = vertInfos[opVert].birthTime;
		pseudoWin.level = 0;
		pseudoSrcQ.push(pseudoWin, pseudoWin.dist);
	} while (curEdge != startEdge);
	vertInfos[sourceVert].birthTime = 0;
	vertInfos[sourceVert].dist = 0.0;
	vertInfos[sourceVert].enterEdge = -1;
}

__device__ void ICH::PropagateWindow(const Window &win)
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
		if (directDist + win.pseudoSrcDist > splitInfos[e0].dist && 
			(directDist + win.pseudoSrcDist) / splitInfos[e0].dist - 1.0 > RELATIVE_ERROR)
		{
			hasLeftChild = splitInfos[e0].x < interX;
			hasRightChild = !hasLeftChild;
			/*cout << "Filter 1 works..." << endl;*/
		}
		else
		{
			splitInfos[e0].dist = directDist + win.pseudoSrcDist;
			splitInfos[e0].x = l0 - interX;

			if (directDist + win.pseudoSrcDist < vertInfos[opVert].dist)
			{
				++vertInfos[opVert].birthTime;
				vertInfos[opVert].dist = directDist + win.pseudoSrcDist;
				vertInfos[opVert].enterEdge = e0;
				if (mesh->angles[opVert] > 2.0 * PI)
				{
					PseudoWindow pseudoWin;
					pseudoWin.vertID = opVert; pseudoWin.dist = vertInfos[opVert].dist;
					pseudoWin.srcId = win.srcID; pseudoWin.pseudoSrcId = win.pseudoSrcId;
					pseudoWin.pseudoBirthTime = vertInfos[opVert].birthTime;
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

__device__ void ICH::GenSubWinsForPseudoSrc(const PseudoWindow &pseudoWin)
{
	unsigned startEdge, endEdge;
	if (mesh->edges[vertInfos[pseudoWin.vertID].enterEdge].verts[0] == pseudoWin.vertID)
		GenSubWinsForPseudoSrcFromPseudoSrc(pseudoWin, startEdge, endEdge);
	else if (mesh->edges[mesh->edges[vertInfos[pseudoWin.vertID].enterEdge].nextEdge].verts[1] == pseudoWin.vertID)
		GenSubWinsForPseudoSrcFromWindow(pseudoWin, startEdge, endEdge);
	else assert(false);

	// generate windows
	while (startEdge != endEdge)
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
	}

	// generate adjacent pseudo sources
	startEdge = mesh->edgeAdjToVert[pseudoWin.vertID];
	int curEdge = startEdge;
	do
	{
		unsigned opVert = mesh->edges[curEdge].verts[1];
		if (mesh->angles[opVert] < 2.0 * PI) continue;
		if (vertInfos[opVert].dist < pseudoWin.dist + mesh->edges[curEdge].edgeLen) continue;

		vertInfos[opVert].dist = pseudoWin.dist + mesh->edges[curEdge].edgeLen;
		++vertInfos[opVert].birthTime;
		vertInfos[opVert].enterEdge = mesh->edges[curEdge].twinEdge;

		PseudoWindow childPseudoWin;
		childPseudoWin.vertID = opVert; childPseudoWin.dist = vertInfos[opVert].dist;
		childPseudoWin.srcId = pseudoWin.srcId; childPseudoWin.pseudoSrcId = pseudoWin.vertID;
		childPseudoWin.pseudoBirthTime = vertInfos[opVert].birthTime;
		childPseudoWin.level = pseudoWin.level;
		pseudoSrcQ.push(childPseudoWin, childPseudoWin.dist);
	} while (curEdge != startEdge);
}

__device__ void ICH::GenSubWinsForPseudoSrcFromWindow(const PseudoWindow &pseudoWin, unsigned &startEdge, unsigned &endEdge)
{
	unsigned e0 = vertInfos[pseudoWin.vertID].enterEdge;
	unsigned e1 = mesh->edges[e0].nextEdge;
	unsigned e2 = mesh->edges[e1].nextEdge;

	double l0 = mesh->edges[e0].edgeLen;
	double l1 = mesh->edges[e1].edgeLen;
	double l2 = mesh->edges[e2].edgeLen;

	unsigned pseudoSrc = pseudoWin.vertID;
	Vector2D enterPoint;
	enterPoint.x = l0 - splitInfos[e0].x;
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
	while (angle0 < PI || curEdge == -1)
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
	while (angle1 < PI || curEdge == -1)
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

__device__ void ICH::GenSubWinsForPseudoSrcFromPseudoSrc(const PseudoWindow &pseudoWin, unsigned &startEdge, unsigned &endEdge)
{
	unsigned pseudoSrc = pseudoWin.vertID;

	// TODO: generate windows using opVert as pseudo sources
	double angle0 = 0.0, angle1 = 0.0;

	startEdge = -1, endEdge = -1;
	// traverse from left
	unsigned curEdge = vertInfos[pseudoWin.vertID].enterEdge;
	while (angle0 < PI || curEdge == -1)
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
	curEdge = mesh->edges[vertInfos[pseudoWin.vertID].enterEdge].twinEdge;
	while (angle1 < PI || curEdge == -1)
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

__device__ bool ICH::IsValidWindow(const Window &win, bool isLeftChild)
{
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
	

	if (win.pseudoSrcDist + (src2D - B).length() > vertInfos[v1].dist + win.b1 && 
		(win.pseudoSrcDist + (src2D - B).length()) / (vertInfos[v1].dist + win.b1) - 1.0 > 0.0)
	{
		/*cout << "Filter 2 works..." << endl;*/
		return false;
	}
	if (win.pseudoSrcDist + (src2D - A).length() > vertInfos[v2].dist + l0 - win.b0 && 
		(win.pseudoSrcDist + (src2D - A).length()) / (vertInfos[v2].dist + l0 - win.b0) - 1.0 > 0.0)
	{
		/*cout << "Filter 2 works..." << endl;*/
		return false;
	}
	if (isLeftChild)
	{
		if (win.pseudoSrcDist + (src2D - A).length() > vertInfos[v3].dist + (p3 - A).length() && 
			(win.pseudoSrcDist + (src2D - A).length()) / (vertInfos[v3].dist + (p3 - A).length()) - 1.0 > 0.0)
		{
			/*cout << "Filter 2 works..." << endl;*/
			return false;
		}
	}
	else
	{
		if (win.pseudoSrcDist + (src2D - B).length() > vertInfos[v3].dist + (p3 - B).length() && 
			(win.pseudoSrcDist + (src2D - B).length()) / (vertInfos[v3].dist + (p3 - B).length()) - 1.0 > RELATIVE_ERROR)
		{
			/*cout << "Filter 2 works..." << endl;*/
			return false;
		}
	}
	return true;
}

__device__ void ICH::BuildWindow(const Window &fatherWin,
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

__device__ double ICH::Intersect(const Vector2D &v0, const Vector2D &v1, const Vector2D &p0, const Vector2D &p1)
{
	double a00 = p0.x - p1.x, a01 = v1.x - v0.x, b0 = v1.x - p1.x;
	double a10 = p0.y - p1.y, a11 = v1.y - v0.y, b1 = v1.y - p1.y;
	return (b0*a11 - b1*a01) / (a00*a11 - a10*a01);
}