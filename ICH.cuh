#ifndef ICH_CUH
#define ICH_CUH

#include "Mesh.cuh"
#include <vector>
#include "PriorityQueue.cuh"
#include <math.h>

#define RELATIVE_ERROR 1e-8
//#define L_RELATIVE_ERROR 1e-3

// This class is exported from the ICH.dll
class ICH {
public:
	struct Window	// 68 bytes
	{
		unsigned edgeID;
		double b0, b1, d0, d1;
		double pseudoSrcDist, minDist;
		unsigned srcID, pseudoSrcId;
		int pseudoSrcBirthTime;
		int level;

		bool operator< (const Window& right) const
		{
			return minDist > right.minDist;
		}

		__host__ __device__ void calcMinDist()
		{
			double wLen = b1 - b0;
			double xProj = (d0*d0 + wLen*wLen - d1*d1) / (2 * wLen);
			if (xProj < 0.0) minDist = d0 + pseudoSrcDist;
			else if (xProj > wLen) minDist = d1 + pseudoSrcDist;
			else {
				minDist = sqrt(fabs(d0*d0 - xProj*xProj)) + pseudoSrcDist;
			}
		}

		__host__ __device__ Vector2D FlatenedSrc() const
		{
			Vector2D src2D;
			double wLen = b1 - b0;
			src2D.x = (d0*d0 + wLen*wLen - d1*d1) / (2.0 * wLen);
			src2D.y = sqrt(fabs(d0*d0 - src2D.x*src2D.x));
			src2D.x += b0;
			return src2D;
		}
	};

	struct PseudoWindow	// 28 bytes
	{
		unsigned vertID;
		double dist;
		unsigned srcId, pseudoSrcId;
		unsigned pseudoBirthTime;
		unsigned level;

		__host__ __device__ bool operator< (const PseudoWindow &right) const {
			return dist > right.dist;
		};
	};

	struct SplitInfo // 20 bytes
	{
		double dist;
		double x;
		unsigned pseudoSrcId;

		__host__ __device__ SplitInfo() { dist = DBL_MAX; x = DBL_MAX; pseudoSrcId = -1; }
	};

	struct VertInfo	// 13 bytes
	{
		char birthTime;
		double dist;
		int enterEdge;

		__host__ __device__ VertInfo() { birthTime = -1; dist = DBL_MAX; enterEdge = -1; }
	};

	struct GeodesicKeyPoint
	{
		bool isVertex;
		unsigned id;
		double pos;
	};

public:
	__host__ __device__ ICH();
	__host__ __device__ ~ICH();

	__host__ __device__ void AssignMesh(Mesh *mesh_);
	__host__ __device__ void AssignBuffers(SplitInfo *splitInfos_, VertInfo *vertInfos_,
		PriorityQueues< Window > winQ_, PriorityQueues< PseudoWindow > pseudoSrcQ_,
		Window* storedWindows_, unsigned *keptFaces_);
	__host__ __device__ void AddSource(unsigned vertId);
	__host__ __device__ void AddSource(unsigned faceId, Vector3D pos);
	__host__ __device__ void AddFacesKeptWindow(unsigned faceId);
	__host__ __device__ void Execute(int totalCalcVertNum_ = -1);
	__host__ __device__ void OutputStatisticInfo();
	__host__ __device__ void BuildGeodesicPathTo(unsigned vertId, unsigned &srcId,
		unsigned &nextToSrcEdge, double &nextToSrcX, unsigned &nextToDstEdge, double &nextToDstX);
	__host__ __device__ void BuildGeodesicPathTo(unsigned faceId, Vector3D pos, unsigned &srcId,
		unsigned &nextToSrcEdge, double &nextToSrcX, unsigned &nextToDstEdge, double &nextToDstX);
	__host__ __device__ double GetDistanceTo(unsigned vertId);
	__host__ __device__ void Clear();

private:
	__host__ __device__ void Initialize();
	__host__ __device__ void PropagateWindow(const Window &win);

	__host__ __device__ void GenSubWinsForPseudoSrc(const PseudoWindow &pseudoWin);
	__host__ __device__ void GenSubWinsForPseudoSrcFromWindow(const PseudoWindow &pseudoWin, unsigned &startEdge, unsigned &endEdge);
	__host__ __device__ void GenSubWinsForPseudoSrcFromPseudoSrc(const PseudoWindow &pseudoWin, unsigned &startEdge, unsigned &endEdge);

	__host__ __device__ bool IsValidWindow(const Window &win, bool isLeftChild);
	__host__ __device__ void BuildWindow(const Window &fatherWin,
		unsigned edge,
		double t0, double t1,
		const Vector2D &v0, const Vector2D &v1,
		Window &win);

	__host__ __device__ double Intersect(const Vector2D &v0, const Vector2D &v1, const Vector2D &p0, const Vector2D &p1);

public:
	// all these fields need to be allocated first (except for the PriorityQueues), 
	// and then to be assigned into
	Mesh *mesh;
	SplitInfo *splitInfos;
	VertInfo *vertInfos;
	PriorityQueues< Window > winQ;
	PriorityQueues< PseudoWindow > pseudoSrcQ;
	unsigned sourceVert;
	unsigned sourcePointFace; Vector3D sourcePointPos;

	Window *storedWindows; unsigned storedWindowsIdx; unsigned storedWindowsSize;
	unsigned *keptFaces; unsigned keptFacesIdx; unsigned keptFacesSize;

	bool pathPassVert;

	// statistics
	int numOfWinGen;
	int maxWinQSize, maxPseudoQSize;
	int totalCalcVertNum;
};

#endif