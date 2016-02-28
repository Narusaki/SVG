#ifndef GPU_PRIORITYQUEUE
#define GPU_PRIORITYQUEUE

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>

template <typename T>
class PriorityQueues
{
public:
	struct PQItem
	{
		double key;
		T item;
	};

public:
	__host__ __device__ PriorityQueues() { maxSize = 0; tail = 1; };
	__host__ __device__ PriorityQueues(int maxSize_) { maxSize = maxSize_; tail = 1; };
	__host__ __device__ ~PriorityQueues();
	__host__ __device__ void AssignMemory(PQItem *d_pqs_, int maxSize_);

	__host__ __device__ void push(T item, double key);
	__host__ __device__ T top();
	__host__ __device__ T pop();
	__host__ __device__ bool empty();
	__host__ __device__ void clear();
	__host__ __device__ int size();

private:
	PQItem *d_pqs;
	int tail;
	int maxSize;
};

template <typename T>
__host__ __device__ PriorityQueues<T>::~PriorityQueues()
{

}

template <typename T>
__host__ __device__ void PriorityQueues<T>::AssignMemory(PQItem *d_pqs_, int maxSize_)
{
	d_pqs = d_pqs_;
	maxSize = maxSize_;
}

template <typename T>
__host__ __device__ void PriorityQueues<T>::push(T item, double key)
{
	if (tail - 1 >= maxSize) return;
	++tail;

	int curIdx = tail - 1;
	int pIdx = curIdx / 2;

	while (pIdx > 0 && d_pqs[pIdx].key > key)
	{
		d_pqs[curIdx] = d_pqs[pIdx];
		curIdx = pIdx;
		pIdx = curIdx / 2;
	}
	d_pqs[curIdx].item = item;
	d_pqs[curIdx].key = key;
}

template <typename T>
__host__ __device__ T PriorityQueues<T>::top()
{
	return d_pqs[1].item;
}

template <typename T>
__host__ __device__ T PriorityQueues<T>::pop()
{
	T ret = d_pqs[1].item;
	--tail;
	PQItem topItem = d_pqs[tail];

	int curIdx = 1;
	int leftChildIdx = curIdx * 2, rightChildIdx = curIdx * 2 + 1;
	int minChildIdx = -1;
	if (rightChildIdx < tail)
		minChildIdx = d_pqs[leftChildIdx].key < d_pqs[rightChildIdx].key ? leftChildIdx : rightChildIdx;
	else if (leftChildIdx < tail)
		minChildIdx = leftChildIdx;

	while (minChildIdx != -1 && d_pqs[minChildIdx].key < topItem.key)
	{
		d_pqs[curIdx] = d_pqs[minChildIdx];
		curIdx = minChildIdx;
		leftChildIdx = curIdx * 2; rightChildIdx = curIdx * 2 + 1;
		if (rightChildIdx < tail)
			minChildIdx = d_pqs[leftChildIdx].key < d_pqs[rightChildIdx].key ? leftChildIdx : rightChildIdx;
		else if (leftChildIdx < tail)
			minChildIdx = leftChildIdx;
		else minChildIdx = -1;
	}

	d_pqs[curIdx] = topItem;

	return ret;
}

template <typename T>
__host__ __device__ bool PriorityQueues<T>::empty()
{
	return tail == 1;
}

template <typename T>
__host__ __device__ void PriorityQueues<T>::clear()
{
	tail = 1;
}

template <typename T>
__host__ __device__ int PriorityQueues<T>::size()
{
	return tail - 1;
}

#endif