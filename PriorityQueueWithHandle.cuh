#ifndef GPU_PRIORITYQUEUEWITHHANDLE
#define GPU_PRIORITYQUEUEWITHHANDLE

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>

template <typename T>
class PriorityQueuesWithHandle
{
public:
	struct PQItem
	{
		double key;
		T item;
		int *backHandle;
	};

public:
	__host__ __device__ PriorityQueuesWithHandle() { maxSize = 0; tail = 1; };
	__host__ __device__ PriorityQueuesWithHandle(int maxSize_) { maxSize = maxSize_; tail = 1; };
	__host__ __device__ ~PriorityQueuesWithHandle();
	__host__ __device__ void AssignMemory(PQItem *d_pqs_, int maxSize_);

	__host__ __device__ void push(T item, int *backHandle, double key);
	__host__ __device__ T top();
	__host__ __device__ T pop();
	__host__ __device__ void decrease(int handle, double newKey);
	__host__ __device__ bool empty();
	__host__ __device__ void clear();
	__host__ __device__ int size();

private:
	PQItem *d_pqs;
	int tail;
	int maxSize;
};

template <typename T>
__host__ __device__ PriorityQueuesWithHandle<T>::~PriorityQueuesWithHandle()
{

}

template <typename T>
__host__ __device__ void PriorityQueuesWithHandle<T>::AssignMemory(PQItem *d_pqs_, int maxSize_)
{
	d_pqs = d_pqs_;
	maxSize = maxSize_;
}

template <typename T>
__host__ __device__ void PriorityQueuesWithHandle<T>::push(T item, int *backHandle, double key)
{
	if (tail - 1 >= maxSize) return;
	++tail;

	int curIdx = tail - 1;
	int pIdx = curIdx / 2;

	while (pIdx > 0 && d_pqs[pIdx].key > key)
	{
		d_pqs[curIdx] = d_pqs[pIdx];
		*(d_pqs[curIdx].backHandle) = curIdx;
		curIdx = pIdx;
		pIdx = curIdx / 2;
	}
	d_pqs[curIdx].item = item;
	d_pqs[curIdx].key = key;
	d_pqs[curIdx].backHandle = backHandle;
	*backHandle = curIdx;
}

template <typename T>
__host__ __device__ T PriorityQueuesWithHandle<T>::top()
{
	return d_pqs[1].item;
}

template <typename T>
__host__ __device__ T PriorityQueuesWithHandle<T>::pop()
{
	T ret = d_pqs[1].item;
	*(d_pqs[1].backHandle) = -1;
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
		*(d_pqs[curIdx].backHandle) = curIdx;

		curIdx = minChildIdx;
		leftChildIdx = curIdx * 2; rightChildIdx = curIdx * 2 + 1;
		if (rightChildIdx < tail)
			minChildIdx = d_pqs[leftChildIdx].key < d_pqs[rightChildIdx].key ? leftChildIdx : rightChildIdx;
		else if (leftChildIdx < tail)
			minChildIdx = leftChildIdx;
		else minChildIdx = -1;
	}

	d_pqs[curIdx] = topItem;
	*(d_pqs[curIdx].backHandle) = curIdx;

	return ret;
}

template <typename T>
__host__ __device__ void PriorityQueuesWithHandle<T>::decrease(int handle, double newKey)
{
	// TODO: decrease key
	int curIdx = handle;
	int pIdx = curIdx / 2;
	PQItem curItem = d_pqs[curIdx];

	while (pIdx > 0 && d_pqs[pIdx].key > newKey)
	{
		d_pqs[curIdx] = d_pqs[pIdx];
		*(d_pqs[curIdx].backHandle) = curIdx;
		curIdx = pIdx;
		pIdx = curIdx / 2;
	}

	d_pqs[curIdx] = curItem;
	d_pqs[curIdx].key = newKey;
	*(d_pqs[curIdx].backHandle) = curIdx;
}

template <typename T>
__host__ __device__ bool PriorityQueuesWithHandle<T>::empty()
{
	return tail == 1;
}

template <typename T>
__host__ __device__ void PriorityQueuesWithHandle<T>::clear()
{
	tail = 1;
}

template <typename T>
__host__ __device__ int PriorityQueuesWithHandle<T>::size()
{
	return tail - 1;
}

#endif