#ifndef MYHASHTABLE_CUH
#define MYHASHTABLE_CUH

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdio>

template <typename T>
class MyHashTable
{
public:
	struct HashItem
	{
		T item;
		int index;
	};

public:
	__host__ __device__ MyHashTable();
	__host__ __device__ MyHashTable(int size_, HashItem* data_);
	__host__ __device__ ~MyHashTable();

	// index operator overloading for writing
	__host__ __device__ T& operator[](unsigned index);

	// index operator only for reading
	__host__ __device__ const T& get(unsigned index);

	// This two methods are used for clearing the data
	__host__ __device__ unsigned Size();
	__host__ __device__ HashItem* Data();

	// used to get the hash key for specific index
	__host__ __device__ unsigned getKey(unsigned index);

private:
	__host__ __device__ unsigned hashFunc(unsigned index);

private:
	HashItem *data;

	unsigned size;
};

template <typename T>
__host__ __device__ MyHashTable<T>::MyHashTable()
{

}

template <typename T>
__host__ __device__ MyHashTable<T>::MyHashTable(int size_, HashItem* data_)
{
	size = size_;
	data = data_;
	for (unsigned i = 0; i < size; ++i)
		data[i].index = -1;
}

template <typename T>
__host__ __device__ MyHashTable<T>::~MyHashTable()
{
	data = NULL;
	size = 0;
}

template <typename T>
__host__ __device__ T& MyHashTable<T>::operator[](unsigned index)
{
	unsigned key = getKey(index);
// #ifndef __CUDA_ARCH__
// 	if (key == size)
// 		printf("Hash table is full when inserting element %d!\n", index);
// #endif
	data[key].index = index;
	return data[key].item;
}

template <typename T>
__host__ __device__ const T& MyHashTable<T>::get(unsigned index)
{
	unsigned key = getKey(index);
	return data[key].item;
}

template <typename T>
__host__ __device__ unsigned MyHashTable<T>::Size()
{
	return size;
}

template <typename T>
__host__ __device__ MyHashTable<T>::HashItem* MyHashTable<T>::Data()
{
	return data;
}

template <typename T>
__host__ __device__ unsigned MyHashTable<T>::getKey(unsigned index)
{
	// linear probing
	unsigned key = hashFunc(index);
	int i = 0;
	for (; i < size; ++i)
	{
		if (data[key].index == index || data[key].index == -1) break;
		key = (key + 1) % size;
	}
	if (i == size) key = size;
	return key;
}

template <typename T>
__host__ __device__ unsigned MyHashTable<T>::hashFunc(unsigned index)
{
	return static_cast<int>((index * 0.618033989 - static_cast<int>(index * 0.618033989)) * size);
}

#endif