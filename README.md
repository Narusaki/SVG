# Saddle Vertex Graph
This is an implementation of Saddle Vertex Graph (SVG) (proposed in [1]), which aims at solving discrete geodesic problem on triangular meshes with large complexities.

The program here provided is implemented on CUDA C and has been tested on Nvidia Tesla 40c. Following is the testing results on some typical meshes (each (t1, t2) pair represents the constructing and solving time, in seconds).

| Mesh          | # of vertices | K = 100            | K = 500        |
| ------------- | ------------- | ------------------ | -------------- |
| Fertility     | 30K           | (7.205, 0.035)     | (57.16, 0.108) |
| Stanford Bunny| 72K           | (14.826, 0.064)    | (99.597, 0.195)|
| Dragon        | 750K          | (1237.63, 0.69)    | Out-of-memory  |

## Note:
1. The time of constructing SVG is influenced by many factors, including mesh complexity, number of invoked threads, efficiency of ICH sub-procedure etc.
2. The "solving" part of the implementation here does NOT follow the paper strictly, i.e. using a naive A-star and Dijkstra implementation separately for SSSD and MSAD problems.

## Reference
[1] Ying X, Wang X, He Y. Saddle vertex graph (SVG): a novel solution to the discrete geodesic problem[J]. ACM Transactions on Graphics (TOG), 2013, 32(6): 170.
