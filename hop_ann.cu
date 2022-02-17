#include<stdio.h>
#include<string.h>
#include<algorithm>
#include<queue>
#include<iostream>

int V, D, E, L, K, A, B, C, M, Q;
int* X_d;
int* X;
int* edges;

int squared_l2_dist(int* x, int* y, int D) {
	int sum2 = 0;
	for (int i = 0; i < D; ++i)
		sum2 += (x[i] - y[i]) * (x[i] - y[i]);
	return sum2;
}

__global__ void computeParallel(int* X, int* query, int D, int* hop, int* id, int* d) {
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int sum2 = 0;
	int* x = X + id[index] * D;

	int* y = query;
	for (int i = 0; i < D; ++i) {
		sum2 += (x[i] - y[i]) * (x[i] - y[i]);
	}

	d[index] = sum2;
}

int nearest_id(int start_point, int max_hop, int* query_data, int* X) {

	//std::cout<<"Edges : ";
	//for (int i = 0; i<V*(L+1); ++i)
	//	std::cout<<edges[i];
	//std::cout<<"\n";

	std::queue<std::pair<int, int>> q;
	q.push(std::make_pair(start_point, 0));
	int min_d = std::numeric_limits<int>::max();
	int min_id = -1;
	do {
		int* id;
		int* hop;
		int siz = q.size();
		//std::cout<<"Size : "<<siz<<"\n";
		cudaMallocManaged(&id, siz * sizeof(int));
		cudaMallocManaged(&hop, siz * sizeof(int));
		int ctr = 0;
		while (!q.empty()) {
			auto now = q.front();
			q.pop();
			id[ctr] = now.first;
			hop[ctr] = now.second;

			//std::cout<<"ID : "<<id[ctr]<<" "<<"HOP : "<<hop[ctr]<<"\n";
			ctr++;
		}
		int* d_d;
		cudaMallocManaged(&d_d, siz * sizeof(int));
		computeParallel << <siz, 1 >> > (X, query_data, D, hop, id, d_d);
		//cudaDeviceSynchronize();

		int* d = new int[siz];
		cudaMemcpy(d, d_d, siz * sizeof(int), cudaMemcpyDeviceToHost);

		for (int i = 0; i < siz; ++i) {
			//std::cout<<"Node rn : "<<id[i]<<" " <<"Hop : "<<hop[i]<<"\n";
			if ((d[i] < min_d) || (d[i] == min_d && id[i] < min_id)) {
				min_d = d[i];
				min_id = id[i];
			}
			if (hop[i] + 1 <= max_hop) {
				//std::cout<<"Nodes being inserted for id : "<<id[i]<< " : ";
				int degree = edges[id[i] * (L + 1)];
				//std::cout<<"Degree : "<<degree<<" : ";
				for (int j = 1; j <= degree; ++j) {
					int v = edges[id[i] * (L + 1) + j];
					//std::cout<<"ID : "<< id[i] << " L+1 : "<<L+1<<" Inter : "<<id[i]*(L+1)<<" Index : "<<id[i] * (L + 1) + i<<" Value : "<<v<<" ";
					q.push(std::make_pair(v, hop[i] + 1));
				}
				//std::cout<<"\n";
			}
		}
		//std::cout<<"\n";
		cudaFree(id);
		cudaFree(hop);
		cudaFree(d);
	} while (!q.empty());
	//std::cout<<"\n\n";
	return min_id;
}

int main(int argc, char** argv) {
	FILE* fin = fopen(argv[1], "r");
	FILE* fout = fopen(argv[2], "w");
	fscanf(fin, "%d%d%d%d%d%d%d%d%d%d", &V, &D, &E, &L, &K, &A, &B, &C, &M, &Q);
	X = new int[V * D];



	for (int i = 0; i < K; ++i)
		fscanf(fin, "%d", &X[i]);

	for (int i = K; i < V * D; ++i)
	{
		X[i] = ((long long)A * X[i - 1] + (long long)B * X[i - 2] + C) % M;
	}


	edges = new int[V * (L + 1)];
	for (int i = 0; i < V; ++i) {
		edges[i * (L + 1)] = 0;
	}
	for (int i = 0; i < E; ++i) {
		int u, v;
		fscanf(fin, "%d%d", &u, &v);
		int degree = edges[u * (L + 1)];
		edges[u * (L + 1) + degree + 1] = v;
		++edges[u * (L + 1)];
	}
	std::cout << "Edges: ";

	for (int i = 0; i < V * (L + 1); ++i) {
		std::cout << edges[i] << " ";
	}
	cudaMallocManaged(&X_d, (V * D) * sizeof(int));
	cudaMemcpy(X_d, X, V * D * sizeof(int), cudaMemcpyHostToDevice);

	//int* query_data = new int[D];
	//int* query_data_d;
	//cudaMallocManaged(&query_data_d, D * sizeof(int));
	//for (int i = 0; i < Q; ++i) {
	//	int start_point, hop;
	//	fscanf(fin, "%d%d", &start_point, &hop);
	//	for (int i = 0; i < D; ++i) {
	//		fscanf(fin, "%d", &query_data[i]);
	//	}
	//	cudaMemcpy(query_data_d, query_data, D * sizeof(int), cudaMemcpyHostToDevice);
	//	fprintf(fout, "%d\n", nearest_id(start_point, hop, query_data_d, X_d));
	//}


	int* query_data = new int[D];

	cudaMallocManaged(&query_data, D * sizeof(int));
	for (int i = 0; i < Q; ++i) {
		int start_point, hop;
		fscanf(fin, "%d%d", &start_point, &hop);
		for (int i = 0; i < D; ++i) {
			fscanf(fin, "%d", &query_data[i]);
		}
		fprintf(fout, "%d\n", nearest_id(start_point, hop, query_data_d, X_d));
	}


	fclose(fin);
	fclose(fout);
	//delete[] X;
	delete[] edges;
	//delete[] query_data;
	cudaFree(X);
	cudaFree(query_data);
	return 0;
}

