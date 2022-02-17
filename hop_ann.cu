#include<stdio.h>
#include<string.h>
#include<algorithm>
#include<queue>

using namespace std;

int V,D,E,L,K,A,B,C,M,Q;
int* X;
int* edges;

int squared_l2_dist(int* x,int* y,int D){
	int sum2 = 0;
	for(int i = 0;i < D;++i)
		sum2 += (x[i] - y[i]) * (x[i] - y[i]);
	return sum2;
}


//__global__ void cuda_squared_l2_dist(int* origin, int* nodes, int* distances) {
//
//	int index = threadIdx.x + blockDim.x * blockIdx.x;
//	int* x = nodes[index];
//
//	int sum2 = 0;
//	for (int i = 0; i < D; ++i)
//		sum2 += (origin[i] - x[i]) * (origin[i] - x[i]);
//
//	distances[index] = sum2;
//}


vector<int> explore(int start_point, int max_hop) {
	queue<pair<int, int>> q;
	vector<int> nodes;
	q.push(make_pair(start_point, 0));
	nodes.push_back(start_point);
	while (!q.empty()) {
		auto now = q.front();
		q.pop();
		int id = now.first;
		int hop = now.second;

		if (hop + 1 <= max_hop) {
			int degree = edges[id * (L + 1)];

			for (int i = 1; i <= degree; ++i) {
				int v = edges[id * (L + 1) + i];

				q.push(std::make_pair(v, hop + 1));
				nodes.push_back(v);
			}

		}
	}

	return nodes;
}

int main(int argc,char** argv){
	FILE* fin = fopen(argv[1],"r");
	FILE* fout = fopen(argv[2],"w");
	fscanf(fin,"%d%d%d%d%d%d%d%d%d%d",&V,&D,&E,&L,&K,&A,&B,&C,&M,&Q);
	X = new int[V * D];
	for(int i = 0;i < K;++i)
		fscanf(fin,"%d",&X[i]);
	for(int i = K;i < V * D;++i)
		X[i] = ((long long)A * X[i - 1] + (long long)B * X[i - 2] + C) % M;
	edges = new int[V * (L + 1)];
	for(int i = 0;i < V;++i){
		edges[i * (L + 1)] = 0;
	}
	for(int i = 0;i < E;++i){
		int u,v;
		fscanf(fin,"%d%d",&u,&v);
		int degree = edges[u * (L + 1)];
		edges[u * (L + 1) + degree + 1] = v;
		++edges[u * (L + 1)];
	}


	int* query_data = new int[D];

	// cuda 
	//cudaMallocManaged(&query_data, D * sizeof(int));
	
	for(int i = 0;i < Q;++i){
		int start_point,hop;
		fscanf(fin,"%d%d",&start_point,&hop);
		for(int i = 0;i < D;++i){
			fscanf(fin,"%d",&query_data[i]);
		}

		// explore for all nodes 
		vector<int> allPossibleNodes = explore(start_point, hop);

		

		// non cuda component 
		// 
		int* distances = new int[allPossibleNodes.size()];

		int* targets = new int[allPossibleNodes.size()*D];
		
		//targets = new int*[allPossibleNodes.size()];
		//for (int j = 0; j < allPossibleNodes.size(); ++j) {
		//	targets[j] = new int[1]; 
		//}


		// cuda 
		//cudaMallocManaged(&targets, D * allPossibleNodes.size() * sizeof(int));
		//cudaMallocManaged(&distances, allPossibleNodes.size() * sizeof(int));
		// cuda 


		for (int j = 0; j < allPossibleNodes.size(); ++j) {
			int p = j * D;

			int* temp = X + allPossibleNodes.at(j) * D;

			for (int k = 0; k < D; k++) {
				targets[p + k] = temp[k];
			}
		}



		// non cuda
		
		for (int j = 0; j < allPossibleNodes.size(); ++j) {
			distances[j] = squared_l2_dist(targets + j * D, query_data, D);
		}

		// non cuda 
		

		//cuda 

		//int threadsPerBlock = 256;
		//int blocksPerGrid = (allPossibleNodes.size() + threadsPerBlock - 1) / threadsPerBlock;

		//cuda_squared_l2_dist << <blocksPerGrid, threadsPerBlock > >> (query_data, targets, distances);
		//cudaDeviceSynchronize();
		// 
		// 
		//cuda 

		// get min 
		int min_d = 2147483647;
		int min_id = 2147483647;
		for (int j = 0; j < allPossibleNodes.size(); ++j) {
			int id = allPossibleNodes.at(j);
			int d = distances[j];

			if (d < min_d || (d == min_d && id < min_id)) {
				min_d = d;
				min_id = id;
			}
		}

		fprintf(fout,"%d\n",min_id);
	}
	fclose(fin);
	fclose(fout);

	delete[] X;
	delete[] edges;
	delete[] query_data;

	return 0;
}

