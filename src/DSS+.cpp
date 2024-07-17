#include <bits/stdc++.h>
using namespace std;
const int INF = 2000000000;

struct Timer {
	double start_time, end_time;
	void start() { start_time = clock(); }
	void end() { end_time = clock(); }
	double time() { return (end_time - start_time) / CLOCKS_PER_SEC; }
};
Timer timer;
template <class T>
struct Set {
	T* nodes; bool* in; int size = -1;
	Set() {}
	Set(int sz) { size = 0; nodes = (T*)malloc(sz * sizeof(T)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); }
	void alloc(int sz) { size = 0; nodes = (T*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); }
	void insert(T x) { nodes[size++] = x; in[x] = true; }
	void clear() { for (int i = 0; i < size; i++) in[nodes[i]] = false; size = 0; }
	~Set() { free(nodes), free(in); }
};
template <class T>
struct Map {
	int* nodes; bool* in; int size = -1; T* value;
	Map() {}
	Map(int sz) { size = 0; nodes = (int*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); value = (T*)malloc(sz * sizeof(T)); memset(value, 0, sz * sizeof(T)); }
	void alloc(int sz) { size = 0; nodes = (int*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); value = (T*)malloc(sz * sizeof(T)); memset(value, 0, sz * sizeof(T)); }
	void freememory() { free(nodes), free(in), free(value); }
	void clear() { for (int i = 0; i < size; i++) in[nodes[i]] = false, value[nodes[i]] = 0; size = 0; }
	T& operator[](int x) { if (!in[x]) nodes[size++] = x, in[x] = true; return value[x]; }
	~Map() { free(nodes), free(in), free(value); }
};
template <class T>
struct Queue {
	T* nodes; int head, tail;
	Queue() {}
	Queue(int size) { nodes = (T*)malloc(size * sizeof(T)); head = tail = 0; }
	void alloc(int sz) { head = tail = 0; nodes = (int*)malloc(sz * sizeof(int)); }
	~Queue() { free(nodes); }
	bool empty() { return head == tail; }
	int pop() { return nodes[head++]; }
	void push(T x) { nodes[tail++] = x; }
	void clear() { head = tail = 0; }
};

inline int read_number(FILE* in) {
	int x = 0; char ch = 0; while (ch < '0' || ch > '9') ch = fgetc(in); while (ch >= '0' && ch <= '9') { x = x * 10 + (ch - '0'); ch = fgetc(in); } return x;
}
inline void check(bool flag, const char* message) {
	if (!flag) {
		printf("!!!!! CHECK ERROR !!!!!\n");
		printf("Error message: %s\n", message);
		assert(0);
	}
}

struct Edge { int u, v, to; };
struct Graph {
	int M, U, V, N, p;
	Edge* e;
	int* undeg, * indeg;
	int** adj;
	void read_graph_from_dataset(char* dataset_address);

	inline bool in_U(int x) { return x < U; }
	inline bool in_V(int x) { return x >= U; }
	inline bool in_S(int x) { return indeg[x] < (in_U(x) ? alpha : beta); }
	inline bool in_T(int x) { return indeg[x] > (in_U(x) ? alpha : beta); }

	void initialize_orientation();
	int alpha, beta;

	Set<int> D;
	void get_dense_subgraph();

	Set<int> S; Queue<int> Q;
	Map<int> parent; Set<int> vis; Map<int> dist, cur;

	bool DinicBFS();
	bool DinicDFS(int x);

	void check_correctness();
	const int OUTPUT_NODES_LIMIT = 10;
	void output_dense_subgraph();
	void display();
};
Graph G;

void Graph::read_graph_from_dataset(char* dataset_address) {
	FILE* in = fopen(dataset_address, "r");
	check(in != NULL, "Can not open file dataset_address\n");

	M = read_number(in); U = read_number(in); V = read_number(in);
	U++, V++; N = U + V;

	e = (Edge*)malloc(M * sizeof(Edge));
	undeg = (int*)malloc(N * sizeof(int)); indeg = (int*)malloc(N * sizeof(int));
	memset(undeg, 0, N * sizeof(int)); memset(indeg, 0, N * sizeof(int));
	adj = (int**)malloc(N * sizeof(int*));


	S.alloc(N); Q.alloc(N); parent.alloc(N); vis.alloc(N); D.alloc(N); dist.alloc(N); cur.alloc(N);

	for (int i = 0; i < M; i++) {
		e[i].u = read_number(in), e[i].v = read_number(in) + U;
		undeg[e[i].u]++, undeg[e[i].v]++;
	}
	for (int i = 0; i < N; i++) adj[i] = (int*)malloc(undeg[i] * sizeof(Edge));
	memset(undeg, 0, N * sizeof(int));
	for (int i = 0; i < M; i++) {
		int u = e[i].u, v = e[i].v;
		adj[u][undeg[u]++] = i;
		adj[v][undeg[v]++] = i;
	}
}
void Graph::initialize_orientation() {
	memset(indeg, 0, N * sizeof(int));
	for (int i = 0; i < M; i++) {
		// Orient to smaller undegree node
		int oriented_to;
		if (undeg[e[i].u] < undeg[e[i].v]) oriented_to = e[i].u;
		else oriented_to = e[i].v;
		e[i].to = oriented_to;
		indeg[oriented_to]++;
	}
}
void Graph::get_dense_subgraph() {
	while (DinicBFS()) {
		for (int x = 0; x < N; x++) {
			if (in_T(x))
				parent[x] = -2, cur[x] = 0, DinicDFS(x);
		}
	}
	D.clear(), Q.clear(), vis.clear();
	for (int x = 0; x < N; x++) {
		if (in_T(x))
			Q.push(x), D.insert(x);
	}
	while (!Q.empty()) {
		int x = Q.pop();
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int from = in_U(ne.to) ? ne.v : ne.u;
			if (D.in[from]) continue;
			Q.push(from), D.insert(from);
		}
	}
	return;
}
bool Graph::DinicBFS() {
	int dist_t = INF;

	Q.clear(), dist.clear(), parent.clear(), cur.clear();
	for (int x = 0; x < N; x++) {
		if (in_T(x))
			dist[x] = 1, Q.push(x);
	}

	bool break_loop = false;
	while (!Q.empty()) {
		int x = Q.pop();
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int from = in_U(ne.to) ? ne.v : ne.u;
			if (in_S(from)) {
				dist_t = dist[x] + 2; break_loop = true; break;
			}
			if (dist.in[from]) continue;
			dist[from] = dist[x] + 1;
			Q.push(from);
		}
		if (break_loop) break;
	}
	return dist_t != INF;
}
bool Graph::DinicDFS(int x) {
	if (in_S(x)) {
		indeg[x]++, indeg[e[parent[x]].to]--, e[parent[x]].to = x;
		return true;
	}
	for (int& j = cur[x]; j < undeg[x]; j++) {
		Edge& ne = e[adj[x][j]];
		if (ne.to != x) continue;
		int from = in_U(ne.to) ? ne.v : ne.u;
		if ((dist[from] != dist[x] + 1) && !in_S(from)) continue;
		parent[from] = adj[x][j];
		if (DinicDFS(from)) {
			if (parent[x] == -2) {
				if (indeg[x] == (in_U(x) ? alpha : beta)) return true;
				continue;
			}
			indeg[x]++, indeg[e[parent[x]].to]--, e[parent[x]].to = x;
			return true;
		}
	}
	return false;
}
void Graph::output_dense_subgraph() {
	printf("- %-20s: %d\n", "Dense subgraph size", D.size);

	if (true) {
		Set<int> output_nodes; output_nodes.alloc(N);
		sort(D.nodes, D.nodes + D.size);
		int i;
		for (i = 0; i < D.size && D.nodes[i] < U; i++)
			output_nodes.insert(D.nodes[i]);
		if (output_nodes.size <= OUTPUT_NODES_LIMIT) {
			printf("- %-20s: ", "Dense subgraph node in U");
			for (int j = 0; j < output_nodes.size; j++)
				printf("%d ", output_nodes.nodes[j]);
			printf("\n");
		}
		else {
			printf("- %-20s: %d\n", "Dense subgraph node in U are too many, with size of", output_nodes.size);
		}

		output_nodes.clear();
		for (; i < D.size; i++)
			output_nodes.insert(D.nodes[i]);
		if (output_nodes.size <= OUTPUT_NODES_LIMIT) {
			printf("- %-20s: ", "Dense subgraph node in V");
			for (int j = 0; j < output_nodes.size; j++)
				printf("%d ", output_nodes.nodes[j]);
			printf("\n");
		}
		else {
			printf("- %-20s: %d\n", "Dense subgraph node in V are too many, with size of", output_nodes.size);
		}
	}
}
void Graph::display() {
	for (int i = 0; i < M; i++) {
		Edge& ne = e[i];
		int from = in_U(ne.to) ? ne.v : ne.u;
		printf("%d %d\n", from, ne.to);
	}
}
void Graph::check_correctness() {
	/*
	* good: all good
	* point: Edges in E_\times(D, V \ D) are not all point to V \ D
	* Dlow: There is a node in D with indeg < alpha or beta
	* Dhigh: There is a node in V \ D with indeg > alpha or beta
	* Ddef: D does not comply with its definition
	* sumD: E(D) is not equal to \sum_{u\in D}indeg[u]
	* Cdef: C does not comply with its definition
	* sandwich1: D is not contained in C1
	* sandwich2: C2 is not contained in D
	*/
	int E1 = 0, E2 = 0;

	// point
	for (int i = 0; i < D.size; i++) {
		int x = D.nodes[i];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int y = in_U(x) ? ne.v : ne.u;
			check(D.in[y], "Edges in E_\\times(D, V \\ D) are not all point to V \\ D)");
		}
	}

	// Dlow and Dhigh
	for (int i = 0; i < N; i++) {
		check(!(D.in[i] && in_S(i)), "There is a node in D with indeg < alpha or beta");
		check(!(!D.in[i] && in_T(i)), "There is a node in V \\ D with indeg > alpha or beta");
	}

	// Ddef
	Q.clear(); vis.clear();
	for (int i = 0; i < N; i++) if (in_T(i)) Q.push(i), vis.insert(i);
	while (!Q.empty()) {
		int x = Q.pop();
		check(D.in[x], "D does not comply with its definition");
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int y = in_U(x) ? ne.v : ne.u;
			if (vis.in[y]) continue;
			Q.push(y); vis.insert(y);
		}
	}

	// sumD
	for (int i = 0; i < D.size; i++) {
		int x = D.nodes[i];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int y = in_U(x) ? ne.v : ne.u;
			if (D.in[y])
				E1++;
		}
	}
	for (int i = 0; i < D.size; i++) E2 += indeg[D.nodes[i]];
	check(E1 == E2, "E(D) is not equal to \\sum_{u\\in D}indeg[u]");

	return;
}

int main(int argc, char** argv) {
	check(argc == 2, "The number of arguments of main are incorrect. ");
	char dataset_address[1000]; strcpy(dataset_address, argv[1]);

	double runtime;
	printf("----------Now processing %s----------\n", dataset_address);
	printf("- %-20s: %s\n", "Algorithm used", "DSS+");

	timer.start();
	G.read_graph_from_dataset(dataset_address);
	timer.end(); runtime = timer.time();
	printf("- %-20s: %d, %d, %d\n", "|E|, |U|, |V|", G.M, G.U, G.V);
	printf("- %-20s: %lf\n", "Read graph time", runtime);

	timer.start();
	G.initialize_orientation();
	timer.end(); runtime = timer.time();
	printf("- %-20s: %lf\n", "Initialization time", runtime);

	vector<int> test_alpha, test_beta; int test_pointer = 0;

	vector<double> runtime_list;
	while (true) {
		printf("----------------------------------\n");
		printf("> Enter alpha and beta, -1 -1 for quit\n");
		if (test_alpha.size() == 0) {
			scanf("%d%d", &G.alpha, &G.beta);
		}
		else {
			G.alpha = test_alpha[test_pointer], G.beta = test_beta[test_pointer];
			test_pointer++;
		}
		if (G.alpha == -1 || G.beta == -1)
			break;
		printf("- %-16s: %d, %d\n", "alpha, beta", G.alpha, G.beta);

		timer.start();
		G.initialize_orientation();
		G.get_dense_subgraph();
		timer.end(); runtime = timer.time();
		printf("- %-20s: %lf\n", "Get D_{a,b} time", runtime); runtime_list.push_back(runtime);
		G.output_dense_subgraph();

		G.check_correctness();
	}
	return 0;
}