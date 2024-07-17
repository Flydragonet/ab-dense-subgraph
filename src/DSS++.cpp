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

	Set<int> C1, C2, C1_minus_C2;
	void get_core();

	Set<int> D;
	void get_dense_subgraph_from_core();

	Set<int> S; Queue<int> Q;
	Map<int> parent; Set<int> vis; Map<int> dist, cur;

	bool DinicBFS();
	bool DinicDFS(int x);

	void check_correctness();
	const int OUTPUT_NODES_LIMIT = 10;
	void output_core(Set<int>& C);
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


	S.alloc(N); Q.alloc(N); parent.alloc(N); vis.alloc(N); D.alloc(N); C1.alloc(N); C2.alloc(N); C1_minus_C2.alloc(N); dist.alloc(N); cur.alloc(N);

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
	for (int i = 0; i < C1_minus_C2.size; i++) {
		indeg[C1_minus_C2.nodes[i]] = 0;
	}
	for (int i = 0; i < C1_minus_C2.size; i++) {
		int x = C1_minus_C2.nodes[i];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			int y = in_U(x) ? ne.v : ne.u;
			if (!C1.in[y])
				ne.to = y;
			else if (C2.in[y]) {
				indeg[x]++; ne.to = x;
			}
			else if (in_U(x)) {
				if (undeg[x] < undeg[y]) ne.to = x, indeg[x]++;
				else ne.to = y, indeg[y]++;
			}
		}
	}
}
void Graph::get_core() {
	int core1_alpha = alpha + 1, core1_beta = beta + 1, core2_alpha = 2 * alpha + 1, core2_beta = 2 * beta + 1;
	C1.clear(), C2.clear();

	// a node --> will_be_deleted_nodes --> deleted_nodes
	Set<int> will_be_deleted_nodes; will_be_deleted_nodes.alloc(N); int pointer = 0;
	Set<int> deleted_nodes; deleted_nodes.alloc(N);
	int* temporary_undeg = (int*)malloc(N * sizeof(int)); memcpy(temporary_undeg, undeg, N * sizeof(int));

	for (int x = 0; x < N; x++)
		if (undeg[x] < (in_U(x) ? core1_alpha : core1_beta))
			will_be_deleted_nodes.insert(x);
	while (pointer < will_be_deleted_nodes.size) {
		int x = will_be_deleted_nodes.nodes[pointer++];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			int y = in_U(x) ? ne.v : ne.u;
			if (deleted_nodes.in[y]) continue;
			if (temporary_undeg[y] == (in_U(y) ? core1_alpha : core1_beta)) {
				will_be_deleted_nodes.insert(y);
			}
			temporary_undeg[x]--;
			temporary_undeg[y]--;
		}
		deleted_nodes.insert(x);
	}
	for (int x = 0; x < N; x++)
		if (!deleted_nodes.in[x])
			C1.insert(x);

	for (int i = 0; i < C1.size; i++) {
		int x = C1.nodes[i];
		if (temporary_undeg[x] < (in_U(x) ? core2_alpha : core2_beta))
			will_be_deleted_nodes.insert(x);
	}
	while (pointer < will_be_deleted_nodes.size) {
		int x = will_be_deleted_nodes.nodes[pointer++];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			int y = in_U(x) ? ne.v : ne.u;
			if (deleted_nodes.in[y]) continue;
			if (temporary_undeg[y] == (in_U(y) ? core2_alpha : core2_beta)) {
				will_be_deleted_nodes.insert(y);
			}
			temporary_undeg[x]--;
			temporary_undeg[y]--;
		}
		deleted_nodes.insert(x);
	}
	for (int i = 0; i < C1.size; i++) {
		int x = C1.nodes[i];
		if (!deleted_nodes.in[x])
			C2.insert(x);
	}

	free(temporary_undeg);
	return;
}
void Graph::get_dense_subgraph_from_core() {
	C1_minus_C2.clear();
	for (int i = 0; i < C1.size; i++) {
		int x = C1.nodes[i];
		if (!C2.in[x])
			C1_minus_C2.insert(x);
	}

	initialize_orientation();

	while (DinicBFS()) {
		for (int i = 0; i < C1_minus_C2.size; i++) {
			int x = C1_minus_C2.nodes[i];
			if (in_T(x))
				parent[x] = -2, cur[x] = 0, DinicDFS(x);
		}
	}

	D.clear(), Q.clear(), vis.clear();
	for (int i = 0; i < C1_minus_C2.size; i++) {
		int x = C1_minus_C2.nodes[i];
		if (in_T(x))
			Q.push(x), D.insert(x);
	}
	while (!Q.empty()) {
		int x = Q.pop();
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int from = in_U(ne.to) ? ne.v : ne.u;
			check(C1.in[from], "C1_minus_C2 error");
			if (C2.in[from]) continue;
			if (D.in[from]) continue;
			Q.push(from), D.insert(from);
		}
	}
	for (int i = 0; i < C2.size; i++) {
		int x = C2.nodes[i];
		check(!D.in[x], "C2 error");
		D.insert(x);
	}
	return;
}
bool Graph::DinicBFS() {
	int dist_t = INF;

	Q.clear(), dist.clear(), parent.clear(), cur.clear();
	for (int i = 0; i < C1_minus_C2.size; i++) {
		int x = C1_minus_C2.nodes[i];
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
			check(C1.in[from], "C1_minus_C2 error");
			if (C2.in[from]) continue;
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
		check(C1.in[from], "C1_minus_C2 error");
		if (C2.in[from]) continue;
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
void Graph::output_core(Set<int>& C) {
	printf("- %-20s: %d\n", "Core size", C.size);

	if (true) {
		Set<int> output_nodes; output_nodes.alloc(N);
		sort(C.nodes, C.nodes + C.size);
		int i;
		for (i = 0; i < C.size && C.nodes[i] < U; i++)
			output_nodes.insert(C.nodes[i]);
		if (output_nodes.size <= OUTPUT_NODES_LIMIT) {
			printf("- %-20s: ", "Core node in U");
			for (int j = 0; j < output_nodes.size; j++)
				printf("%d ", output_nodes.nodes[j]);
			printf("\n");
		}
		else {
			printf("- %-20s: %d\n", "Core node in U are too many, with size of", output_nodes.size);
		}

		output_nodes.clear();
		for (; i < C.size; i++)
			output_nodes.insert(C.nodes[i]);
		if (output_nodes.size <= OUTPUT_NODES_LIMIT) {
			printf("- %-20s: ", "Core node in V");
			for (int j = 0; j < output_nodes.size; j++)
				printf("%d ", output_nodes.nodes[j]);
			printf("\n");
		}
		else {
			printf("- %-20s: %d\n", "Core node in V are too many, with size of", output_nodes.size);
		}
	}
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
		if (!C1_minus_C2.in[x])
			continue;
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int y = in_U(x) ? ne.v : ne.u;
			check(D.in[y], "Edges in E_\\times(D, V \\ D) are not all point to V \\ D)");
		}
	}

	// Dlow and Dhigh
	/*
	for (int i = 0; i < N; i++) {
		check(!(D.in[i] && in_S(i)), "There is a node in D with indeg < alpha or beta");
		check(!(!D.in[i] && in_T(i)), "There is a node in V \\ D with indeg > alpha or beta");
	}*/
	for (int i = 0; i < C1_minus_C2.size; i++) {
		int x = C1_minus_C2.nodes[i];
		check(!(D.in[x] && in_S(x)), "There is a node in D with indeg < alpha or beta");
		check(!(!D.in[x] && in_T(x)), "There is a node in V \\ D with indeg > alpha or beta");
	}

	// Ddef
	/*
	Q.clear(); vis.clear();
	for (int i = 0; i < N; i++) if (in_T(i)) Q.push(i), vis[i] = true;
	while (!Q.empty()) {
		int x = Q.pop();
		check(D.in[x], "D does not comply with its definition");
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int y = in_U(x) ? ne.v : ne.u;
			if (vis[y]) continue;
			Q.push(y); vis[y] = true;
		}
	}*/

	// sumD
	/*
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
	check(E1 == E2, "E(D) is not equal to \\sum_{u\\in D}indeg[u]");*/

	// Cdef
	for (int i = 0; i < C1.size; i++) {
		int x = C1.nodes[i], count_undeg = 0;
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			int y = in_U(x) ? ne.v : ne.u;
			if (C1.in[y])
				count_undeg++;
		}
		check(count_undeg >= (in_U(x) ? alpha + 1 : beta + 1), "C1 does not comply with its definition");
	}
	for (int i = 0; i < C2.size; i++) {
		int x = C2.nodes[i], count_undeg = 0;
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			int y = in_U(x) ? ne.v : ne.u;
			if (C2.in[y])
				count_undeg++;
		}
		check(count_undeg >= (in_U(x) ? 2 * alpha + 1 : 2 * beta + 1), "C2 does not comply with its definition");
	}

	// sandwich1
	for (int i = 0; i < D.size; i++) {
		int x = D.nodes[i];
		check(C1.in[x], "D is not contained in C");
	}

	// sandwich2
	for (int i = 0; i < C2.size; i++) {
		int x = C2.nodes[i];
		check(D.in[x], "C2 is not contained in D");
	}

	return;
}

int main(int argc, char** argv) {
	check(argc == 2, "The number of arguments of main are incorrect. ");
	char dataset_address[1000]; strcpy(dataset_address, argv[1]);

	double runtime;
	printf("----------Now processing %s----------\n", dataset_address);
	printf("- %-20s: %s\n", "Algorithm used", "DSS++");

	timer.start();
	G.read_graph_from_dataset(dataset_address);
	timer.end(); runtime = timer.time();
	printf("- %-20s: %d, %d, %d\n", "|E|, |U|, |V|", G.M, G.U, G.V);
	printf("- %-20s: %lf\n", "Read graph time", runtime);

	vector<int> test_alpha, test_beta; int test_pointer = 0;

	vector<double> runtime_list, runtime_list2;
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
		G.get_core();
		timer.end(); runtime = timer.time();
		printf("- %-20s: %lf\n", "Get C_1 C_2 time", runtime); runtime_list2.push_back(runtime);
		G.output_core(G.C1);
		G.output_core(G.C2);

		timer.start();
		G.get_dense_subgraph_from_core();
		timer.end(); runtime = timer.time();
		printf("- %-20s: %lf\n", "Get D_{a,b} time", runtime); runtime_list.push_back(runtime);
		G.output_dense_subgraph();

		G.check_correctness();
	}
	return 0;
}