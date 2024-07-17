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
	int M, U, V, N;
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

	Set<int> C1, C2;
	Set<int> C1_minus_C2;
	void get_core(int core_alpha, int core_beta, Set<int>& C);

	int* sorted; Set<int> D;
	bool* computed; int* position, * number_of_edges, * number_of_nodes;
	void change_sorted(int position_l, int position_u, Set<int>& D);

	int alpha_max, beta_max, d_max;
	void get_dense_subgraph_from_core();
	void ReTest(int _alpha, int _beta, bool the_first_iteration);
	void Divide(int D_l_rank, int D_u_rank, bool the_first_iteration);

	int* r;

	Set<int> S; Queue<int> Q;
	Map<int> parent; Set<int> vis; Map<int> dist, cur;

	bool DinicBFS(int D_l_rank, int D_u_rank);
	bool DinicDFS(int x, int D_l_rank, int D_u_rank);

	void analyze_r();
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
	r = (int*)malloc(N * sizeof(int));

	S.alloc(N); Q.alloc(N); parent.alloc(N); vis.alloc(N); C1.alloc(N); C2.alloc(N); C1_minus_C2.alloc(N); dist.alloc(N); cur.alloc(N); D.alloc(N);

	for (int i = 0; i < M; i++) {
		e[i].u = read_number(in), e[i].v = read_number(in) + U;
		undeg[e[i].u]++, undeg[e[i].v]++;
	}
	d_max = alpha_max = beta_max = 0;
	for (int u = 0; u < U; u++) alpha_max = max(alpha_max, undeg[u]), adj[u] = (int*)malloc(undeg[u] * sizeof(int));
	for (int v = U; v < N; v++) beta_max = max(beta_max, undeg[v]), adj[v] = (int*)malloc(undeg[v] * sizeof(int));
	d_max = max(alpha_max, beta_max) + 2;
	sorted = (int*)malloc(N * sizeof(int)), computed = (bool*)malloc(d_max * sizeof(bool)), position = (int*)malloc(d_max * sizeof(int)), number_of_edges = (int*)malloc(d_max * sizeof(int)), number_of_nodes = (int*)malloc(d_max * sizeof(int));

	memset(undeg, 0, N * sizeof(int));
	for (int i = 0; i < M; i++) {
		int u = e[i].u, v = e[i].v;
		adj[u][undeg[u]++] = i;
		adj[v][undeg[v]++] = i;
	}
}
void Graph::initialize_orientation() {
	for (int i = 0; i < M; i++) {
		// Arbitrarily orient
		int oriented_to;
		if (i % 2 == 0) oriented_to = e[i].u;
		else oriented_to = e[i].v;
		e[i].to = oriented_to;
		indeg[oriented_to]++;
	}
}
void Graph::get_core(int core_alpha, int core_beta, Set<int>& C) {
	C.clear();

	// a node --> will_be_deleted_nodes --> deleted_nodes
	Set<int> will_be_deleted_nodes; will_be_deleted_nodes.alloc(N); int pointer = 0;
	Set<int> deleted_nodes; deleted_nodes.alloc(N);
	int* temporary_undeg = (int*)malloc(N * sizeof(int)); memcpy(temporary_undeg, undeg, N * sizeof(int));

	for (int x = 0; x < N; x++)
		if (undeg[x] < (in_U(x) ? core_alpha : core_beta))
			will_be_deleted_nodes.insert(x);

	while (pointer < will_be_deleted_nodes.size) {
		int x = will_be_deleted_nodes.nodes[pointer++];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			int y = in_U(x) ? ne.v : ne.u;
			if (deleted_nodes.in[y]) continue;
			if (temporary_undeg[y] == (in_U(y) ? core_alpha : core_beta)) {
				will_be_deleted_nodes.insert(y);
			}
			temporary_undeg[x]--;
			temporary_undeg[y]--;
		}
		deleted_nodes.insert(x);
	}

	for (int x = 0; x < N; x++)
		if (!deleted_nodes.in[x])
			C.insert(x);

	free(temporary_undeg);
	return;
}
int max_rank;
void Graph::change_sorted(int position_l, int position_u, Set<int>& D) {
	position_u--;
	while (position_l <= position_u) {
		while (position_l <= position_u && !D.in[sorted[position_l]]) position_l++;
		while (position_l <= position_u && D.in[sorted[position_u]]) position_u--;
		if (position_l <= position_u)
			swap(sorted[position_l], sorted[position_u]);
	}
}
void Graph::Divide(int D_l_rank, int D_u_rank, bool the_first_iteration) {
	// printf("LOG: enter Dividea(%d, %d)\n", D_l_rank, D_u_rank);
	int E_all = number_of_edges[D_l_rank] - number_of_edges[D_u_rank], l = D_l_rank, u = D_u_rank;
	while (u > l) {
		int mid = (u + l + 1) / 2;
		if (the_first_iteration) ReTest(alpha, mid, true);
		else ReTest(mid, beta, false);
		if (number_of_edges[D_l_rank] - number_of_edges[mid] < E_all / 2)
			l = mid;
		else
			u = mid - 1;
	}
	int m = l;
	if (m - D_l_rank >= 2 && number_of_nodes[D_l_rank] != number_of_nodes[m])
		Divide(D_l_rank, m, the_first_iteration);
	m++;
	if (D_u_rank - m >= 2 && number_of_nodes[D_u_rank] != number_of_nodes[m])
		Divide(m, D_u_rank, the_first_iteration);
	return;
}
void Graph::ReTest(int _alpha, int _beta, bool the_first_iteration) {
	alpha = _alpha, beta = _beta;
	// printf("LOG: ReTest(%d, %d)\n", alpha, beta);
	int D_l_rank, D_u_rank;
	if (the_first_iteration) {
		if (computed[beta]) return;
		D_l_rank = beta - 1, D_u_rank = beta + 1;
		while (!computed[D_l_rank]) D_l_rank--;
		while (!computed[D_u_rank]) D_u_rank++;
	}
	else {
		if (computed[alpha]) return;
		D_l_rank = alpha - 1, D_u_rank = alpha + 1;
		while (!computed[D_l_rank]) D_l_rank--;
		while (!computed[D_u_rank]) D_u_rank++;
	}

	C1_minus_C2.clear(); int D_l_position = position[D_l_rank], D_u_position = position[D_u_rank];
	for (int i = D_l_position; i < D_u_position; i++) {
		C1_minus_C2.insert(sorted[i]);
	}

	while (DinicBFS(D_l_rank, D_u_rank)) {
		for (int i = 0; i < C1_minus_C2.size; i++) {
			int x = C1_minus_C2.nodes[i];
			if (in_T(x))
				parent[x] = -2, cur[x] = 0, DinicDFS(x, D_l_rank, D_u_rank);
		}
	}

	int index = the_first_iteration ? beta : alpha;
	Set<int>& ans = D;
	ans.clear(), Q.clear(), vis.clear();
	for (int i = 0; i < C1_minus_C2.size; i++) {
		int x = C1_minus_C2.nodes[i];
		if (in_T(x))
			Q.push(x), ans.insert(x);
	}
	while (!Q.empty()) {
		int x = Q.pop();
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int from = in_U(ne.to) ? ne.v : ne.u;
			check(r[from] >= D_l_rank, "r[from] < index");
			if (r[from] >= D_u_rank) continue;
			if (ans.in[from]) continue;
			Q.push(from), ans.insert(from);
		}
	}
	for (int i = 0; i < ans.size; i++) {
		int x = ans.nodes[i];
		check(r[x] == D_l_rank, "r[x] != D_l_rank");
		r[x] = index;
	}
	computed[index] = true;
	number_of_nodes[index] = ans.size + number_of_nodes[D_u_rank]; position[index] = D_u_position - ans.size;
	change_sorted(D_l_position, D_u_position, ans);
	int& ans_edges = the_first_iteration ? number_of_edges[beta] : number_of_edges[alpha];
	ans_edges = number_of_edges[D_u_rank];
	for (int i = 0; i < ans.size; i++) ans_edges += indeg[ans.nodes[i]];
	if (ans.size != 0) max_rank = max(max_rank, index);
	return;
}
void Graph::get_dense_subgraph_from_core() {
	// printf("----- First iteration -----\n");
	for (alpha = 0; ; alpha++) {
		max_rank = 0; memset(r, -1, N * sizeof(int));
		for (int i = 0; i < N; i++) sorted[i] = i;
		for (int i = 0; i <= beta_max; i++) computed[i] = false, number_of_edges[i] = 0, number_of_nodes[i] = 0;
		// Get D[0] and D[alpha_max]
		D.clear();
		for (int u = 0; u < U; u++)
			if (undeg[u] > alpha) {
				D.insert(u);
				number_of_edges[0] += undeg[u];
				for (int j = 0; j < undeg[u]; j++) {
					Edge& ne = e[adj[u][j]];
					if (!D.in[ne.v])
						D.insert(ne.v);
				}
			}
		for (int i = 0; i < D.size; i++) {
			int x = D.nodes[i];
			r[x] = 0;
			if (in_U(x)) continue;
			for (int j = 0; j < undeg[x]; j++) {
				Edge& ne = e[adj[x][j]];
				if (ne.to != ne.u)
					indeg[ne.v]--, indeg[ne.u]++, ne.to = ne.u;
			}
		}
		change_sorted(0, N, D), computed[0] = true, position[0] = N - D.size, number_of_nodes[0] = D.size;
		computed[beta_max] = true, position[beta_max] = N, number_of_edges[beta_max] = 0, number_of_nodes[beta_max] = 0;
		Divide(0, beta_max, true);
		// printf("- %-20s: %d, %d\n", "max non-empty D_{alpha, beta}", alpha, max_rank);
		// analyze_r();
		if (alpha >= max_rank)
			break;
	}
	
	// printf("----- Second iteration -----\n");
	for (beta = 0; ; beta++) {
		max_rank = 0; memset(r, -1, N * sizeof(int));
		for (int i = 0; i < N; i++) sorted[i] = i;
		for (int i = 0; i <= alpha_max; i++) computed[i] = false, number_of_edges[i] = 0, number_of_nodes[i] = 0;
		// Get D[0] and D[alpha_max]
		D.clear();
		for (int v = U; v < N; v++)
			if (undeg[v] > beta) {
				D.insert(v);
				number_of_edges[0] += undeg[v];
				for (int j = 0; j < undeg[v]; j++) {
					Edge& ne = e[adj[v][j]];
					if (!D.in[ne.u])
						D.insert(ne.u);
				}
			}
		for (int i = 0; i < D.size; i++) {
			int x = D.nodes[i];
			r[x] = 0;
			if (!in_U(x)) continue;
			for (int j = 0; j < undeg[x]; j++) {
				Edge& ne = e[adj[x][j]];
				if (ne.to != ne.v)
					indeg[ne.u]--, indeg[ne.v]++, ne.to = ne.v;
			}
		}
		change_sorted(0, N, D), computed[0] = true, position[0] = N - D.size, number_of_nodes[0] = D.size;
		computed[alpha_max] = true, position[alpha_max] = N, number_of_edges[alpha_max] = 0, number_of_nodes[alpha_max] = 0;
		Divide(0, alpha_max, false);
		// printf("- %-20s: %d, %d\n", "max non-empty D_{alpha, beta}", max_rank, beta);
		// analyze_r();
		if (max_rank <= beta)
			break;
	}


}
bool Graph::DinicBFS(int D_l_rank, int D_u_rank) {
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
			check(r[from] >= D_l_rank, "r[from] < D_l_rank 2");
			if (r[from] >= D_u_rank) continue;
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
bool Graph::DinicDFS(int x, int D_l_rank, int D_u_rank) {
	if (in_S(x)) {
		indeg[x]++, indeg[e[parent[x]].to]--, e[parent[x]].to = x;
		return true;
	}
	for (int& j = cur[x]; j < undeg[x]; j++) {
		Edge& ne = e[adj[x][j]];
		if (ne.to != x) continue;
		int from = in_U(ne.to) ? ne.v : ne.u;
		check(r[from] >= D_l_rank, "r[from] < D_l_rank 3");
		if (r[from] >= D_u_rank) continue;
		if ((dist[from] != dist[x] + 1) && !in_S(from)) continue;
		parent[from] = adj[x][j];
		if (DinicDFS(from, D_l_rank, D_u_rank)) {
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
bool analyze_r_cmp(int x, int y) {
	return G.r[x] > G.r[y];
}
void Graph::analyze_r() {
	int* sorted_nodes = (int*)malloc(N * sizeof(int));
	for (int i = 0; i < N; i++)	sorted_nodes[i] = i;
	sort(sorted_nodes, sorted_nodes + N, analyze_r_cmp);
	for (int i = 1; i < N; i++) {
		if (r[sorted_nodes[i]] < r[sorted_nodes[i - 1]]) {
			// printf("- %-20s: %d, %d\n", "sorted_nodes[i] size", r[sorted_nodes[i - 1]], i);
		}
	}
	return;
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
	/*
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
	}*/
}
void Graph::display() {
	for (int i = 0; i < M; i++) {
		Edge& ne = e[i];
		int from = in_U(ne.to) ? ne.v : ne.u;
		printf("%d %d\n", from, ne.to);
	}
}
void Graph::check_correctness() {
	return;
}

int main(int argc, char** argv) {
	// argc = 2;
	// strcpy(argv[1], "../BiGraphs/marvel.txt");

	check(argc == 2, "The number of arguments of main are incorrect. ");
	char dataset_address[1000]; strcpy(dataset_address, argv[1]);

	double runtime;
	printf("----------Now processing %s----------\n", dataset_address);
	printf("- %-20s: %s\n", "Algorithm used", "DD+");

	timer.start();
	G.read_graph_from_dataset(dataset_address);
	timer.end(); runtime = timer.time();
	printf("- %-20s: %d, %d, %d\n", "|E|, |U|, |V|", G.M, G.U, G.V);
	printf("- %-20s: %lf\n", "Read graph time", runtime);

	timer.start();
	G.initialize_orientation();
	timer.end(); runtime = timer.time();
	printf("- %-20s: %lf\n", "Initialization time", runtime);

	timer.start();
	G.get_dense_subgraph_from_core();
	timer.end(); runtime = timer.time();
	printf("- %-20s: %lf\n", "Get all layers time", runtime);

	G.output_dense_subgraph();
	G.check_correctness();

	return 0;
}