#include "Entry.h"
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <set>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <chrono>
#include <iostream>
#include <tuple> // Add at the top of the file
#include <utility> // For std::pair
#include <limits>
#include <random>
#include <numeric>  // for accumulate
#include <vector>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <utility>
#include <algorithm>
#include <iostream>

using namespace std;

// Use pii (pair<int, int>) for grid cell coordinates.
using pii = pair<int, int>;

struct pair_hash {
    size_t operator()(const pii& p) const {
        return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
    }
};


enum class PivotMethod { ED, HE };
enum class HeuristicCombination { SUM, MAX };

struct HeuristicConfig {
    int fmDims;
    int dhDims;
    PivotMethod pivotMethod;
    HeuristicCombination combMethod;
};

// ---------------- Helper: Bounds Check ----------------
bool is_in_bounds(int r, int c, int rows, int cols) {
    return (r >= 0 && r < rows && c >= 0 && c < cols);
}

// ---------------- Helper: Get 8-connected Neighbors ----------------
vector<pii> get_neighbors_8(const pii& current, const pair<int,int>& grid_size,
                              const unordered_set<pii, pair_hash>& obstacles) {
    vector<pii> neighbors;
    int r = current.first, c = current.second;
    int rows = grid_size.first, cols = grid_size.second;
    const int DX[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    const int DY[8] = {-1,  0,  1,-1, 1,-1, 0, 1};
    for (int i = 0; i < 8; ++i) {
        int nr = r + DX[i], nc = c + DY[i];
        if (!is_in_bounds(nr, nc, rows, cols)) continue;
        if (obstacles.count({nr, nc})) continue;
        neighbors.emplace_back(nr, nc);
    }
    // Remove diagonal neighbors that would cut corners.
    const int DX_d[4] = {-1, -1, 1, 1};
    const int DY_d[4] = {-1,  1,-1, 1};
    for (int i = 0; i < 4; ++i) {
        int nr = r + DX_d[i], nc = c + DY_d[i];
        if (!is_in_bounds(nr, nc, rows, cols)) continue;
        if (obstacles.count({r, nc}) || obstacles.count({nr, c})) {
            auto it = find(neighbors.begin(), neighbors.end(), make_pair(nr, nc));
            if (it != neighbors.end())
                neighbors.erase(it);
        }
    }
    return neighbors;
}

// ---------------- Movement Cost ----------------
float move_cost(const pii& a, const pii& b) {
    return (a.first != b.first && a.second != b.second) ? sqrt(2.0) : 1.0;
}

// ---------------- Path Reconstruction for A* ----------------
vector<pii> reconstruct_path(unordered_map<pii, pii, pair_hash>& came_from, pii current) {
    vector<pii> path = { current };
    while (came_from.count(current)) {
        current = came_from[current];
        path.push_back(current);
    }
    reverse(path.begin(), path.end());
    return path;
}

// -------------- FM PREPROCESSING FUNCTIONS --------------

// Structure representing an edge in the FM graph.
struct Edge {
    int u, v;
    double w;
};

// Type alias for an adjacency list.
using AdjList = vector<vector<pair<int, double>>>;

// Build an adjacency list from a given edge list.
AdjList buildAdjList(const vector<Edge>& edges, int n) {
    AdjList adj(n);
    for (const auto &edge : edges) {
        adj[edge.u].push_back({edge.v, edge.w});
        adj[edge.v].push_back({edge.u, edge.w});
    }
    return adj;
}

// Dijkstra's algorithm used by FM to compute shortest-path distances.
vector<double> dijkstra(const AdjList &adj, int src, int n) {
    vector<double> dist(n, numeric_limits<double>::max());
    priority_queue<pair<double,int>, vector<pair<double,int>>, greater<pair<double,int>>> pq;
    dist[src] = 0.0;
    pq.push({0.0, src});
    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (d > dist[u])
            continue;
        for (auto &neighbor : adj[u]) {
            int v = neighbor.first;
            double w = neighbor.second;
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                pq.push({dist[v], v});
            }
        }
    }
    return dist;
}

// Identify a farthest pair of nodes in the graph using τ iterations.
pair<int, int> getFarthestPair(const AdjList &adj, int n, int tau = 10) {
    int a = rand() % n;
    int b = a;
    for (int t = 0; t < tau; t++) {
        vector<double> distA = dijkstra(adj, a, n);
        double maxDist = -1.0;
        for (int i = 0; i < n; i++) {
            if (distA[i] < numeric_limits<double>::max() && distA[i] > maxDist) {
                maxDist = distA[i];
                b = i;
            }
        }
        vector<double> distB = dijkstra(adj, b, n);
        maxDist = -1.0;
        for (int i = 0; i < n; i++) {
            if (distB[i] < numeric_limits<double>::max() && distB[i] > maxDist) {
                maxDist = distB[i];
                a = i;
            }
        }
    }
    return {a, b};
}

// Get the farthest pair and their distances.
std::tuple<int, int, vector<double>, vector<double>> getFarthestPairAndDistances(const AdjList &adj, int n, int tau = 10) {
    int a = rand() % n;
    int b = a;
    vector<double> distA, distB;

    for (int t = 0; t < tau; t++) {
        distA = dijkstra(adj, a, n);
        double maxDist = -1.0;
        for (int i = 0; i < n; i++) {
            if (distA[i] < std::numeric_limits<double>::max() && distA[i] > maxDist) {
                maxDist = distA[i];
                b = i;
            }
        }

        distB = dijkstra(adj, b, n);
        maxDist = -1.0;
        for (int i = 0; i < n; i++) {
            if (distB[i] < std::numeric_limits<double>::max() && distB[i] > maxDist) {
                maxDist = distB[i];
                a = i;
            }
        }
    }

    return std::make_tuple(a, b, distA, distB);
}


// FastMap algorithm: returns an embedding (a vector of coordinate vectors for each node).
// This version uses the paper's indexing: Start with K = 1 and decrement the available budget (Kmax) each iteration.
vector<vector<double>> fastMap(vector<Edge> edges, int n, int Kmax, double epsilon) {
    // Initialize: we work on a copy of the edges; w0 represents the current (residual) edge weights.
    vector<Edge> w0 = edges;
    
    // Build the initial residual graph from w0.
    AdjList workAdj = buildAdjList(w0, n);
    
    // Initialize the embedding for each node; each entry will collect one coordinate per iteration.
    vector<vector<double>> embedding(n);
    
    // Set the starting dimension index K. (K = 1 in the paper)
    int K = 1;
    
    // The loop will iterate until we run out of dimension budget (Kmax > 0)
    while (Kmax > 0) {
        // Build the current residual graph G0 = (V, E, w0) using the current edge weights.
        workAdj = buildAdjList(w0, n);
        
        // Line 4: Get farthest pair in G0 (na, nb)
        pair<int, int> farPair = getFarthestPair(workAdj, n);
        int na = farPair.first, nb = farPair.second;
        
        // Line 5: Compute shortest path trees from na and nb.
        vector<double> distA = dijkstra(workAdj, na, n);
        vector<double> distB = dijkstra(workAdj, nb, n);
        
        // Get the distance between na and nb.
        double dab = distA[nb];
        // Line 6-7: If the farthest distance is below threshold epsilon, stop.
        if (dab < epsilon)
            break;
        
        // Line 8-9: For each node v, compute the new coordinate for this dimension.
        //  [p_v]_K = (d(na, v) + d(na, nb) - d(nb, v)) / 2
        const double INF = std::numeric_limits<double>::max();  // from <limits>
        double coord;

        for (int v = 0; v < n; v++) {
            // double coord = (distA[v] + dab - distB[v]) / 2.0;

            if (distA[v] == INF || distB[v] == INF || dab == INF) {
                coord = 0.0; // fallback if any component is invalid
            } else {
                coord = (distA[v] + dab - distB[v]) / 2.0;
            }
            embedding[v].push_back(coord);
        }
        
        // Line 10-11: For each edge (u, v), update its weight with the residual component.
        for (auto &edge : w0) {
            double diff = fabs(embedding[edge.u].back() - embedding[edge.v].back());
            edge.w = max(0.0, edge.w - diff);
        }
        
        // Line 12: Update indexing.
        // Increment the current dimension K and reduce the remaining dimension budget Kmax.
        K = K + 1;
        Kmax = Kmax - 1;
    }
    
    // The final embedding contains K-1 coordinates for each node.
    return embedding;
}




enum class FastMapNorm { L1, L2 };


// Enhanced FastMap
// - edges:      original graph edges
// - n:          number of nodes
// - Kmax:       number of dimensions to build
// - epsilon:    stop when pivot‐pair distance < epsilon
// - normMode:   use L1 or L2 for both heuristic and residual update
// - pivotsOut:  if non‐null, filled with the chosen pivot pairs (na, nb)
vector<vector<double>> fastMapEnhanced(
    const vector<Edge>& edges,
    int n,
    int Kmax,
    double epsilon,
    FastMapNorm normMode = FastMapNorm::L1,
    vector<pair<int,int>>* pivotsOut = nullptr
) {
    // 1) Guards
    if (n <= 0) return {};         // no nodes
    if (Kmax <= 0) return vector<vector<double>>(n);
    const double INF = numeric_limits<double>::max();

    // 2) Prepare
    vector<Edge> w0 = edges;             // residual weights
    vector<vector<double>> embedding(n);  // embedding[v][dim]
    vector<pair<int,int>> pivots;         // record (na,nb)

    // 3) Main loop
    while (Kmax-- > 0) {
        // Build residual graph
        AdjList adj = buildAdjList(w0, n);

        // Select pivot pair
        auto [na, nb] = getFarthestPair(adj, n);
        pivots.emplace_back(na, nb);

        // Compute distances in residual graph
        vector<double> distA = dijkstra(adj, na, n);
        vector<double> distB = dijkstra(adj, nb, n);

        double dab = distA[nb];
        if (dab == INF || dab < epsilon) break;

        // 4) Embed each node
        for (int v = 0; v < n; v++) {
            double coord;
            if (distA[v] == INF || distB[v] == INF) {
                coord = 0.0;  // fallback if unreachable
            } else {
                coord = (distA[v] + dab - distB[v]) / 2.0;
            }
            embedding[v].push_back(coord);
        }

        // 5) Residual normalization: subtract the SAME norm from each edge
        for (auto &e : w0) {
            double x_u = embedding[e.u].back();
            double x_v = embedding[e.v].back();
            double diff = (normMode == FastMapNorm::L2)
                          ? sqrt((x_u - x_v)*(x_u - x_v))
                          : fabs(x_u - x_v);
            e.w = max(0.0, e.w - diff);
        }
    }

    // 6) Output pivot pairs if requested
    if (pivotsOut) *pivotsOut = pivots;

    return embedding;
}


// FM-based heuristic function: returns L1 distance between two embedding vectors.
float fm_heuristic(const vector<double>& embA, const vector<double>& embB) {
    float distance = 0.0f;
    for (size_t i = 0; i < embA.size(); i++) {
        distance += fabs(embA[i] - embB[i]); // L1 distance
    }
    return distance;
}



unordered_map<pii, vector<double>, pair_hash> fmEmbedding_HE;
unordered_map<pii, int, pair_hash> cellToIndex_E;
vector<pii> indexToCell_E;

// Assume an Edge is defined as follows:
struct Edge_E {
    int u, v;
    double w;
};

// Define the adjacency list type (you likely already have this)
typedef vector<vector<pair<int, double>>> AdjList;


// Octile distance function (using grid coordinates)
double octile(const pii &a, const pii &b) {
    int dr = abs(a.first - b.first);
    int dc = abs(a.second - b.second);
    return max(dr, dc) + (sqrt(2.0) - 1.0) * min(dr, dc);
}

// Build an adjacency list from a given edge list.
AdjList buildAdjList_E(const vector<Edge_E>& edges, int n) {
    AdjList adj(n);
    for (const auto &edge : edges) {
        adj[edge.u].push_back({edge.v, edge.w});
        adj[edge.v].push_back({edge.u, edge.w});
    }
    return adj;
}
//------------------------------------------------------------
// HE Pivot Selection Routine
pair<int,int> getFarthestPair_HE(const AdjList &workAdj, int n) {
    if (n == 0) {
        std::cerr << "[ERROR] getFarthestPair_HE called with n = 0." << std::endl;
        exit(1);
    }

    int t = rand() % n;
    if (t < 0 || t >= indexToCell_E.size()) {
        std::cerr << "[ERROR] indexToCell_E out-of-bounds access for t = " << t << std::endl;
        exit(1);
    }

    vector<double> distT = dijkstra(workAdj, t, n);

    double maxVal = -numeric_limits<double>::infinity();
    int p1 = t;

    for (int v = 0; v < n; v++) {
        if (v < 0 || v >= indexToCell_E.size()) {
            std::cerr << "[ERROR] indexToCell_E out-of-bounds access for v = " << v << std::endl;
            continue;
        }

        pii cell_t = indexToCell_E[t];
        pii cell_v = indexToCell_E[v];

        double d_val = distT[v];
        if (d_val == numeric_limits<double>::max()) continue; // unreachable

        double h_val = octile(cell_t, cell_v);
        double value = 3.0 * d_val - 2.0 * h_val;
        if (value > maxVal) {
            maxVal = value;
            p1 = v;
        }
    }

    vector<double> distP1 = dijkstra(workAdj, p1, n);
    maxVal = -numeric_limits<double>::infinity();
    int p2 = p1;

    for (int v = 0; v < n; v++) {
        if (v < 0 || v >= indexToCell_E.size()) continue;

        pii cell_p1 = indexToCell_E[p1];
        pii cell_v = indexToCell_E[v];

        double d_val = distP1[v];
        if (d_val == numeric_limits<double>::max()) continue;

        double h_val = octile(cell_p1, cell_v);
        double value = 3.0 * d_val - 2.0 * h_val;
        if (value > maxVal) {
            maxVal = value;
            p2 = v;
        }
    }

    return {p1, p2};
}


//------------------------------------------------------------
// FastMap_HE: FastMap embedding with HE pivot selection.
vector<vector<double>> fastMap_HE(vector<Edge_E> edges, int n, int Kmax, double epsilon) {
    // Work on a copy of the edges; w0 holds the current (residual) edge weights.
    vector<Edge_E> w0 = edges;
    
    // Build the initial residual graph.
    AdjList workAdj = buildAdjList_E(w0, n);
    
    // Initialize the embedding for each node.
    vector<vector<double>> embedding(n);
    
    // Loop until dimension budget runs out.
    while (Kmax > 0) {
        workAdj = buildAdjList_E(w0, n);
        
        // Use HE pivot selection to get the farthest pair.
        pair<int,int> farPair = getFarthestPair_HE(workAdj, n);
        int na = farPair.first, nb = farPair.second;
        
        // Compute shortest path distances from the pivots.
        vector<double> distA = dijkstra(workAdj, na, n);
        vector<double> distB = dijkstra(workAdj, nb, n);
        
        double dab = distA[nb]; // distance between pivots
        if (dab < epsilon)
            break; // Terminate if the farthest distance is below threshold.
        
        // For each node v, compute the coordinate for this dimension.
        // [coord]_K = (d(na,v) + d(na,nb) - d(nb,v)) / 2
        const double INF = std::numeric_limits<double>::max();  // from <limits>
        double coord;
        for (int v = 0; v < n; v++) {
            // double coord = (distA[v] + dab - distB[v]) / 2.0;

            if (distA[v] == INF || distB[v] == INF || dab == INF) {
                coord = 0.0; // fallback if any component is invalid
            } else {
                coord = (distA[v] + dab - distB[v]) / 2.0;
            }
            embedding[v].push_back(coord);
        }
        
        // Update residual edge weights: subtract the captured distance.
        for (auto &edge : w0) {
            double diff = fabs(embedding[edge.u].back() - embedding[edge.v].back());
            edge.w = max(0.0, edge.w - diff);
        }
        
        // One embedding dimension is built; decrement the budget.
        Kmax--;
    }
    
    return embedding;
}


// ---------------- Global Storage ----------------
unordered_set<pii, pair_hash> global_obstacles;
int global_width = 0, global_height = 0;



// Build a graph from the grid: free (non-obstacle) cells become nodes; edges connect 8-connected free cells.
void buildGraphFromGrid_E(const pair<int, int>& grid_size, const unordered_set<pii, pair_hash>& obstacles,
                        vector<Edge_E>& edges) {
    cellToIndex_E.clear();
    indexToCell_E.clear();
    int n = 0;
    int rows = grid_size.first, cols = grid_size.second;
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            pii cell = {r, c};
            if (obstacles.count(cell) == 0) { // free cell
                cellToIndex_E[cell] = n++;
                indexToCell_E.push_back(cell);
            }
        }
    }
    // For each free cell, add an edge to each free 8-connected neighbor (only once).
    for (auto &kv : cellToIndex_E) {
        pii cell = kv.first;
        int idx = kv.second;
        vector<pii> nbrs = get_neighbors_8(cell, grid_size, obstacles);
        for (auto &nb : nbrs) {
            if (cellToIndex_E.find(nb) != cellToIndex_E.end()) {
                int nb_idx = cellToIndex_E[nb];
                if (idx < nb_idx) { // avoid duplicates
                    float cost = move_cost(cell, nb);
                    edges.push_back({idx, nb_idx, cost});
                }
            }
        }
    }
}



void build_fmEmbedding_HE(int width, int height, int Kmax_HE, double epsilon) {
    // int Kmax_HE = 9;           // Number of FM dimensions
    // double epsilon = 1e-2;     // Termination threshold
    vector<Edge_E> edges;

    // Ensure cellToIndex_E and indexToCell_E are aligned
    if (cellToIndex_E.empty() || indexToCell_E.empty()) {
        std::cerr << "[FATAL] cellToIndex_E or indexToCell_E is empty — cannot build FM embedding." << std::endl;
        exit(1);
    }

    int numFreeCells = cellToIndex_E.size();

    if (indexToCell_E.size() != numFreeCells) {
        std::cerr << "[FATAL] indexToCell_E and cellToIndex_E sizes do not match." << std::endl;
        exit(1);
    }

    if (width <= 0 || height <= 0) {
        std::cerr << "[FATAL] Invalid grid dimensions: " << width << "x" << height << std::endl;
        exit(1);
    }

    // Build edge list from grid
    pair<int, int> grid_size = {height, width};
    buildGraphFromGrid_E(grid_size, global_obstacles, edges);

    if (edges.empty()) {
        std::cerr << "[FATAL] Edge list is empty. Graph construction failed." << std::endl;
        exit(1);
    }

    // Run FastMap with HE pivot selection
    vector<vector<double>> embedding_HE = fastMap_HE(edges, numFreeCells, Kmax_HE, epsilon);

    if (embedding_HE.size() != numFreeCells) {
        std::cerr << "[FATAL] FastMap embedding size mismatch." << std::endl;
        exit(1);
    }

    // Store embedding into global map
    fmEmbedding_HE.clear();
    for (int i = 0; i < numFreeCells; i++) {
        if (i < 0 || i >= indexToCell_E.size()) {
            std::cerr << "[ERROR] indexToCell_E out of bounds at i=" << i << std::endl;
            continue;
        }
        fmEmbedding_HE[indexToCell_E[i]] = embedding_HE[i];
    }

    // std::cerr << "[INFO] FM_HE embedding built for " << numFreeCells << " free cells with "
    //           << Kmax_HE << " dimensions." << std::endl;
}



// Global variables for FM.
unordered_map<pii, vector<double>, pair_hash> fmEmbedding;
unordered_map<pii, int, pair_hash> cellToIndex;
vector<pii> indexToCell;
vector<double> dhVector;  // or vector<float> dhVector;
vector<vector<double>> dhVectors;  // Each inner vector is one landmark's Dijkstra result

struct FMContext {
    long long heuristic_micro_time;  // or double heuristic_seconds;
    // You can add more things here later (e.g., Kmax, embedding stats, etc.)
};


vector<vector<double>> dhVectors_E;  // Each inner vector is one landmark's Dijkstra result
vector<int> expansion_log;


// Build a graph from the grid: free (non-obstacle) cells become nodes; edges connect 8-connected free cells.
void buildGraphFromGrid(const pair<int, int>& grid_size, const unordered_set<pii, pair_hash>& obstacles,
                        vector<Edge>& edges) {
    cellToIndex.clear();
    indexToCell.clear();
    int n = 0;
    int rows = grid_size.first, cols = grid_size.second;
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            pii cell = {r, c};
            if (obstacles.count(cell) == 0) { // free cell
                cellToIndex[cell] = n++;
                indexToCell.push_back(cell);
            }
        }
    }
    // For each free cell, add an edge to each free 8-connected neighbor (only once).
    for (auto &kv : cellToIndex) {
        pii cell = kv.first;
        int idx = kv.second;
        vector<pii> nbrs = get_neighbors_8(cell, grid_size, obstacles);
        for (auto &nb : nbrs) {
            if (cellToIndex.find(nb) != cellToIndex.end()) {
                int nb_idx = cellToIndex[nb];
                if (idx < nb_idx) { // avoid duplicates
                    float cost = move_cost(cell, nb);
                    edges.push_back({idx, nb_idx, cost});
                }
            }
        }
    }
}


// Build a graph from the grid: free (non-obstacle) cells become nodes; edges connect 8-connected free cells.

enum HeuristicMode {
    DH,
    FM,
    FM_DH_ED,
    MAX_DH_FM_ED,
    MAX_DH_FMDH_ED,
    FM_DH_HE,
    MAX_DH_FM_HE,
    MAX_DH_FMDH_HE
};


// -------------- A* Algorithm using FM-Based Heuristic --------------
vector<pii> a_star(const pii& start, const pii& goal, const pair<int, int>& grid_size,
                   const unordered_set<pii, pair_hash>& obstacles) {
    int expansions = 0;

    using QueueElement = pair<float, pii>;
    priority_queue<QueueElement, vector<QueueElement>, greater<>> open_list;
    open_list.emplace(0.0f, start);

    unordered_set<pii, pair_hash> closed_list;
    unordered_map<pii, float, pair_hash> g;
    unordered_map<pii, float, pair_hash> f;
    unordered_map<pii, pii, pair_hash> came_from;
    
    HeuristicMode selected_mode = MAX_DH_FM_ED;  // change this to try different options
    std::function<float(const pii&, const pii&)> get_h;
    switch (selected_mode) {

        case FM: {
        // FM heuristic
            get_h = [&](const pii& a, const pii& b) -> float {
                // int dr = abs(a.first - b.first);
                // int dc = abs(a.second - b.second);
                // float octile = max(dr, dc) + (sqrt(2.0f) - 1.0f) * min(dr, dc);
                if (fmEmbedding.count(a) && fmEmbedding.count(b) &&
                    !fmEmbedding[a].empty() && !fmEmbedding[b].empty()) {
                        float fm_h = fm_heuristic(fmEmbedding[a], fmEmbedding[b]);
                        // if (fm_h > octile) {
                        //     std::cerr << "[FM] FM heuristic is greater than octile distance." << "FM = " << fm_h << " || octile = " << octile << std::endl;
                        // }
                    return fm_h;
                } else {
                    std::cerr << "[FM] FM heuristic not available, using octile." << std::endl;
                    int dr = abs(a.first - b.first);
                    int dc = abs(a.second - b.second);
                    float octile = max(dr, dc) + (sqrt(2.0f) - 1.0f) * min(dr, dc);
                    return octile;
                }
            }; 
        break;
        }

        case DH: {
            // DH heuristic
            get_h = [&](const pii &a, const pii &b) -> float {
                int idx_a = a.first * global_width + a.second;
                int idx_b = b.first * global_width + b.second;
                float max_h = 0.0f;
                for (const auto& dh : dhVectors) {
                    double dist_a = dh[idx_a];
                    double dist_b = dh[idx_b];
                    // Check for unreachable cells; if either is unreachable, skip this pivot.
                    if(dist_a == numeric_limits<double>::max() || dist_b == numeric_limits<double>::max())
                        continue;
                    float h = fabs(dist_a - dist_b);
                    if (h > max_h)
                        max_h = h;
                }
                // Compute Octile distance (admissible for 8-connected grid)
                int dr = abs(a.first - b.first);
                int dc = abs(a.second - b.second);
                float octile = max(dr, dc) + (sqrt(2.0f) - 1) * min(dr, dc);

                // Return the maximum of the two heuristics
                if (octile > max_h)         
                    max_h = octile;
                return max_h;
            };
            break;
        }
    
        case FM_DH_ED:{
            // FM + DH with ED
            get_h = [&](const pii &a, const pii &b) -> float {
                float fm_val = 0.0f;
                if (fmEmbedding.count(a) && fmEmbedding.count(b) &&
                    !fmEmbedding[a].empty() && !fmEmbedding[b].empty()) {
                    fm_val = fm_heuristic(fmEmbedding[a], fmEmbedding[b]);
                } else {
                    // Fallback: use octile distance if the FM embedding is missing.
                    int dr = abs(a.first - b.first);
                    int dc = abs(a.second - b.second);
                    fm_val = max(dr, dc) + (sqrt(2.0) - 1) * min(dr, dc);
                }

                int idx_a = a.first * global_width + a.second;
                int idx_b = b.first * global_width + b.second;
            
                float max_h = 0.0f;
                for (const auto& dh : dhVectors) {
                    double dist_a = dh[idx_a];
                    double dist_b = dh[idx_b];
                    // Check for unreachable cells; if either is unreachable, skip this pivot.
                    if(dist_a == numeric_limits<double>::max() || dist_b == numeric_limits<double>::max())
                        continue;
                    float h = fabs(dist_a - dist_b);
                    if (h > max_h)
                        max_h = h;
                }
                float dh_val = max_h;
                return fm_val + dh_val;
            };
            break;
        }

        case MAX_DH_FM_ED: {
            // max[DH5, FM5] using ED pivots
            get_h = [&](const pii &a, const pii &b) -> float {
                float fm_val = 0.0f;
                if (fmEmbedding.count(a) && fmEmbedding.count(b) &&
                    !fmEmbedding[a].empty() && !fmEmbedding[b].empty()) {
                    fm_val = fm_heuristic(fmEmbedding[a], fmEmbedding[b]);
                } else {
                    // Fallback: use octile distance
                    int dr = abs(a.first - b.first);
                    int dc = abs(a.second - b.second);
                    fm_val = max(dr, dc) + (sqrt(2.0f) - 1.0f) * min(dr, dc);
                }
        
                int idx_a = a.first * global_width + a.second;
                int idx_b = b.first * global_width + b.second;
        
                float dh_val = 0.0f;
                for (const auto& dh : dhVectors) {
                    double dist_a = dh[idx_a];
                    double dist_b = dh[idx_b];
                    if (dist_a == numeric_limits<double>::max() || dist_b == numeric_limits<double>::max())
                        continue;
                    float h = fabs(dist_a - dist_b);
                    if (h > dh_val)
                        dh_val = h;
                }
        
                return std::max(fm_val, dh_val);  // max[DH5, FM5]
            };
            break;
        }
        
        case MAX_DH_FMDH_ED: {
            get_h = [&](const pii &a, const pii &b) -> float {
                // --- FM4 ---
                float fm_val = 0.0f;
                if (fmEmbedding.count(a) && fmEmbedding.count(b) &&
                    !fmEmbedding[a].empty() && !fmEmbedding[b].empty()) {
                    fm_val = fm_heuristic(fmEmbedding[a], fmEmbedding[b]);
                } else {
                    // Fallback to octile distance
                    int dr = abs(a.first - b.first);
                    int dc = abs(a.second - b.second);
                    fm_val = max(dr, dc) + (sqrt(2.0f) - 1.0f) * min(dr, dc);
                }
        
                // --- DH[0] (first vector) ---
                int idx_a = a.first * global_width + a.second;
                int idx_b = b.first * global_width + b.second;
        
                const auto& dhEmbedding0 = dhVectors[0];
                float dh1_val = fabs(dhEmbedding0[idx_a] - dhEmbedding0[idx_b]);
        
                float hybrid_val = fm_val + dh1_val;
        
                // --- DH5 (max over all landmarks) ---
                float dh5_val = 0.0f;
                for (const auto& dh : dhVectors) {
                    double dist_a = dh[idx_a];
                    double dist_b = dh[idx_b];
                    if (dist_a == numeric_limits<double>::max() || dist_b == numeric_limits<double>::max())
                        continue;
                    float h = fabs(dist_a - dist_b);
                    if (h > dh5_val)
                        dh5_val = h;
                }
        
                // Return max[FM4 + DH1, DH5]
                return std::max(hybrid_val, dh5_val);
            };
            break;
        }
        
        case FM_DH_HE:{
            // FM + DH with HE
            get_h = [&](const pii &a, const pii &b) -> float {
                // Hybrid: FMk+DH computed using HE-based FM.
                // FM
                float fm4_val = 0.0f;
                // Octile distance
                float octile = 0.0f;
                int dr = abs(a.first - b.first);
                int dc = abs(a.second - b.second);
                octile = max(dr, dc) + (sqrt(2.0f) - 1) * min(dr, dc);

                if (fmEmbedding_HE.count(a) && fmEmbedding_HE.count(b) &&
                    !fmEmbedding_HE[a].empty() && !fmEmbedding_HE[b].empty()) {
                    fm4_val = fm_heuristic(fmEmbedding_HE[a], fmEmbedding_HE[b]);
                } else {
                    fm4_val = octile;
                }
                
                int idx_a = a.first * global_width + a.second;
                int idx_b = b.first * global_width + b.second;

                // // DH1
                // const auto& dhEmbedding_HE = dhVectors[0];  // Assuming you used computeDHVectors(1)
                // float dh1_val = fabs(dhEmbedding_HE[idx_a] - dhEmbedding_HE[idx_b]);  // DH for hybrid computed with HE-based FM.

                float max_h = 0.0f;
                for (const auto& dh_E : dhVectors_E) {
                    double dist_a = dh_E[idx_a];
                    double dist_b = dh_E[idx_b];
                    // Check for unreachable cells; if either is unreachable, skip this pivot.
                    if(dist_a == numeric_limits<double>::max() || dist_b == numeric_limits<double>::max())
                        continue;
                    float h = fabs(dist_a - dist_b);
                    if (h > max_h)
                        max_h = h;
                }
                float dh__E_val = max_h;

                return fm4_val + dh__E_val;
            };
            break;
        }
            
        case MAX_DH_FMDH_HE:{
            get_h = [&](const pii &a, const pii &b) -> float {
                // Hybrid: FM4+DH computed using HE-based FM.
                float fm4_val = 0.0f;
                float octile = 0.0f;
                int dr = abs(a.first - b.first);
                int dc = abs(a.second - b.second);
                octile = max(dr, dc) + (sqrt(2.0f) - 1) * min(dr, dc);

                if (fmEmbedding_HE.count(a) && fmEmbedding_HE.count(b) &&
                    !fmEmbedding_HE[a].empty() && !fmEmbedding_HE[b].empty()) {
                    fm4_val = fm_heuristic(fmEmbedding_HE[a], fmEmbedding_HE[b]);
                } else {
                    fm4_val = octile;
                }
                
                int idx_a = a.first * global_width + a.second;
                int idx_b = b.first * global_width + b.second;

                const auto& dhEmbedding_HE = dhVectors[0];  // Assuming you used computeDHVectors(1)
                float dh1_val = fabs(dhEmbedding_HE[idx_a] - dhEmbedding_HE[idx_b]);  // DH for hybrid computed with HE-based FM.
            
                float hybrid_val = fm4_val + dh1_val;
                // return hybrid_val;
            
                // Standard DHk using ED or HE for DH might be the same.
                float dh5_val = 0.0f;
                for (const auto& dh : dhVectors) {
                    double dist_a = dh[idx_a];
                    double dist_b = dh[idx_b];
                    if (dist_a == numeric_limits<double>::max() || dist_b == numeric_limits<double>::max())
                        continue;
                    float h = fabs(dist_a - dist_b);
                    if (h > dh5_val)
                        dh5_val = h;
                }

                // max[FM4+DH,DH5]
                return std::max(hybrid_val, dh5_val);
            };
            break;
        }
    
    }

    g[start] = 0.0;
    f[start] = g[start] + get_h(start, goal);

    while (!open_list.empty()) {
        pii current = open_list.top().second;
        open_list.pop();

        if (current == goal) {
            expansion_log.push_back(expansions); // ✅ log before returning
            return reconstruct_path(came_from, current);
        }

        // ✅ Count expanded node
        expansions++;
        closed_list.insert(current);

        for (const pii& nb : get_neighbors_8(current, grid_size, obstacles)) {
            if (closed_list.count(nb))
                continue;
            float tentative_g = g[current] + move_cost(current, nb);
            if (!g.count(nb) || tentative_g < g[nb]) {
                came_from[nb] = current;
                g[nb] = tentative_g;
                f[nb] = tentative_g + get_h(nb, goal);
                open_list.emplace(f[nb], nb);
            }
        }
    }
    expansion_log.push_back(expansions); // still log it for failure cases
    return {};  // No path found.
}


// Returns a vector<double> of size (width * height). The index for cell (r, c) is: index = r * width + c.
vector<double> dijkstraForLandmark(const pii& landmark,
                                   int width,
                                   int height,
                                   const unordered_set<pii, pair_hash>& obstacles) {
    int totalCells = width * height;
    if (width == 0 || height == 0) {
        std::cerr << "[ERROR] Grid width or height is zero." << std::endl;
        exit(1);
    }
    
    vector<double> distance(totalCells, numeric_limits<double>::max());
    pair<int, int> grid_size = {height, width};
    
    auto index = [width](const pii& cell) {
        return cell.first * width + cell.second;
    };
    
    // std::cerr << "Landmark: " << landmark.first << ", " << landmark.second << std::endl;
    // Initialize the landmark's own distance to 0.
    distance[index(landmark)] = 0.0;
    
    // Priority queue: pair of (distance, cell)
    using QueueElement = pair<double, pii>;
    priority_queue<QueueElement, vector<QueueElement>, greater<QueueElement>> pq;
    pq.push({0.0, landmark});
    
    while (!pq.empty()) {
        auto current = pq.top();
        pq.pop();
        double curDist = current.first;
        pii curCell = current.second;
        
        // If we already found a shorter path, skip.
        if (curDist > distance[index(curCell)])
            continue;
        
        vector<pii> neighbors = get_neighbors_8(curCell, grid_size, obstacles);
        for (const pii& neighbor : neighbors) {
            double stepCost = move_cost(curCell, neighbor);
            double newDist = curDist + stepCost;
            if (newDist < distance[index(neighbor)]) {
                distance[index(neighbor)] = newDist;
                pq.push({newDist, neighbor});
            }
        }
    }

    // double max_dist = 0.0;
    // for (double d : distance) {
    //     if (d < numeric_limits<double>::max()) {
    //         max_dist = max(max_dist, d);
    //     }
    // }
    // std::cerr << "Landmark max reachable distance: " << max_dist << std::endl;

    // int reachable = 0;
    // for (double d : distance) {
    //     if (d < numeric_limits<double>::max()) reachable++;
    // }
    // std::cerr << "Reachable cells from this landmark: " << reachable << std::endl;


    return distance;
}



void computeDHVectors(int numLandmarks) {
    dhVectors.clear();

    // Step 1: Collect all walkable cells
    vector<pii> free_cells;
    for (const auto& entry : cellToIndex) {
        free_cells.push_back(entry.first);
    }

    // Step 2: Choose the first pivot randomly
    vector<pii> pivots;
    std::mt19937 rng(42);  // or time(0) if you want nondeterministic
    std::uniform_int_distribution<int> dist(0, free_cells.size() - 1);
    pivots.push_back(free_cells[dist(rng)]);

    // Step 3: Iteratively pick the farthest cell from all previous pivots
    while (pivots.size() < numLandmarks) {
        float best_min_dist = -1.0;
        pii best_candidate;

        for (const pii& candidate : free_cells) {
            float min_dist = std::numeric_limits<float>::max();

            for (const pii& p : pivots) {
                int dr = abs(candidate.first - p.first);
                int dc = abs(candidate.second - p.second);
                float d = std::max(dr, dc) + (sqrt(2.0f) - 1) * std::min(dr, dc);  // Octile
                min_dist = std::min(min_dist, d);
            }

            if (min_dist > best_min_dist) {
                best_min_dist = min_dist;
                best_candidate = candidate;
            }
        }

        pivots.push_back(best_candidate);
    }

    // Step 4: Compute DH vectors for each pivot
    for (const pii& landmark : pivots) {
        dhVectors.push_back(dijkstraForLandmark(landmark, global_width, global_height, global_obstacles));
    }
}

void computeDHVectors_E(int numLandmarks) {
    dhVectors_E.clear();

    // Step 1: Collect all walkable cells
    vector<pii> free_cells;
    for (const auto& entry : cellToIndex_E) {
        free_cells.push_back(entry.first);
    }

    if (free_cells.empty()) {
        std::cerr << "[ERROR] No free cells available for DH pivot selection." << std::endl;
        exit(1);  // Or return, or throw, depending on your setup
    }
    
    // Step 2: Choose the first pivot randomly
    vector<pii> pivots;
    std::mt19937 rng(42);  // or time(0) if you want nondeterministic
    std::uniform_int_distribution<int> dist(0, free_cells.size() - 1);
    pivots.push_back(free_cells[dist(rng)]);

    // Step 3: Iteratively pick the farthest cell from all previous pivots
    while (pivots.size() < numLandmarks) {
        float best_min_dist = -1.0;
        pii best_candidate;

        for (const pii& candidate : free_cells) {
            float min_dist = std::numeric_limits<float>::max();

            for (const pii& p : pivots) {
                int dr = abs(candidate.first - p.first);
                int dc = abs(candidate.second - p.second);
                float d = std::max(dr, dc) + (sqrt(2.0f) - 1) * std::min(dr, dc);  // Octile
                min_dist = std::min(min_dist, d);
            }

            if (min_dist > best_min_dist) {
                best_min_dist = min_dist;
                best_candidate = candidate;
            }
        }

        pivots.push_back(best_candidate);
    }

    // Step 4: Compute DH vectors for each pivot
    for (const pii& landmark : pivots) {
        dhVectors_E.push_back(dijkstraForLandmark(landmark, global_width, global_height, global_obstacles));
    }
}

void computeCapturedHeuristic(
    const unordered_map<pii, int, pair_hash>& cellToIndex,
    const unordered_map<pii, vector<double>, pair_hash>& fmEmbedding,
    const vector<vector<double>>& dhVectors,
    int k,
    const unordered_set<pii, pair_hash>& obstacles,
    int width,
    int height
) {
    double captured_fm_k = 0.0;
    double captured_fm_km1_dh = 0.0;

    for (const auto& entry : cellToIndex) {
        const pii& u = entry.first;

        // Get all 8-connected walkable neighbors
        vector<pii> neighbors = get_neighbors_8(u, {height, width}, obstacles);
        for (const pii& v : neighbors) {
            if (!cellToIndex.count(v)) continue;

            // --- FM[k] ---
            const auto& emb_u = fmEmbedding.at(u);
            const auto& emb_v = fmEmbedding.at(v);
            float fm_k_val = 0.0f;
            for (int i = 0; i < k; ++i) {
                fm_k_val += fabs(emb_u[i] - emb_v[i]);  // L1 norm
            }
            captured_fm_k += fm_k_val;

            // --- FM[k−1] + DH[0] ---
            float fm_km1_val = 0.0f;
            for (int i = 0; i < k - 1; ++i) {
                fm_km1_val += fabs(emb_u[i] - emb_v[i]);
            }

            int idx_u = u.first * width + u.second;
            int idx_v = v.first * width + v.second;
            float dh_val = fabs(dhVectors[0][idx_u] - dhVectors[0][idx_v]);

            captured_fm_km1_dh += fm_km1_val + dh_val;
        }
    }

    std::cerr << " Captured Heuristic for FM[" << k << "]: "
              << captured_fm_k << std::endl;

    std::cerr << " Captured Heuristic for FM[" << (k - 1) << "] + DH[0]: "
              << captured_fm_km1_dh << std::endl;
}



void summarize_expansions() {
    if (expansion_log.empty()) return;

    int n = expansion_log.size();

    // Mean
    double mean = accumulate(expansion_log.begin(), expansion_log.end(), 0.0) / n;

    // Median
    vector<int> sorted = expansion_log;
    sort(sorted.begin(), sorted.end());
    double median = (n % 2 == 0)
        ? (sorted[n/2 - 1] + sorted[n/2]) / 2.0
        : sorted[n/2];

    // 95% confidence interval for the mean
    double sum_sq = 0.0;
    for (int x : expansion_log) sum_sq += (x - mean) * (x - mean);
    double stdev = sqrt(sum_sq / (n - 1));
    double ci_95 = 1.96 * stdev / sqrt(n);  // 95% CI assuming normal distribution

    cout << "Mean: " << mean << endl;
    cout << "Median: " << median << endl;
    cout << "95% CI on mean: " << ci_95 << endl;
}

// ===================================================================================
// ---------------- PreprocessMap (optional) ----------------
void PreprocessMap(const vector<bool> &bits, int width, int height, const string &filename) {
    // Additional pre-processing can be done here if needed.
}


long long heuristic_micro_time_global;

// ---------------- PrepareForSearch ----------------
void* PrepareForSearch(const vector<bool>& bits, int width, int height, const string& filename) {
    global_obstacles.clear();
    global_width = width;
    global_height = height;

    // Build obstacles.
    for (int r = 0; r < height; ++r) {
        for (int c = 0; c < width; ++c) {
            int idx = r * width + c;
            if (!bits[idx]) {
                global_obstacles.insert({r, c});
            }
        }
    }

    pair<int, int> grid_size = {height, width};
    vector<Edge> edges;
    buildGraphFromGrid(grid_size, global_obstacles, edges);
    int numFreeCells = cellToIndex.size();

    auto* ctx = new FMContext(); // dynamically allocate context

    auto start_time = std::chrono::high_resolution_clock::now();

    // FM and FM+DH
    int Kmax = 5;           // 10 -default Maximum number of dimensions.
    double epsilon = 1e-2;     // 1e-6 - Threshold for minimal distance.
    vector<vector<double>> embedding = fastMap(edges, numFreeCells, Kmax, epsilon);

    // // L1 norm
    // vector<vector<double>> embedding = fastMapEnhanced(edges, numFreeCells, Kmax, epsilon);
    // // L2 norm
    // vector<vector<double>> embedding = fastMapEnhanced(edges, numFreeCells, Kmax, epsilon, FastMapNorm::L2);
    

    // Populate the global FM embedding mapping.
    fmEmbedding.clear();
    for (int i = 0; i < numFreeCells; i++) {
        fmEmbedding[indexToCell[i]] = embedding[i];
    }

    cellToIndex_E.clear();
    indexToCell_E.clear();
    int id = 0;

    for (int r = 0; r < height; ++r) {
        for (int c = 0; c < width; ++c) {
            int idx = r * width + c;
            if (bits[idx]) { // walkable cell
                cellToIndex_E[{r, c}] = id;
                indexToCell_E.push_back({r, c});
                id++;
            }
        }
    }

    int Kmax_E = 9;           // 10 -default Maximum number of dimensions.
    double epsilon_E = 1e-2;     // 1e-6 - Threshold for minimal distance.
    build_fmEmbedding_HE(width, height, Kmax_E, epsilon_E); // Build the FM embedding with HE pivot selection.

    int numLandmarks_E = 1; // Number of landmarks for DH
    computeDHVectors_E(numLandmarks_E); // Compute DH vectors for the first landmark.

    // DH
    // === Select multiple landmarks for DH ===
    int numLandmarks = 5;  // or 20 or more, based on performance/memory tradeoff
    computeDHVectors(numLandmarks);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    ctx->heuristic_micro_time = duration.count(); 
    heuristic_micro_time_global = ctx->heuristic_micro_time;
    std::cerr << "Heuristic computed in " << ctx->heuristic_micro_time << " ms " << std::endl;


    // // Compute captured heuristic
    // for (int k = 1; k <= 10; ++k) {
    //     computeDHVectors(1);  // always use 1 DH for FM[k−1] + DH
    
    //     vector<vector<double>> embedding = fastMap(edges, numFreeCells, k, epsilon);
    //     fmEmbedding.clear();
    //     for (int i = 0; i < numFreeCells; i++) {
    //         fmEmbedding[indexToCell[i]] = embedding[i];
    //     }
    
    //     computeCapturedHeuristic(
    //         cellToIndex,
    //         fmEmbedding,
    //         dhVectors,
    //         k,
    //         global_obstacles,
    //         global_width,
    //         global_height
    //     );
    // }
    
    return ctx;  
}

// ---------------- GetPath ----------------
bool GetPath(void *data, xyLoc s, xyLoc g, vector<xyLoc> &path) {
    auto* ctx = static_cast<FMContext*>(data);  // Cast the context

    pair<int, int> grid_size = {global_height, global_width};
    vector<pii> raw_path = a_star({s.y, s.x}, {g.y, g.x}, grid_size, global_obstacles);

    path.clear();
    for (const pii& p : raw_path) {
        // Convert from (row, col) to (x, y); note the coordinate swap.
        path.push_back(xyLoc{ static_cast<short>(p.second), static_cast<short>(p.first) });
    }

    return true;
}

// ---------------- GetName ----------------
string GetName() {
    return "FM_A*_8-connected";
}
