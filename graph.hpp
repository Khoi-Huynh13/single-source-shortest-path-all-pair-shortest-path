#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <iostream>
#include <fstream>
#include <utility>
#include <functional>
#include <vector>
#include <string>
#include <queue>
#include <unordered_map>
#include <limits>

template <typename T>
class Graph {
 private:
  std::vector<std::unordered_map<int, T> > adjList {};
  int numVertices {};

 public:
  // empty graph with N vertices
  explicit Graph(int N);

  // construct graph from edge list in filename
  explicit Graph(const std::string& filename);

  // add an edge directed from vertex i to vertex j with given weight
  void addEdge(int i, int j, T weight);

  // removes edge from vertex i to vertex j
  void removeEdge(int i, int j);

  // is there an edge from vertex i to vertex j?
  bool isEdge(int i, int j) const;

  // return weight of edge from i to j
  // will throw an exception if there is no edge from i to j
  T getEdgeWeight(int i, int j) const;

  // returns number of vertices in the graph
  int size() const;

  // return iterator to a particular vertex
  const std::unordered_map<int, T>& neighbours(int a) const {
    return adjList.at(a);
  }

  // Add the temporary vertext in johnsonAPSP
  void addTempVertex(); 

  // Remove the temporary vertex in johnsonAPSP
  void removeTempVertex();

  // Copy assignment operator
  Graph& operator=(Graph other);

  // Copy constructor
  Graph(const Graph& other) {
    numVertices = other.numVertices;
    adjList = other.adjList;
  }
};

template <typename T>
Graph<T>& Graph<T>::operator=(Graph other) {
  std::swap(numVertices, other.numVertices);
  std::swap(adjList, other.adjList);
  return *this;
}

template <typename T>
void Graph<T>::addTempVertex() {
  numVertices++;
  adjList.resize(numVertices);
}

template <typename T>
void Graph<T>::removeTempVertex() {
  numVertices--;
  adjList.resize(numVertices);
}

template <typename T>
Graph<T>::Graph(int N) : adjList(N), numVertices {N} {}

template <typename T>
Graph<T>::Graph(const std::string& inputFile) {
  std::ifstream infile {inputFile};
  if (!infile) {
    std::cerr << inputFile << " could not be opened\n";
    return;
  }
  // first line has number of vertices
  infile >> numVertices;
  adjList.resize(numVertices);
  int i {};
  int j {};
  double weight {};
  // assume each remaining line is of form
  // origin dest weight
  while (infile >> i >> j >> weight) {
    addEdge(i, j, static_cast<T>(weight));
  }
}

template <typename T>
int Graph<T>::size() const {
  return numVertices;
}

template <typename T>
void Graph<T>::addEdge(int i, int j, T weight) {
  if (i < 0 or i >= numVertices or j < 0 or j >= numVertices) {
    throw std::out_of_range("invalid vertex number");
  }
  adjList[i].insert({j, weight});
}

template <typename T>
void Graph<T>::removeEdge(int i, int j) {
  // check if i and j are valid
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    adjList[i].erase(j);
  }
}

template <typename T>
bool Graph<T>::isEdge(int i, int j) const {
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    return adjList.at(i).contains(j);
  }
  return false;
}

template <typename T>
T Graph<T>::getEdgeWeight(int i, int j) const {
  return adjList.at(i).at(j);
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Graph<T>& G) {
  for (int i = 0; i < G.size(); ++i) {
    out << i << ':';
    for (const auto& [neighbour, weight] : G.neighbours(i)) {
      out << " (" << i << ", " << neighbour << ")[" << weight << ']';
    }
    out << '\n';
  }
  return out;
}


// APSP functions
// Use this function to return an "infinity" value
// appropriate for the type T
template <typename T>
T infinity() {
  if (std::numeric_limits<T>::has_infinity) {
    return std::numeric_limits<T>::infinity();
  } else {
    return std::numeric_limits<T>::max();
  }
}

// implement an algorithm for determining if G
// has a negative weight cycle here
template <typename T>
bool existsNegativeCycle(const Graph<T>& G) {
  std::vector<T> bestDistanceTo(G.size(), 0); // A zero vector instead of a vector where source = 0 and every other vertex is infinity to achieve the same affect as adding an additional vertex to the graph;
  T new_distance;

  for (int i = 0; i < G.size(); i++) { 
    for (int start_vertex = 0; start_vertex < G.size(); start_vertex++) { // Relax the edges G.size() times and not G.size() - 1 like normal since we want to detect negative weight cycle
      for (const auto& [neighbour, weight] : G.neighbours(start_vertex)) {
        new_distance = bestDistanceTo.at(start_vertex) + G.getEdgeWeight(start_vertex, neighbour);
        if (bestDistanceTo.at(neighbour) > new_distance) {
          if (i == G.size() - 1) { // If still able to relax an edge at the G.size() iteration then there is a negative weighted cycle
            return true;
          } else {
            bestDistanceTo.at(neighbour) = new_distance;
          }
        }
      }
    }
  }

  return false;
}

// implement Johnson's APSP algorithm here
template <typename T>
std::vector<std::vector<T>> floydWarshallAPSP(const Graph<T>& G) {
  std::vector<std::vector<T>> distanceMatrix(G.size(), std::vector<T>(G.size(), infinity<T>()));
  T new_distance;

  // Set up distanceMatrix for best distance between start_vertex and end_vertex WITHOUT any intermediate vertices
  for (int start_vertex = 0; start_vertex < G.size(); start_vertex++) {
    for (int end_vertex = 0; end_vertex < G.size(); end_vertex++) {
      if (start_vertex == end_vertex) {
        distanceMatrix[start_vertex][end_vertex] = T {};
      } else if (G.isEdge(start_vertex, end_vertex)) {
        distanceMatrix[start_vertex][end_vertex] = G.getEdgeWeight(start_vertex, end_vertex);
      }
    }
  }

  // Iteratively improve on the distances by ocnsidering intermediate vertices
  for (int intermediate_vertex = 0; intermediate_vertex < G.size(); intermediate_vertex++) {
    for (int start_vertex = 0; start_vertex < G.size(); start_vertex++) {
      for (int end_vertex = 0; end_vertex < G.size(); end_vertex++) {
        if (distanceMatrix[start_vertex][intermediate_vertex] != infinity<T>() && distanceMatrix[intermediate_vertex][end_vertex] != infinity<T>()) { // Check if the intermediate vertices bridge between start vertex to end vertex
          new_distance = distanceMatrix[start_vertex][intermediate_vertex] + distanceMatrix[intermediate_vertex][end_vertex];
          if (distanceMatrix[start_vertex][end_vertex] > new_distance) { // If path through intermediate vertices is better than current path, update it
            distanceMatrix[start_vertex][end_vertex] = new_distance;
          }
        }
      }
    }
  }

  return distanceMatrix;
}

// implement the Floyd-Warshall APSP algorithm here
template <typename T>
std::vector<std::vector<T>> johnsonAPSP(const Graph<T>& G) {
  std::vector<std::vector<T>> result(G.size());
  int s = G.size();
  Graph<T> G_new = G; // Set up modified graph, introducing dummy vertex and connecting dummy vertex with all other vertices with edge weight of 0
  G_new.addTempVertex();

  for (int vertex = 0; vertex < G_new.size() - 1; vertex++) {
    G_new.addEdge(s, vertex, 0);
  } 
 
  // Find shortest path lenghs from dummy vertex to every other vertices
  std::vector<T> h_v = bellman_ford(G_new, s);
  T original_weight {};

  // Modify the edge weights using h(v)
  for (int start_vertex = 0; start_vertex < G_new.size() - 1; start_vertex++) {
    for (const auto& [neighbour, weight] : G.neighbours(start_vertex)) {
      original_weight = G_new.getEdgeWeight(start_vertex, neighbour);
      G_new.removeEdge(start_vertex, neighbour);
      G_new.addEdge(start_vertex, neighbour, original_weight + h_v[start_vertex] - h_v[neighbour]); // Update edge weight with w(u, v) + h(u) + h(v)
    }
  }

  G_new.removeTempVertex();
  std::vector<T> bestDistanceTo(G_new.size());

  // Run Dijkstra starting from each vertex
  for (int start_vertex = 0; start_vertex < G_new.size(); start_vertex++) {
    bestDistanceTo = dijkstra(G_new, start_vertex);
    for (int end_vertex = 0; end_vertex < G_new.size(); end_vertex++) {
      if (bestDistanceTo[end_vertex] != infinity<T>()) { // Only able to revert weight of path that is reachable
        bestDistanceTo[end_vertex] += h_v[end_vertex] - h_v[start_vertex]; // Transform the weight of the path back to corresponding original weight
      }
    }
    result[start_vertex] = bestDistanceTo;
  }

  return result;
}

template <typename T>
std::vector<T> dijkstra(const Graph<T>& G, int source) {
  using DistAndVertex = std::pair<T, int>;
  using minPQ = std::priority_queue<DistAndVertex, std::vector<DistAndVertex>, std::greater<DistAndVertex>>;
  // Basic set up
  minPQ queue {};
  queue.push({0.0, source});
  std::vector<T> bestDistanceTo(G.size(), infinity<T>());
  bestDistanceTo.at(source) = 0;
  std::vector<bool> visited(G.size());

  while (!queue.empty()) {
    auto [distance, current] = queue.top(); // Keep track of currently processing vertex
    queue.pop();

    if (visited[current] == true) { // If vertex has already been processed, skip
      continue;
    } else {
      visited[current] = true;

      for (const auto& [neighbour, weight] : G.neighbours(current)) {
        if (bestDistanceTo[neighbour] > bestDistanceTo[current] + weight) { // Update best path from start vertex to end vertex
          bestDistanceTo[neighbour] = bestDistanceTo[current] + weight;
          queue.push({bestDistanceTo[neighbour], neighbour}); // Keep track of TO BE processed vertex, choose the TO BE processed vertex with the smallest distance from start vertex using priority queue
        }
      }
    }
  }

  return bestDistanceTo;
}

template <typename T>
std::vector<T> bellman_ford(const Graph<T>& G, int source) {
  std::vector<T> bestDistanceTo(G.size(), infinity<T>()); // Basic set up
  bestDistanceTo.at(source) = 0;
  T new_distance;

  for (int i = 0; i < G.size() - 1; i++) { // Only relax the edges G.size() - 1 times because we don't care about negative weighted cycle here
    for (int start_vertex = 0; start_vertex < G.size(); start_vertex++) {
      if (bestDistanceTo[start_vertex] != infinity<T>()) { // Only able to relax edges from vertex with defined best distance
        for (const auto& [neighbour, weight] : G.neighbours(start_vertex)) {
          new_distance = bestDistanceTo[start_vertex] + G.getEdgeWeight(start_vertex, neighbour);
          if (bestDistanceTo[neighbour] > new_distance) { // Update best path from start vertex to end vertex 
            bestDistanceTo[neighbour] = new_distance;
          }
        }
      }
    }
  }

  return bestDistanceTo;
}

#endif      // GRAPH_HPP_
