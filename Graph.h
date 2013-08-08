// ----------------------
// projects/graph/Graph.h
// Copyright (C) 2013
// Glenn P. Downing
// ----------------------

#ifndef Graph_h
#define Graph_h

// --------
// includes
// --------

#include <cassert> // assert
#include <cstddef> // size_t
#include <utility> // make_pair, pair
#include <vector>  // vector
#include <stack>   // stack
#include <set>     // set
#include "boost/graph/exception.hpp"
// -----
// Graph
// -----

using namespace std;

class Graph {
public:
    // --------
    // typedefs
    // --------
    
    typedef int vertex_descriptor;
    typedef std::pair <vertex_descriptor, vertex_descriptor> edge_descriptor;
    
    typedef std::vector<vertex_descriptor>::const_iterator vertex_iterator;
    typedef std::set<edge_descriptor>::iterator edge_iterator;
    typedef std::vector<vertex_descriptor>::const_iterator adjacency_iterator;
    
    typedef std::size_t vertices_size_type;
    typedef std::size_t edges_size_type;
    
public:
    // --------
    // add_edge
    // --------
    
    /**
     * Add a given edge (v1, v2) to the graph 'gr', if the graph already contains an edge between v1 and v2, then return the existing edge and set the bool false.
     */
    friend std::pair<edge_descriptor, bool> add_edge (vertex_descriptor v1, vertex_descriptor v2, Graph& gr) {
        edge_descriptor edge = make_pair(v1,v2);
        bool edgeExists = gr.edges.insert(edge).second;
        if(edgeExists)
            gr.g[v1].push_back(v2);
        return std::make_pair(edge, edgeExists);}
    
    // ----------
    // add_vertex
    // ----------
    
    /**
     * Adds a vertex to the graph and returns the vertex descriptor for the new vertex.
     */
    friend vertex_descriptor add_vertex (Graph& gr) {
        vector<vertex_descriptor> newVertex;
        gr.g.push_back(newVertex);
        vertex_descriptor v;
        v = gr.g.size()-1;
        gr.vertices.push_back(v);
        return v;}
    
    // -----------------
    // adjacent_vertices
    // -----------------
    
    /**
     * Returns an iterator-range providing access to the vertices adjacent to vertex u in graph g. For example, if u -> v is an edge in the graph, then v will be in this iterator-range.
     */
    friend std::pair<adjacency_iterator, adjacency_iterator> adjacent_vertices (vertex_descriptor u, const Graph& gr) {
        adjacency_iterator b = gr.g[u].begin();
        adjacency_iterator e = gr.g[u].end();
        return std::make_pair(b, e);}
    
    // ----
    // edge
    // ----
    
    /**
     * If an edge from vertex v1 to vertex v2 exists, return a pair containing one such edge and true. If there are no edges between u and v, return a pair with an arbitrary edge descriptor and false.
     */
    friend std::pair<edge_descriptor, bool> edge (vertex_descriptor v1, vertex_descriptor v2, const Graph& gr) {
        bool exists = true;
        edge_descriptor edge = make_pair(v1,v2);
        edge_iterator itr = gr.edges.find(edge);
        if(itr == gr.edges.end()) {
            exists = false;
            return std::make_pair(edge, exists);
        }
        edge = *itr;
        return std::make_pair(edge, exists);}
    
    // -----
    // edges
    // -----
    
    /**
     * Returns an iterator-range providing access to the edge set of graph gr.
     */
    friend std::pair<edge_iterator, edge_iterator> edges (const Graph& gr) {
        edge_iterator b = gr.edges.begin();
        edge_iterator e = gr.edges.end();
        return std::make_pair(b, e);}
    
    // ---------
    // num_edges
    // ---------
    
    /**
     * Returns the number of edges in the graph gr.
     */
    friend edges_size_type num_edges (const Graph& gr) {
        return gr.edges.size();}
    
    // ------------
    // num_vertices
    // ------------
    
    /**
     * Returns the number of vertices in the graph gr.
     */
    friend vertices_size_type num_vertices (const Graph& gr) {
        return gr.vertices.size();}
    
    // ------
    // source
    // ------
    
    /**
     * Returns the source vertex of edge e.
     */
    friend vertex_descriptor source (edge_descriptor edge, const Graph& gr) {
        return edge.first;}
    
    // ------
    // target
    // ------
    
    /**
     * Returns the target vertex of edge e.
     */
    friend vertex_descriptor target (edge_descriptor edge, const Graph& gr) {
        return edge.second;}
    
    // ------
    // vertex
    // ------
    
    /**
     * Returns the nth vertex in the graph's vertex list.
     */
    friend vertex_descriptor vertex (vertices_size_type n, const Graph& gr) {
        return gr.vertices[n];}
    
    // --------
    // vertices
    // --------
    
    /**
     * Returns an iterator-range providing access to the vertex set of graph g.
     */
    friend std::pair<vertex_iterator, vertex_iterator> vertices (const Graph& gr) {
        vertex_iterator b = gr.vertices.begin();
        vertex_iterator e = gr.vertices.end();
        return std::make_pair(b, e);}
    
private:
    // ----
    // data
    // ----
    
    std::vector< std::vector<vertex_descriptor> > g; // something like this
    std::vector<vertex_descriptor> vertices;
    std::set<edge_descriptor> edges;
    // -----
    // valid
    // -----
    
    /**
     * <your documentation>
     */
    bool valid () const {
        return true;}
    
public:
    // ------------
    // constructors
    // ------------
    
    /**
     * <your documentation>
     */
    Graph () {
        assert(valid());}
    
    // Default copy, destructor, and copy assignment
    // Graph  (const Graph<T>&);
    // ~Graph ();
    // Graph& operator = (const Graph&);
};




// ---------
// DFS_visit
// ---------
/**
 * White marks vertices that have yet to be discovered(=0).
 * Gray marks a vertex that is discovered but still has vertices adjacent to it that are undiscovered(=1).
 * A black vertex is discovered vertex that is not adjacent to any white vertices(=2).
 */

template <typename V, typename P = std::vector<V>, typename C = vector<int>, typename G>
bool DFS_visit(V u, P& p, C& color, const G& gr) {
    std::pair<typename G::adjacency_iterator, typename G::adjacency_iterator> adjacenctVerticesIterators= adjacent_vertices(u, gr);
    typename G::adjacency_iterator b = adjacenctVerticesIterators.first;
    typename G::adjacency_iterator e = adjacenctVerticesIterators.second;
    typename G::vertex_descriptor v;
    
    color[u] = 1;
    while (b != e) {
        v = *b;
        if (color[v] == 1) {
            return true;
        }
        else if(color[v] == 0) {
            p[v] = u;
            if(DFS_visit(v, p, color, gr)) {
                return true;
            }
        }
        ++b;
    }
    
    color[u] = 2;
    return false;
}

// ---------
// has_cycle
// ---------
/**
 * depth-first traversal
 * three colors
 * determines where a given graph gr contains a circle or not.
 * returns true if a circle is found, false otherwise.
 */

template <typename G>
bool has_cycle (const G& gr) {
    std::pair<typename G::vertex_iterator, typename G::vertex_iterator> verticesIterators= vertices(gr);
    typename G::vertex_iterator b = verticesIterators.first;
    typename G::vertex_iterator e = verticesIterators.second;
    std::vector<int> color(num_vertices(gr), 0);
    std::vector<typename G::vertex_descriptor> p;//predecessors
    typename G::vertex_descriptor v;
    
    for(unsigned int i = 0; i<num_vertices(gr); i++) {
        p.push_back(i);
    }
    
    while (b!= e) {
        v = *b;
        if(color[v] == 0) {
            if(DFS_visit<typename G::vertex_descriptor, std::vector<typename G::vertex_descriptor>, vector<int>, G>(v, p, color, gr))
                return true;
        }
        ++b;
    }
    return false;}


template <typename V, typename B = std::vector<int>, typename C = std::stack<int>, typename G>
void topological_sort_helper(V u, B& color, C& s, const G& gr) {
    std::pair<typename G::adjacency_iterator, typename G::adjacency_iterator> adjacenctVerticesIterators= adjacent_vertices(u, gr);
    typename G::adjacency_iterator b = adjacenctVerticesIterators.first;
    typename G::adjacency_iterator e = adjacenctVerticesIterators.second;
    typename G::vertex_descriptor v;
    
    color[u] = 1;
    while (b != e) {
        v = *b;
        if (color[v] == 1) {
            throw boost::not_a_dag();
        }
        else if(color[v] == 0) {
            topological_sort_helper(v, color, s, gr);
        }
        ++b;
    }
    
    color[u] = 2;
    s.push(u);
}


// ----------------
// topological_sort
// ----------------

/**
 * depth-first traversal
 * two colors
 * <your documentation>
 * @throws Boost's not_a_dag exception if has_cycle()
 */

template <typename G, typename OI>
void topological_sort (const G& gr, OI x) {
    if(has_cycle(gr))
        throw boost::not_a_dag();
    
    std::pair<typename G::vertex_iterator, typename G::vertex_iterator> verticesIterators= vertices(gr);
    typename G::vertex_iterator b = verticesIterators.first;
    typename G::vertex_iterator e = verticesIterators.second;
    std::vector<int> color(num_vertices(gr), 0);
    std::stack<int> s;
    typename G::vertex_descriptor v;
    
    while(b!=e) {
        v = *b;
        if(color[v] == 0) {
            topological_sort_helper<typename G::vertex_descriptor, std::vector<int>, std::stack<int>, G> (v, color, s, gr);
        }
        ++b;
    }
    while(!s.empty()) {
        *x = s.top();
        ++x;
        s.pop();
    }
}
#endif // Graph_h


