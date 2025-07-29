// graph.cpp

/**
 * Implementación de funciones para construir y analizar el grafo de contactos.
 * Permite construir un grafo no dirigido a partir de contactos entre particulas
 * (con una fuerza superior a cierto umbral) y detectar las componentes conexas, que equivalen
 * a las cadenas de fuerza en el sistema granular.
 */

#include "graph.hpp"
#include <unordered_set>
#include <queue>

/**
 * Construye un grafo a partir de una lista de contactos.
 * Agrega una arista entre dos partículas si la magnitud de la fuerza del contacto supera
 * el umbral especificado.
 */
Graph build_graph(const std::vector<Contact>& contacts, double threshold) {
    Graph g;
    for (auto& c : contacts) {
        if (c.fuerza() >= threshold) {
            g[c.id1].push_back(c.id2);
            g[c.id2].push_back(c.id1);
        }
    }
    return g;
}

/**
 * Encuentra las componentes conexas del grafo.
 * Recorre el grafo utilizando BFS para detectar conjuntos de partuculas conectadas entre si,
 * que corresponden a las cadenas de fuerza en la simulacion.
 */
std::vector<std::vector<int>> find_connected_components(const Graph& graph) {
    std::unordered_set<int> vis;
    std::vector<std::vector<int>> comps;

    for (auto& kv : graph) {
        int u = kv.first;
        if (vis.count(u)) continue;

        std::vector<int> comp;
        std::queue<int> q;
        q.push(u);
        vis.insert(u);

        while (!q.empty()) {
            int x = q.front(); q.pop();
            comp.push_back(x);

            for (int v : graph.at(x)) {
                if (!vis.count(v)) {
                    vis.insert(v);
                    q.push(v);
                }
            }
        }
        comps.push_back(std::move(comp));
    }
    return comps;
}

/**
 * Extrae todas las aristas unicas del grafo.
 * Para evitar duplicados, se incluyen solo aquellas aristas donde id1 < id2.
 */
std::vector<std::pair<int,int>> get_edges(const Graph& g) {
    std::vector<std::pair<int,int>> out;
    for (auto& kv : g) {
        int u = kv.first;
        for (int v : kv.second) {
            if (u < v) {
                out.emplace_back(u, v);
            }
        }
    }
    return out;
}
