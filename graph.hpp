// graph.hpp

/**
 * Declaraciones de funciones y tipos para la construccion y analisis del grafo de contactos.
 * Proporciona utilidades para construir un grafo a partir de contactos entre particulas
 * y detectar componentes conexas (cadenas de fuerza).
 */

#pragma once
#include "contact.hpp"
#include <vector>
#include <unordered_map>
#include <utility>

/**
 * Representa un grafo no dirigido mediante un mapa de listas de adyacencia,
 * donde cada nodo tiene una lista de vecinos conectados.
 */
using Graph = std::unordered_map<int, std::vector<int>>;

/**
 * Construye un grafo de contactos a partir de los datos de fuerza.
 */
Graph build_graph(const std::vector<Contact>& contacts, double threshold);

/**
 * Encuentra las componentes conexas del grafo.
 * Cada componente representa una cadena de fuerza dentro del sistema granular.
 */
std::vector<std::vector<int>> find_connected_components(const Graph& graph);

/**
 *Extrae todas las aristas unicas del grafo.
 */
std::vector<std::pair<int,int>> get_edges(const Graph& g);
