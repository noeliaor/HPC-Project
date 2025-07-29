#ifndef IO_UTILS_HPP
#define IO_UTILS_HPP

/// Factor umbral relativo usado para filtrar contactos en la deteccion de cadenas.
/// threshold = UMBRAL_FACTOR * fuerza_promedio
#define UMBRAL_FACTOR 0.5

#include <string>
#include <vector>
#include "contact.hpp"

/**
 *  Lista todos los archivos con extension .dat en un directorio.
 */
std::vector<std::string> obtener_archivos_dat(const std::string& dir);

/**
 * Extrae el timestep desde el nombre de un archivo.
 */
int extraer_timestep(const std::string& fname);

/**
 * Lee los contactos desde un archivo en formato ESyS-Particle RAW WITH POS ID.
 * Usa la funcion `load_contacts()` implementada en utils.cpp.
 */
std::vector<Contact> leer_contactos(const std::string& path);

/**
 * Detecta componentes conexas (cadenas) usando fuerza relativa como umbral.
 * Construye un grafo a partir de los contactos que superan el umbral y luego identifica las componentes conexas.
 */
std::vector<std::vector<int> >
detectar_componentes(const std::vector<Contact>& contactos, int rank, int np);

/**
 * Une las componentes detectadas por todos los procesos en el proceso maestro (rank 0). 
 */
std::vector<std::vector<int> >
unir_componentes(const std::vector<std::vector<int> >& locales, int rank, int np);

#endif
