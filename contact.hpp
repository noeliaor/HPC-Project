// contact.hpp
#pragma once
#include <cmath>

/**
 * Estructura que representa un contacto entre 2 particulas de la simulacion granular.
 * Almacena la info. necesaria para construir el grafo de contactos
 * y calcular las fuerzas involucradas.
 */
struct Contact {
    // Identificadores unicos de las 2 particulas en contacto
    int id1, id2;

    // Posiciones de las particulas involucradas en el contacto.
    double x1, y1, z1; 
    double x2, y2, z2;

    // Componentes de la fuerza ejercida en el contacto
    double fx, fy, fz;

    /**
     * Calcula la magnitud de la fuerza de contacto.
     *  Valor escalar de la magnitud de la fuerza.
     */
    double fuerza() const {
        return std::sqrt(fx * fx + fy * fy + fz * fz);
    }
};
