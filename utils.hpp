#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <array>
#include <mpi.h>
#include "contact.hpp"

/**
 * Crea un tipo de dato MPI para transmitir objetos Contact.
 * Este tipo se utiliza para enviar/recibir vectores de Contact usando MPI,
 * tratándolos como bloques contiguos de bytes.
 */
MPI_Datatype create_mpi_contact_type();

/**
 * Calcula el centroide de una cadena de particulas a partir de sus contactos.
 */
std::array<double,3> calcular_centroide(const std::vector<int>& ids,
                                        const std::vector<Contact>& contactos);

/**
 * Calcula la distancia euclidea entre dos puntos 3D.
 */
double calcular_distancia(const std::array<double,3>& a,
                          const std::array<double,3>& b);

/**
 * Calcula el índice de Jaccard entre dos conjuntos de IDs.
 */
double calcular_jaccard(const std::vector<int>& A,
                        const std::vector<int>& B);

/**
 * Calcula la fuerza media de un conjunto de contactos.
 */
double compute_average_force(const std::vector<Contact>& contacts);

/**
 * Carga desde archivo una lista de contactos entre particulas.
 */
std::vector<Contact> load_contacts(const std::string& filename);

#endif
