// utils.cpp  (C++11)
// Funciones utilitarias: parsing de contactos, centroides, metricas, MPI dtype
// ============================================================================

#include "utils.hpp"
#include "contact.hpp"   // struct Contact
#include <mpi.h>

#include <cmath>
#include <cctype>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <array>
#include <unordered_set>
/**
* Crear un tipo de dato personalizado para MPI que permite enviar y recibir
* estructuras Contact como un bloque contiguo de bytes. Se usa MPI_BYTE
* para evitar definir manualmente el layout.
*/ 
MPI_Datatype create_mpi_contact_type() {
    MPI_Datatype t;
    MPI_Type_contiguous(sizeof(Contact), MPI_BYTE, &t);
    MPI_Type_commit(&t);
    return t;
}

/** Calcula el centroide de una cadena de fuerza a partir de los IDs de
* las particulas que la componen. El centroide se estima como el promedio
* de los puntos medios de los contactos que involucran alguno de esos IDs.
*/ 

std::array<double,3> calcular_centroide(
    const std::vector<int>& ids,
    const std::vector<Contact>& contacts
) {
    std::array<double,3> c = {{0.0, 0.0, 0.0}};
    if (ids.empty() || contacts.empty()) return c;

    std::unordered_set<int> idset(ids.begin(), ids.end());

    int n = 0;
    for (std::size_t k = 0; k < contacts.size(); ++k) {
        const Contact& ct = contacts[k];
        if (idset.count(ct.id1) || idset.count(ct.id2)) {
            c[0] += 0.5 * (ct.x1 + ct.x2);
            c[1] += 0.5 * (ct.y1 + ct.y2);
            c[2] += 0.5 * (ct.z1 + ct.z2);
            ++n;
        }
    }

    if (n > 0) {
        c[0] /= n;
        c[1] /= n;
        c[2] /= n;
    }
    return c;
}

/** Calcula la distancia euclidea entre dos centroides en 3D.
*/
double calcular_distancia(
    const std::array<double,3>& a,
    const std::array<double,3>& b
) {
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    double dz = a[2] - b[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

/** Calcula el indice de Jaccard entre dos cadenas de fuerza representadas
*como vectores de IDs de particulas.
*/
double calcular_jaccard(
    const std::vector<int>& A,
    const std::vector<int>& B
) {
    if (A.empty() && B.empty()) return 1.0;

    std::unordered_set<int> s(A.begin(), A.end());
    int inter = 0;
    for (std::size_t i = 0; i < B.size(); ++i) {
        if (s.count(B[i])) ++inter;
    }
    int uni = static_cast<int>(A.size() + B.size() - inter);
    return (uni > 0) ? double(inter) / double(uni) : 0.0;
}

/** Calcular la fuerza media (magnitud) de todos los contactos en un vector.
*/ 
double compute_average_force(const std::vector<Contact>& contacts) {
    if (contacts.empty()) return 0.0;
    double sum = 0.0;
    for (std::size_t k = 0; k < contacts.size(); ++k)
        sum += contacts[k].fuerza();
    return sum / contacts.size();
}

/**
*  PARSER DE ARCHIVOS partForce.*.dat  (12 o 14 columnas)
* Formatos segun manual de ESyS-Particle (RAW WITH POS ID):
* Escalar (12 columnas):
*   id1 id2 p1x p1y p1z p2x p2y p2z ipx ipy ipz v

* Vectorial (14 columnas):
* id1 id2 p1x p1y p1z p2x p2y p2z ipx ipy ipz fx fy fz
*/
std::vector<Contact> load_contacts(const std::string& filename)
{
    std::vector<Contact> contacts;

    std::ifstream f(filename.c_str());
    if (!f.is_open()) {
        std::cerr << "[load_contacts] No se pudo abrir: " << filename << "\n";
        return contacts;
    }

    std::string line;
    std::size_t lineno = 0;
    std::size_t good = 0;
    std::size_t bad  = 0;

    while (std::getline(f, line)) {
        ++lineno;

        // Saltar líneas vacías o comentarios
        std::size_t p = 0;
        while (p < line.size() && std::isspace(static_cast<unsigned char>(line[p]))) ++p;
        if (p == line.size()) continue;
        if (line[p] == '#') continue;

        // Parsear tokens
        std::istringstream iss(line);
        std::vector<double> vals;
        vals.reserve(16);
        double tmp;
        while (iss >> tmp) vals.push_back(tmp);

        if (vals.size() < 12) {
            ++bad;
            continue;
        }

        Contact c;
        c.id1 = static_cast<int>(vals[0]);
        c.id2 = static_cast<int>(vals[1]);
        c.x1 = vals[2];  c.y1 = vals[3];  c.z1 = vals[4];
        c.x2 = vals[5];  c.y2 = vals[6];  c.z2 = vals[7];

        if (vals.size() >= 14) {
            // Formato vectorial
            c.fx = vals[11];
            c.fy = vals[12];
            c.fz = vals[13];
        } else {
            // Formato escalar
            c.fx = vals[11];
            c.fy = 0.0;
            c.fz = 0.0;
        }

        contacts.push_back(c);
        ++good;
    }

#if 0
    std::cerr << "[load_contacts] " << filename
              << " good=" << good << " bad=" << bad
              << " totalLines=" << lineno << "\n";
#endif

    return contacts;
}
