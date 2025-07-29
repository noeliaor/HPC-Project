// io_utils.cpp --------------------------------------------------
// Utilidades de E/S y wrappers para detección de componentes.

#include "io_utils.hpp" 
#include "graph.hpp"  
#include "contact.hpp"   
#include "utils.hpp"      

#include <dirent.h>
#include <mpi.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/**
 * Lista los archivos con extension .dat en el directorio especificado.
 */
std::vector<std::string> obtener_archivos_dat(const std::string& dir)
{
    std::vector<std::string> out;
    DIR* d = opendir(dir.c_str());
    if (!d) return out;

    struct dirent* e;
    while ((e = readdir(d))) {
        std::string n(e->d_name);
        if (n.size() > 4 && n.substr(n.size() - 4) == ".dat")
            out.push_back(n);
    }
    closedir(d);
    return out;
}

/**
 * @brief Extrae el nro. de timestep desde el nombre del archivo.
 * Busca y concatena todos los digitos presentes en el nombre.
 */
int extraer_timestep(const std::string& fname)
{
    std::string num;
    num.reserve(fname.size());
    for (size_t i = 0; i < fname.size(); ++i) {
        unsigned char uc = static_cast<unsigned char>(fname[i]);
        if (std::isdigit(uc)) num.push_back(static_cast<char>(uc));
    }
    return num.empty() ? -1 : std::stoi(num);
}

/**
 * Lee los contactos desde un archivo .dat (formato ESyS RAW WITH POS ID).
 * Usa la función load_contacts definida en utils.cpp e imprime
 * advertencias si el archivo esta vacio o es invalido.
 */
std::vector<Contact> leer_contactos(const std::string& path)
{
    std::vector<Contact> v = load_contacts(path);

#ifdef IOUTILS_VERBOSE_LOAD
    std::cerr << "[leer_contactos] " << path
              << " -> " << v.size() << " contactos.\n";
#endif

    if (v.empty()) {
        std::cerr << "[io_utils] Archivo vacío o malformado: " << path << "\n";
    }
    return v;
}

/**
 * Detecta componentes conexas usando un umbral relativo a la fuerza media.
 * Calcula la fuerza promedio y define un umbral como thr = UMBRAL_FACTOR * fuerza_media.
 * Construye un grafo con los contactos que superan el umbral y detecta las
 * componentes conexas (cadenas de fuerza).
 * Se genera   un archivo de log por rank en results/debug_rankX.log.
 */
std::vector<std::vector<int>> detectar_componentes(
    const std::vector<Contact>& cts,
    int rank,
    int /*np*/)
{
    std::vector<std::vector<int>> empty_ret;
    if (cts.empty()) return empty_ret;

    double sum = 0.0;
    for (size_t k = 0; k < cts.size(); ++k)
        sum += cts[k].fuerza();
    double avg = sum / cts.size();

    double thr = UMBRAL_FACTOR * avg;

    // Conteo para debug
    int above = 0;
    for (size_t k = 0; k < cts.size(); ++k)
        if (cts[k].fuerza() >= thr) ++above;

    // Logging de depuracion por proceso
    {
        std::ostringstream oss;
        oss << "results/debug_rank" << rank << ".log";
        std::ofstream dbg(oss.str().c_str(), std::ios::app);
        if (dbg) {
            dbg << "call_detectar_componentes "
                << "avgF=" << avg
                << " thr=" << thr
                << " aboveThr=" << above
                << " total=" << cts.size()
                << " frac=" << (cts.empty() ? 0.0 : double(above) / cts.size())
                << "\n";
        }
    }

    Graph g = build_graph(cts, thr);
    return find_connected_components(g);
}

/**
 * Une los resultados de componentes detectadas en todos los procesos MPI.
 * El proceso 0 recolecta los resultados de los otros procesos, los combina
 * y retorna el vector global de componentes. Los demas procesos solo envian
 * sus datos y retornan un vector vacio.
 */
std::vector<std::vector<int>> unir_componentes(
    const std::vector<std::vector<int>>& local,
    int rank,
    int np)
{
    if (np == 1) return local;

    if (rank == 0) {
        std::vector<std::vector<int>> global = local;
        for (int src = 1; src < np; ++src) {
            int ncomp = 0;
            MPI_Recv(&ncomp, 1, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < ncomp; ++i) {
                int sz = 0;
                MPI_Recv(&sz, 1, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::vector<int> buf(sz);
                if (sz > 0)
                    MPI_Recv(buf.data(), sz, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                global.push_back(buf);
            }
        }
        return global;
    } else {
        int ncomp = static_cast<int>(local.size());
        MPI_Send(&ncomp, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        for (size_t k = 0; k < local.size(); ++k) {
            int sz = static_cast<int>(local[k].size());
            MPI_Send(&sz, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            if (sz > 0) {
                MPI_Send(const_cast<int*>(local[k].data()), sz, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
        }
        return std::vector<std::vector<int>>();
    }
}

