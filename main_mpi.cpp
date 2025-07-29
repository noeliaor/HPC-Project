//  detectar_cadenas_mpi  – master–worker (colectivo) + tracking global
//  Procesamiento paralelo de archivos de contacto procedentes de simulaciones DEM.
// Este codigo:
//  - Distribuye archivos entre procesos MPI
//  - Detecta componentes conexas por archivo
//  - Realiza el tracking local y global de cadenas entre timesteps
//  - Usa /tmp si hay espacio, para mejorar I/O
//  - Agrega metricas de rendimiento y resumen global
//  - Concatena archivos de tracking
//
//  Estructura del codigo:
//  1. Inicializacion de MPI y argumentos
//  2. Distribucion de archivos entre procesos
//  3. Procesamiento de archivos por rank
//     - Lectura, deteccion, escritura parcial, tracking local
//  4. Intercambio de frontera para tracking entre procesos consecutivos
//  5. Reduccion y escritura de metricas globales
//  6. Escritura del archivo de tracking global

#include <mpi.h>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/statvfs.h>
#include <limits>

#include "contact.hpp"
#include "utils.hpp"
#include "graph.hpp"
#include "io_utils.hpp"

// rutas dinamicas por corrida
static std::string kDirBase;
static std::string kDirTimings;
static std::string kDirChains;
static std::string kDirResults;

// compara ficheros por timestep
struct CmpTS {
    bool operator()(const std::string &a, const std::string &b) const {
        return extraer_timestep(a) < extraer_timestep(b);
    }
};

// Funcion para obtener espacio disponible en /tmp (en bytes)
size_t espacio_tmp() {
    const char* tmp = std::getenv("SLURM_JOB_ID");
    if (!tmp) return 0;
    struct statvfs buf;
    std::string path = "/tmp/" + std::string(tmp);
    if (statvfs(path.c_str(), &buf) == 0)
        return buf.f_bavail * buf.f_bsize;
    return 0;
}

// Suma el tamaño de los archivos que le toca procesar a este rank
size_t suma_archivos(const std::vector<std::string>& files, const std::string& folder) {
    size_t total = 0;
    for (const auto& f : files) {
        std::ifstream in((folder + "/" + f).c_str(), std::ios::binary | std::ios::ate);
        if (in) total += in.tellg();
        in.close();
    }
    return total;
}

// Para salida agrupada: agregar encabezado si es el primer write
void write_chain_block(std::ofstream &out, int ts, const std::vector<std::vector<int>> &comps) {
    out << "timestep: " << ts << '\n';
    int idx = 1;
    for (auto &c : comps) {
        out << idx++;
        for (int n : c) out << ' ' << n;
        out << '\n';
    }
}

void write_resultado_block(std::ofstream &out, int ts,
    const std::vector<std::vector<int>> &comps,
    const std::vector<Contact> &contacts)
{
    std::unordered_map<int, int> nodo2cadena;
    for (size_t i = 0; i < comps.size(); ++i)
        for (int n : comps[i]) nodo2cadena[n] = i;
    struct Acc { int n; double fx, fy, fz, x, y, z; };
    std::unordered_map<int, Acc> acc;
    for (auto &c : contacts) {
        int cid = -1;
        auto it1 = nodo2cadena.find(c.id1);
        auto it2 = nodo2cadena.find(c.id2);
        if (it1 != nodo2cadena.end()) cid = it1->second;
        else if (it2 != nodo2cadena.end()) cid = it2->second;
        if (cid < 0) continue;
        Acc &a = acc[cid]; if (a.n == 0) a = {0, 0, 0, 0, 0, 0, 0};
        a.n++; a.fx += c.fx; a.fy += c.fy; a.fz += c.fz;
        a.x += 0.5 * (c.x1 + c.x2);
        a.y += 0.5 * (c.y1 + c.y2);
        a.z += 0.5 * (c.z1 + c.z2);
    }
    for (size_t id = 0; id < comps.size(); ++id) {
        auto &a = acc[id]; if (a.n == 0) continue;
        double fmean = std::sqrt(a.fx*a.fx + a.fy*a.fy + a.fz*a.fz) / a.n;
        out << ts << ',' << id + 1 << ',' << a.n << ',' << fmean << ','
            << a.x / a.n << ',' << a.y / a.n << ',' << a.z / a.n << '\n';
    }
}

// tracking local (por rank) y frontera
void save_tracking_file(const std::string &fname, int ts_prev, int ts,
    const std::vector<std::vector<int>> &prev, const std::vector<std::vector<int>> &curr)
{
    std::ofstream trak(fname, std::ios::app);
    if (!trak) {
        std::cerr << "[ERROR] No se pudo abrir archivo de tracking para escritura: " << fname << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 10);
    }
    if (trak.tellp() == 0) {
        trak << "timestep_prev,cadena_prev,timestep_curr,cadena_curr,inter,jaccard\n";
    }
    std::unordered_set<int> used;
    for (size_t i = 0; i < curr.size(); ++i) {
        int best = -1, inter_max = 0;
        for (size_t j = 0; j < prev.size(); ++j) {
            if (used.count(j)) continue;
            int cnt = 0;
            for (int n : curr[i])
                if (std::find(prev[j].begin(), prev[j].end(), n) != prev[j].end()) cnt++;
            if (cnt > inter_max) { inter_max = cnt; best = j; }
        }
        double jac = best >= 0 && best < int(prev.size())
            ? double(inter_max) / (curr[i].size() + prev[best].size() - inter_max)
            : 0.0;
        trak << ts_prev << ',' << best << ',' << ts << ',' << i << ','
             << inter_max << ',' << jac << '\n';
        used.insert(best);
    }
}

// reconstruye el vector de componentes
static std::vector<std::vector<int>> rebuild_components(
    const std::vector<int> &sizes,
    const std::vector<int> &flat)
{
    std::vector<std::vector<int>> comps;
    int pos = 0;
    for (size_t i = 0; i < sizes.size(); ++i) {
        comps.emplace_back(flat.begin() + pos, flat.begin() + pos + sizes[i]);
        pos += sizes[i];
    }
    return comps;
}

// serializa componentes para frontera (solo el primer/ultimo timestep local)
std::vector<std::vector<int>> serializar_componentes(const std::string &fname, int ts)
{
    std::vector<std::vector<int>> comps;
    std::ifstream in(fname);
    if (!in) {
        std::cerr << "[ERROR] No se pudo abrir archivo de componentes para frontera: " << fname << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 20);
    }
    std::string line;
    bool in_ts = false;
    while (std::getline(in, line)) {
        if (line.find("timestep:") != std::string::npos) {
            int t = -1;
            std::istringstream(line.substr(9)) >> t;
            in_ts = (t == ts);
            continue;
        }
        if (!in_ts || line.empty() || !std::isdigit(line[0])) continue;
        std::istringstream ss(line);
        int dummy; ss >> dummy;
        std::vector<int> comp;
        int id;
        while (ss >> id) comp.push_back(id);
        comps.push_back(comp);
    }
    in.close();
    return comps;
}

std::string get_tmp_dir(int rank) {
    const char* tmp = std::getenv("SLURM_JOB_ID");
    struct stat st;
    if (tmp && stat(("/tmp/" + std::string(tmp)).c_str(), &st) == 0) {
        return "/tmp/" + std::string(tmp) + "/proc_" + std::to_string(rank);
    }
    return "";
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    if (argc < 3) {
        if (rank == 0)
            std::cerr << "Uso: mpirun -np N ./detectar_cadenas_mpi <carpeta_data> <nro_repeticion>\n";
        MPI_Finalize(); return 1;
    }
    std::string folder = argv[1];
    int rep = std::stoi(argv[2]);

    kDirBase    = "results/np" + std::to_string(size) + "_rep" + std::to_string(rep);
    kDirTimings = kDirBase + "/timings";
    kDirChains  = kDirBase + "/chains";
    kDirResults = kDirBase + "/results";

    std::string tracking_fname;
    std::string tmpdir = get_tmp_dir(rank);
    if (rank == 0) {
        if (std::system(("mkdir -p " + kDirChains + " " + kDirResults + " " + kDirTimings).c_str()) != 0) {
            std::cerr << "[ERROR] No se pudo crear directorios de salida." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 100);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (!tmpdir.empty()) {
        if (std::system(("mkdir -p " + tmpdir).c_str()) != 0) {
            std::cerr << "[ERROR] Rank " << rank << " no pudo crear directorio en /tmp." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 101);
        }
        tracking_fname = tmpdir + "/tracking_rank" + std::to_string(rank) + ".txt";
    } else {
        tracking_fname = kDirBase + "/tracking_rank" + std::to_string(rank) + ".txt";
    }

    std::vector<std::string> all_files = obtener_archivos_dat(folder);
    std::sort(all_files.begin(), all_files.end(), CmpTS());
    int F = all_files.size();
    int base = F / size, rem = F % size;
    int my_count = base + (rank < rem);
    int offset = rank * base + std::min(rank, rem);
    std::vector<std::string> my_files(
        all_files.begin() + offset,
        all_files.begin() + offset + my_count
    );

    bool usar_tmp = false;
    size_t bytes_needed = suma_archivos(my_files, folder);
    size_t tmp_avail = 0;
    if (!tmpdir.empty()) {
        tmp_avail = espacio_tmp();
        usar_tmp = (tmp_avail > 2 * bytes_needed);
    }

    if (usar_tmp) {
        if (rank == 0) std::cout << "[INFO] Rank " << rank << " copiará archivos y procesará en /tmp (disponible: "
                                 << (tmp_avail / 1024 / 1024) << " MB, necesita: "
                                 << (bytes_needed / 1024 / 1024) << " MB)\n";
        for (const auto& f : my_files) {
            std::string src = folder + "/" + f;
            std::string dst = tmpdir + "/" + f;
            std::ifstream in(src, std::ios::binary);
            if (!in) {
                std::cerr << "[ERROR] Rank " << rank << " no pudo abrir archivo: " << src << "\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            std::ofstream out(dst, std::ios::binary);
            if (!out) {
                std::cerr << "[ERROR] Rank " << rank << " no pudo crear archivo en /tmp: " << dst << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 2);
            }
            out << in.rdbuf();
            in.close(); out.close();
        }
    } else if (!tmpdir.empty()) {
        if (rank == 0)
            std::cerr << "[WARN] Rank " << rank << " espacio en /tmp insuficiente ("
                      << (tmp_avail / 1024 / 1024) << " MB disponibles, "
                      << (bytes_needed / 1024 / 1024) << " MB requeridos)\n";
    }

    std::string chains_fname, results_fname;
    if (usar_tmp) {
        chains_fname = tmpdir + "/cadenas_rank" + std::to_string(rank) + ".ids";
        results_fname = tmpdir + "/resultados_rank" + std::to_string(rank) + ".csv";
    } else {
        chains_fname = kDirChains + "/cadenas_rank" + std::to_string(rank) + ".ids";
        results_fname = kDirResults + "/resultados_rank" + std::to_string(rank) + ".csv";
    }
    std::ofstream out_chains(chains_fname.c_str());
    if (!out_chains) {
        std::cerr << "[ERROR] Rank " << rank << " no pudo abrir archivo de salida: " << chains_fname << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 3);
    }
    std::ofstream out_results(results_fname.c_str());
    if (!out_results) {
        std::cerr << "[ERROR] Rank " << rank << " no pudo abrir archivo de salida: " << results_fname << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 4);
    }
    out_results << "timestep,cadena_id,size,f_promedio,x_prom,y_prom,z_prom\n";

    double t0_all = MPI_Wtime();
    double tot_read = 0, tot_comp = 0, tot_io = 0, tot_comm = 0;
    int tareas_procesadas = 0;

    int ts_prev = -1;
    std::vector<std::vector<int>> prev_comps;
    int ts_first = -1, ts_last = -1;
    int total_cadenas_local = 0;
    int max_cadena_local = 0;
    int min_cadena_local = std::numeric_limits<int>::max();
    long long suma_tamanios_local = 0;
    long long cantidad_cadenas_local = 0;

    for (size_t idx = 0; idx < my_files.size(); ++idx) {
        int ts = extraer_timestep(my_files[idx]);
        if (idx == 0) ts_first = ts;
        if (idx == my_files.size() - 1) ts_last = ts;

        double t;
        std::string dat_path;
        if (usar_tmp) dat_path = tmpdir + "/" + my_files[idx];
        else          dat_path = folder + "/" + my_files[idx];
        t = MPI_Wtime(); auto contacts = leer_contactos(dat_path); tot_read += MPI_Wtime() - t;
        t = MPI_Wtime(); auto comps = detectar_componentes(contacts, rank, size);     tot_comp += MPI_Wtime() - t;

        t = MPI_Wtime();
        write_chain_block(out_chains, ts, comps);
        write_resultado_block(out_results, ts, comps, contacts);
        save_tracking_file(tracking_fname, ts_prev, ts, prev_comps, comps);
        tot_io += MPI_Wtime() - t;
        tareas_procesadas++;
        ts_prev = ts;
        prev_comps = comps;

        total_cadenas_local += comps.size();
        cantidad_cadenas_local += comps.size();
        for (const auto& c : comps) {
            int sz = c.size();
            suma_tamanios_local += sz;
            if (sz > max_cadena_local) max_cadena_local = sz;
            if (sz < min_cadena_local) min_cadena_local = sz;
        }
    }
    if (cantidad_cadenas_local == 0) min_cadena_local = 0;

    out_chains.close();
    out_results.close();

    if (size > 1) {
        std::vector<std::vector<int>> first_comps, last_comps;
        if (my_count > 0) {
            first_comps = serializar_componentes(chains_fname, ts_first);
            last_comps  = serializar_componentes(chains_fname, ts_last);
        }
        if (rank < size - 1 && my_count > 0) {
            int nc = last_comps.size();
            std::vector<int> sizes(nc), flat;
            for (int i = 0; i < nc; ++i) {
                sizes[i] = last_comps[i].size();
                flat.insert(flat.end(), last_comps[i].begin(), last_comps[i].end());
            }
            if (MPI_Send(&ts_last, 1, MPI_INT, rank + 1, 100, MPI_COMM_WORLD) != MPI_SUCCESS ||
                MPI_Send(&nc, 1, MPI_INT, rank + 1, 101, MPI_COMM_WORLD) != MPI_SUCCESS ||
                MPI_Send(sizes.data(), nc, MPI_INT, rank + 1, 102, MPI_COMM_WORLD) != MPI_SUCCESS ||
                MPI_Send(flat.data(), flat.size(), MPI_INT, rank + 1, 103, MPI_COMM_WORLD) != MPI_SUCCESS) {
                std::cerr << "[ERROR] Rank " << rank << " falló enviando datos de frontera.\n";
                MPI_Abort(MPI_COMM_WORLD, 5);
            }
        }
        if (rank > 0 && my_count > 0) {
            int ts_ant, nc_ant;
            if (MPI_Recv(&ts_ant, 1, MPI_INT, rank - 1, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS ||
                MPI_Recv(&nc_ant, 1, MPI_INT, rank - 1, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
                std::cerr << "[ERROR] Rank " << rank << " falló recibiendo datos de frontera.\n";
                MPI_Abort(MPI_COMM_WORLD, 6);
            }
            std::vector<int> sizes_ant(nc_ant);
            if (MPI_Recv(sizes_ant.data(), nc_ant, MPI_INT, rank - 1, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
                std::cerr << "[ERROR] Rank " << rank << " falló recibiendo tamaños frontera.\n";
                MPI_Abort(MPI_COMM_WORLD, 7);
            }
            int total_ant = std::accumulate(sizes_ant.begin(), sizes_ant.end(), 0);
            std::vector<int> flat_ant(total_ant);
            if (MPI_Recv(flat_ant.data(), total_ant, MPI_INT, rank - 1, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
                std::cerr << "[ERROR] Rank " << rank << " falló recibiendo datos frontera.\n";
                MPI_Abort(MPI_COMM_WORLD, 8);
            }
            auto comps_ant = rebuild_components(sizes_ant, flat_ant);
            save_tracking_file(tracking_fname, ts_ant, ts_first, comps_ant, first_comps);
        }
    }

    if (usar_tmp) {
        if (std::system(("cp " + chains_fname + " " + kDirChains + "/").c_str()) != 0 ||
            std::system(("cp " + results_fname + " " + kDirResults + "/").c_str()) != 0 ||
            std::system(("cp " + tracking_fname + " " + kDirBase + "/").c_str()) != 0) {
            std::cerr << "[ERROR] Rank " << rank << " no pudo copiar archivos de resultados desde /tmp.\n";
            MPI_Abort(MPI_COMM_WORLD, 9);
        }
    }

    double sum_read = 0, sum_comp = 0, sum_io = 0, sum_comm = 0;
    int sum_tasks = 0;
    MPI_Reduce(&tot_read, &sum_read, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tot_comp, &sum_comp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tot_io, &sum_io, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tot_comm, &sum_comm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tareas_procesadas, &sum_tasks, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::ofstream out(kDirTimings + "/metrics_np" + std::to_string(size) + ".txt");
        if (!out) {
            std::cerr << "[ERROR] No se pudo escribir métricas globales.\n";
            MPI_Abort(MPI_COMM_WORLD, 11);
        }
        out << "read," << sum_read << "\n";
        out << "compute," << sum_comp << "\n";
        out << "io," << sum_io << "\n";
        out << "comm," << sum_comm << "\n";
        out << "tasks," << sum_tasks << "\n";
        out << "total_all," << (MPI_Wtime() - t0_all) << "\n";
        out.close();
    }

    int total_cadenas_global = 0;
    MPI_Reduce(&total_cadenas_local, &total_cadenas_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    int max_cadena_global = 0;
    int min_cadena_global = 0;
    MPI_Reduce(&max_cadena_local, &max_cadena_global, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&min_cadena_local, &min_cadena_global, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

    long long suma_tamanios_global = 0;
    long long cantidad_cadenas_global = 0;
    MPI_Reduce(&suma_tamanios_local, &suma_tamanios_global, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&cantidad_cadenas_local, &cantidad_cadenas_global, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    int cantidad_por_rank = total_cadenas_local;
    std::vector<int> all_cant(size, 0);
    MPI_Gather(&cantidad_por_rank, 1, MPI_INT, all_cant.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    int total_cadenas_all = 0;
    MPI_Allreduce(&total_cadenas_local, &total_cadenas_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        std::ofstream resumen((kDirTimings + "/resumen_global.txt").c_str());
        if (!resumen) {
            std::cerr << "[ERROR] No se pudo escribir resumen global.\n";
            MPI_Abort(MPI_COMM_WORLD, 12);
        }
        resumen << "Cantidad total de cadenas: " << total_cadenas_global << "\n";
        resumen << "Tamaño máximo de cadena: " << max_cadena_global << "\n";
        resumen << "Tamaño mínimo de cadena: " << min_cadena_global << "\n";
        if (cantidad_cadenas_global > 0)
            resumen << "Promedio de tamaño de cadena: "
                    << (double)suma_tamanios_global / cantidad_cadenas_global << "\n";
        resumen << "Distribución de cadenas por rank:\n";
        for (int i = 0; i < size; ++i)
            resumen << "Rank " << i << ": " << all_cant[i] << "\n";
        resumen.close();
        std::cout << "[INFO] Resumen global guardado en: " << kDirTimings << "/resumen_global.txt\n";
    }

    if (rank == 0) {
        std::vector<std::string> tracking_files;
        for (int i = 0; i < size; ++i) {
            tracking_files.push_back(kDirBase + "/tracking_rank" + std::to_string(i) + ".txt");
        }
        std::vector<std::string> all_lines;
        std::string header;
        for (const auto& fname : tracking_files) {
            std::ifstream in(fname);
            if (!in) {
                std::cerr << "[ERROR] No se pudo abrir archivo de tracking para concatenar: " << fname << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 13);
            }
            std::string line;
            bool first = true;
            while (std::getline(in, line)) {
                if (first && header.empty() && !line.empty()) {
                    header = line;
                } else if (!line.empty() && line.find("timestep_prev") == std::string::npos) {
                    all_lines.push_back(line);
                }
                first = false;
            }
            in.close();
        }
        std::sort(all_lines.begin(), all_lines.end(), [](const std::string& a, const std::string& b) {
            std::istringstream sa(a), sb(b);
            std::string token;
            int ta = -1, tb = -1;
            for (int i = 0; i < 2; ++i) std::getline(sa, token, ',');
            std::getline(sa, token, ','); ta = std::stoi(token);
            for (int i = 0; i < 2; ++i) std::getline(sb, token, ',');
            std::getline(sb, token, ','); tb = std::stoi(token);
            return ta < tb;
        });
        std::ofstream out(kDirBase + "/tracking_global.txt");
        if (!out) {
            std::cerr << "[ERROR] No se pudo crear archivo de tracking global.\n";
            MPI_Abort(MPI_COMM_WORLD, 14);
        }
        out << header << "\n";
        for (const auto& l : all_lines) out << l << "\n";
        out.close();
        std::cout << "[INFO] Archivo de tracking global generado: " << kDirBase << "/tracking_global.txt\n";
    }

    MPI_Finalize();
    return 0;
}
