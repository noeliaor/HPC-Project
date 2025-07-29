Detección y análisis de cadenas de fuerza en simulaciones granulares 3D

Este repositorio contiene el código de un proyecto de análisis en paralelo de cadenas de fuerza en simulaciones de medios granulares 3D, ejecutable en un entorno de cómputo de alto rendimiento (HPC) utilizando MPI.
El objetivo principal es detectar cadenas de fuerza a partir de archivos de contactos generados por simulaciones realizadas ESyS-Particle, analizar su evolución temporal y determinar métricas físicas (tamaño, fuerza, persistencia, etc.) así como métricas de rendimiento para el código (tiempos por fase, escalabilidad, etc.).

Estructura del repositorio
main_mpi.cpp: Programa principal con paralelización mediante MPI.

contact.{cpp,hpp}: Definición de la estructura Contact y funciones auxiliares.

graph.{cpp,hpp}: Construcción del grafo de contactos y detección de componentes conexas.

io_utils.{cpp,hpp}: Funciones de entrada/salida para leer archivos .dat y escribir resultados.

utils.{cpp,hpp}: Utilidades generales, como cálculo de centroides de la interacción, similitud de Jaccard y carga de contactos.

Makefile: Archivo para compilar el proyecto con make.

Requisitos
Compilador C++11 (por ejemplo, g++)
Implementación de MPI (por ejemplo, OpenMPI o MPICH)
Sistema operativo Linux (en este caso se utilizó el entorno HPC con SLURM de ClusterUY)

Compilación
Para compilar el proyecto, usar:
make
Esto generará el ejecutable detectar_cadenas_mpi.

Ejecución
El programa se ejecuta con mpirun o srun, dependiendo del entorno. Por ejemplo:

mpirun -np 8 ./detectar_cadenas_mp ruta_a_archivos_partForce
Se recomienda correrlo en un clúster y asignar un nodo por proceso para maximizar la eficiencia de E/S.

Entrada
Archivos de contacto en formato partForce.*.dat generados por ESyS-Particle.
Cada archivo contiene la información de los contactos entre partículas en un timestep específico.

Salida
Archivos .csv con estadísticas por cadena.
Archivos .ids con los identificadores de partículas que conforman cada cadena.
Archivos de tracking para seguir la evolución de cadenas entre timesteps.
Métricas de rendimiento por fase (lectura, cómputo, escritura, comunicación).

Tracking de cadenas
El programa incluye unaetapa de tracking donde se asocian cadenas de fuerza entre snapshots sucesivos utilizando el índice de Jaccard. Esto permite analizar la persistencia temporal de las cadenas y su evolución estructural.
