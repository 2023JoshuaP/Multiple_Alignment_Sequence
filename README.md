# Multiple Sequence Alignment (MSA) - TIF

Este repositorio contiene la implementación desde cero de algoritmos clásicos de Alineamiento Múltiple de Secuencias (MSA), desarrollados como parte de un Trabajo Integral Final (TIF).

El proyecto está diseñado de forma puramente modular, separando las matemáticas subyacentes de la orquestación de los algoritmos de alto nivel, y estructurado siguiendo estrictamente las literaturas científicas originales (Robert Edgar, 2004 para MUSCLE, etc).

## Estructura del Proyecto

* **`core/`**: El "corazón" matemático del proyecto. Contiene las herramientas lógicas base:
  * `pairwise_nw.py`: Algoritmo Needleman-Wunsch para alineamiento global.
  * `distance_matrix.py`: Matrices de distancias tradicionales a partir de alineamientos.
  * `kmer_distance.py`: Conteo fraccional rápido por palabras (K-mers) para aceleración.
  * `kimura_distance.py`: Correcciones evolutivas de divergencia.
  * `upgma.py` / `neighbor_joining.py`: Construcción matemática de árboles filogenéticos.
  * `sp_score.py`: Función objetivo matemática (Sum-of-Pairs) para calificar iteraciones.
  * `profile_alignment.py`: Alineamiento progresivo entre bloques/perfiles.
  * `tree_plotter.py`: Generador de gráficos de cladogramas visuales usando Matplotlib.
  * `fasta_parser.py`: Lector y codificador de archivos FASTA estandarizados.
  * `sa_moves.py`: Operadores de mutación estocástica (gap shifting, etc.) para exploración heurística.

* **`algorithms/`**: Orquestadores de alto nivel.
  * `clustalw.py`: Implementación del algoritmo **Progresivo** (ClustalW) guiado por NJ.
  * `muscle.py`: Implementación del algoritmo **Iterativo** (MUSCLE) con sus 3 Fases completas de refinamiento.
  * `simulated_annealing.py`: Implementación del algoritmo **Heurístico** basado en recocido termodinámico estocástico para evitar trampas locales.

* **`benchmarks/`**: Carpeta autogenerada donde se exportan los árboles visuales (`.png`) tras la ejecución.

* **`papers/`**: Literatura científica y artículos base.

## Ejecución e Interfaz

El proyecto cuenta con una Interfaz Gráfica (GUI) construida en `Tkinter`, que permite cargar archivos de secuencias `.fasta`, seleccionar el algoritmo a utilizar y visualizar los resultados del análisis en tiempo real. 

Para ejecutar la herramienta:

```bash
python3 main.py
```

*Nota: Requiere tener instalada la librería gráfica `matplotlib` (`pip install matplotlib`).*

## Estado del Proyecto

* [x] Alineamiento Progresivo (ClustalW)
* [x] Alineamiento Iterativo (MUSCLE)
* [x] Alineamiento Heurístico (Simulated Annealing)