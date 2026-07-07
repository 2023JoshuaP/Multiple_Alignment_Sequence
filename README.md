# Multiple Sequence Alignment (MSA) - TIF

Este repositorio contiene la implementación desde cero de algoritmos clásicos de Alineamiento Múltiple de Secuencias (MSA), desarrollados como parte de un Trabajo Integral Final (TIF).

El proyecto está diseñado de forma modular, separando las matemáticas puras de la orquestación de los algoritmos de alto nivel.

## Estructura del Proyecto

* **`core/`**: El "corazón" matemático del proyecto. Contiene las implementaciones base:
  * `pairwise_nw.py`: Algoritmo de Needleman-Wunsch para alineamiento global de pares.
  * `distance_matrix.py`: Cálculo de matrices de distancias a partir de alineamientos.
  * `upgma.py`: Construcción de árboles filogenéticos mediante UPGMA.
  * `neighbor_joining.py`: Construcción de árboles filogenéticos mediante Neighbor-Joining.
  * `profile_alignment.py`: Alineamiento progresivo de perfiles.
* **`algorithms/`**: Orquestadores de alto nivel.
  * `clustalw.py`: Implementación del algoritmo Progresivo ClustalW.
* **`papers/`**: Literatura científica y artículos base de donde se han extraído los algoritmos (ClustalW, MUSCLE, Recocido Simulado, etc.).

## Ejecución

Para ejecutar el algoritmo ClustalW con secuencias de prueba:

```bash
python3 algorithms/clustalw.py
```

## Próximos Pasos (En Desarrollo)

* [x] Alineamiento Progresivo (ClustalW)
* [ ] Alineamiento Iterativo (MUSCLE)
* [ ] Alineamiento Heurístico (Simulated Annealing)
