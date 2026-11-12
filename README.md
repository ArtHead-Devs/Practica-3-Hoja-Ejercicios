# Práctica 2: Hoja de Ejercicios 2

Resolución de la Hoja de Ejercicios 2 de la asignatura **Bioinformática en laboratorio** (alineamientos de secuencias). Los ejercicios se resuelven en un único cuaderno de Jupyter con [Biopython](https://biopython.org/), y el ejercicio 2 trabaja con proteínas reales descargadas del NCBI.

## Autores

- Pablo Cabeza Lantigua
- Fabio Nesta Arteaga Cabrera

## Contenido de los ejercicios

| Ejercicio | Tema | Qué hace |
|---|---|---|
| 1 | Grid Search de penalizaciones de hueco | Recorre con dos bucles anidados 1000 combinaciones de `open_gap` y `extend_gap`, alinea dos secuencias de ADN con `PairwiseAligner` y guarda en un DataFrame el número de alineamientos óptimos empatados. Dibuja un mapa de calor, demuestra matemáticamente dónde se pasa de empates a una única solución y muestra los alineamientos junto a esa frontera. |
| 2 | Similitud entre especies | Descarga el citocromo c de cinco especies (humano, chimpancé, ratón, pez cebra y mosca), las lee con `Bio.SeqIO` y calcula el score de alineamiento global de todas contra todas con `BLOSUM62` (25 puntuaciones). Muestra el resultado en un mapa de calor, lo normaliza y repite el experimento con `BLOSUM45`, `BLOSUM80`, `PAM30` y `PAM250`. |
| 3 | Escalabilidad del alineador | Genera secuencias de ADN aleatorias de 100 a 5000 bases con `random`, mide con `time.perf_counter` cuánto tarda `PairwiseAligner` en alinear cada par, grafica longitud frente a tiempo y estima el exponente de crecimiento y el tiempo para 1 millón de bases. |

## Estructura del repositorio

```
.
├── datos/
│   └── citocromo_c.fasta                  Cinco secuencias de citocromo c (ejercicio 2)
├── .gitignore                             
├── .python-version                        Versión de Python
├── Ejercicios BIO-GCID Práctica 2.pdf     Enunciado de la práctica
├── LICENSE
├── README.md                              
├── ejercicios.ipynb                       Cuaderno con la resolución de todos los ejercicios
├── pyproject.toml                         Dependencias del proyecto
└── uv.lock                                Versiones fijadas de las dependencias

```

La descarga del ejercicio 2 solo se hace si `datos/citocromo_c.fasta` no existe, así que se puede trabajar sin conexión una vez descargado.

## Requisitos

- Python 3.12 o superior.
- [uv](https://docs.astral.sh/uv/) como gestor de entornos y dependencias.
- Conexión a internet solo para la primera descarga de las proteínas.
- Dependencias, definidas en `pyproject.toml`: `biopython`, `pandas`, `matplotlib`, `seaborn` y `python-dotenv`. Para el cuaderno también se usa `ipykernel`, incluido en el grupo de desarrollo.

## Instalación y ejecución

1. Clona el repositorio y entra en la carpeta del proyecto.
2. Instala las dependencias:

   ```bash
   uv sync
   ```

3. Crea un fichero `.env` en la raíz del proyecto con tu correo. El NCBI exige identificarse en cada consulta de Entrez, y el cuaderno lee la variable `ENTREZ_EMAIL` una sola vez al principio. El fichero está en `.gitignore`, así que no se sube al repositorio.

   ```
   ENTREZ_EMAIL=tu_correo@ejemplo.com
   ```

4. Abre `ejercicios.ipynb` en VS Code o en Jupyter, selecciona el kernel del entorno `.venv` y ejecuta las celdas en orden.

## Decisiones y notas por ejercicio

### Ejercicio 1

- **Intervalo de búsqueda.** El enunciado da los límites y deja elegir el paso: `open_gap` de -1 a -20 con paso 1 (20 valores) y `extend_gap` de -0.1 a -5.0 con paso 0.1 (50 valores), en total 1000 combinaciones.
- **Puntuaciones de match y mismatch.** El enunciado no las fija. Se usan `match = 5` y `mismatch = -4`, los valores de los ejemplos del profesor, en modo global. Biopython usa por defecto `1` y `0`.
- **Punto de colapso.** Las 8 bases de la secuencia del paciente aparecen en orden en la sana, así que solo hay que decidir cómo se reparten 6 huecos en bloques. El score es `40 + 6·extend + k·(open − extend)`, donde `k` es el número de bloques. Hay una única solución cuando `|extend| < |open|` y varios empates cuando `|extend| ≥ |open|`. La frontera es `open = extend`, y la fórmula se comprueba contra Biopython en las 1000 combinaciones.
- **Dependencia de match y mismatch.** El apartado 8 repite la búsqueda con `5/-4`, `1/0` y `2/-1`. El punto de colapso es el mismo en los tres. La fórmula del score solo es exacta con `5/-4`: con `1/0` un fallo apenas cuesta y en parte de la zona aparecen otros alineamientos.

### Ejercicio 2

- **Proteína elegida.** Citocromo c, muy conservada, descargada por su código de UniProt con `Entrez.efetch` sobre la base de datos `protein`:

  | Especie | Código |
  |---|---|
  | Humano | `P99999` |
  | Chimpancé | `P99998` |
  | Ratón | `P62897` |
  | Pez cebra | `Q6IQM2` |
  | Mosca | `P84029` |

- **Penalizaciones de hueco.** El enunciado no las fija. Se usan `open_gap = -10` y `extend_gap = -0.5`, valores habituales en proteínas.
- **Normalización.** El score bruto depende de la longitud de cada secuencia, así que para comparar parejas y matrices se divide entre la media geométrica de los scores de las dos secuencias consigo mismas: `score(A, B) / √(score(A, A) · score(B, B))`.
- **Resultado.** Con las cinco matrices el orden de parecido respecto al humano es el mismo: chimpancé (idéntico), ratón, pez y mosca. Cambiar de matriz modifica la escala, no la conclusión.

### Ejercicio 3

- **Medición.** Una sola medición por tamaño, con el alineador con sus parámetros por defecto y una semilla fija en `random` para que las secuencias sean reproducibles.
- **Resultados variables.** Los tiempos cambian un poco en cada ejecución. El exponente estimado sale entre 1,8 y 1,9, es decir, un crecimiento aproximadamente cuadrático, O(n²).
- **Extrapolación.** Con ese ajuste, un solo par de secuencias de 1 millón de bases tardaría más de una hora, además de necesitar una matriz de unos 10¹² celdas. Si lo que crece es el número de pares de secuencias cortas, el coste es lineal y asumible.

## Licencia

Consulta el fichero [LICENSE](LICENSE).
