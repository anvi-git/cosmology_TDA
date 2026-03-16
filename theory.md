## 3. Mathematical Theory

### 3.1 Simplicial Complexes

A simplicial complex is built from
- **0-simplices (vertices)**: Individual points (halos)
- **1-simplices (edges)**: Connections between pairs of points
- **2-simplices (triangles)**: Connections between triples
- **3-simplices (tetrahedra)**: Connections between quadruples

**Alpha Complex**: Vertices connected if within distance r (Voronoi-based)
$$\text{Alpha complex} = \{\text{simplices from Delaunay triangulation within radius }r\}$$

### 3.2 Filtration Methods

#### Method 0: Alpha Filtration
Connect points based on **Euclidean distance**
$$\alpha_r = \{(i,j) : |x_i - x_j| < r\}$$

#### Method 1: Sublevel Filtration
Based on field values (cubical complex)
$$\text{Filtration based on height function values}$$

#### Method 2: Alpha-DTM Filtration
Combines alpha complex with **Distance-To-Measure (DTM)**:
$$\text{DTM}_m(x) = \sqrt{\frac{1}{k} \sum_{i=1}^{k} d_i(x)^2}$$
where $d_i(x)$ is distance to i-th nearest neighbor, $k = \lfloor m \cdot N \rfloor$

#### Method 3: Alpha-DTM with Length-DTM Mixing
Enhanced version combining:
- **Density information** (DTM at query points)
- **Edge weighting** (distance-dependent filtration values)

Edge filtration value:
$$\text{filt}(e_{ij}) = \begin{cases}
\sqrt{\frac{((f_i + f_j)^2 + d^2) \cdot ((f_i - f_j)^2 + d^2)}{4d^2}} & \text{if } p=2 \\
\frac{f_i + f_j + d}{2} & \text{if } p=1 \\
\max(f_i, f_j, d/2) & \text{if } p=\infty
\end{cases}$$
where $f_i, f_j$ are DTM values and $d$ is edge length.

### 3.3 Homology & Persistence

**Persistent Homology** tracks topology across scales:

$$H_p(\text{Complex}_\alpha) = \text{p-th homology group at scale }\alpha$$

For each feature (cycle):
- **Birth** ($b$): Scale where feature appears
- **Death** ($d$): Scale where feature merges/disappears
- **Persistence** ($p = d - b$): Lifetime of feature

**Interpretation**:
- **High persistence**: Robust topological feature (signal)
- **Low persistence**: Noise

### 3.4 Persistence Diagrams

Graphically represent persistence as (birth, death) points:
- Points far from diagonal: strong, persistent features
- Points near diagonal: weak, short-lived features

Three types by dimension:
- **H₀** (0-dimensional): Connected components/clusters
- **H₁** (1-dimensional): Holes/loops (void structure)
- **H₂** (2-dimensional): Cavities (enclosed regions)