# Persistent Homology for Large Scale Structure (PHLSS)
## Comprehensive Guide: Theory, Implementation, and Applications

**Authors:** Matteo Biagetti and Alex Cole  
**Reference:** [The Persistence of Large Scale Structure I: Primordial non-Gaussianity](https://arxiv.org/abs/2009.04819)  
**Public Release:** v1.0

---
## Table of Contents
1. Introduction to Persistent Homology
2. Large Scale Structure & Cosmology
3. Mathematical Theory
4. Implementation Architecture
5. Workflow & Usage
6. Code Components Analysis
7. Practical Examples
---
## 1. Introduction to Persistent Homology

### What is Persistent Homology?

Persistent Homology is a method from **Topological Data Analysis (TDA)** that studies the **topological features** of point clouds across multiple scales.

Key concepts:
- **Homology**: Studies hole-like structures in data (0-dim: components/clusters, 1-dim: loops/voids, 2-dim: cavities)
- **Persistence**: Tracks how long these features exist as scale parameter changes
- **Filtration**: Parametrized family of nested complexes at increasing scales

### Why It Matters for Cosmology

The large-scale structure of the universe exhibits:
- **Clusters** of dark matter halos (filaments, clusters, nodes)
- **Voids** between these structures
- **Topological structures** encoded in cosmological parameters

Persistent homology captures these features to:
1. Constrain cosmological parameters (e.g., primordial non-Gaussianity)
2. Extract scale-dependent information
3. Avoid relying only on two-point statistics (power spectrum)
---
## 2. Large Scale Structure & Cosmology

### Background

- **N-body simulations**: Evolve billions of dark matter particles from early universe to present day
- **Halo catalogs**: Extract self-bound structures (halos) representing galaxy clusters
- **Positions & Masses**: Halos have 3D spatial coordinates and mass estimates

### Box Geometry

The code operates on periodic cubic boxes:
- **Main box**: 2000 Mpc/h (comoving megaparsec/h) - full simulation volume
- **Sub-boxes**: 1000 Mpc/h - subdomain partitions for computational efficiency
- **Redshift**: Analysis at redshift z=1 (intermediate cosmic epoch)

### Halo Selection

```
Mass range: 9.189×10^12 to 10^20 Msun/h
Total halos per subbox: 252,321
Subsampling: For computational tractability
```
---
## 3. Mathematical Theory

### 3.1 Simplicial Complexes

A simplicial complex is built from:
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
---
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
---
### 4.2 Key Classes

#### `DIAG` (pdiagram.py)
**Purpose**: Compute persistence diagrams from halo catalog

**Main methods**:
- `__init__`: Initialize with data, split into sub-boxes
- `parallelPD()`: Parallel computation across sub-boxes
- `computeSavePDalpha()`: Alpha filtration method
- `computeSavePDalphaDTM()`: Alpha-DTM filtration method
- `AlphaWeightedFiltration()`: DTM-weighted edge values
- `DTM()`: Distance-to-measure computation
- `Filtration_value()`: Complex edge weight calculation

**Output**: Persistence diagrams saved as CSV files with format:
```
dimension, birth, death
0, 0.5, 1.2
1, 0.3, 0.8
...
```

#### `IMAG` (pimages.py)
**Purpose**: Convert persistence diagrams to persistence images (2D images)

**Main methods**:
- `parallelPI()`: Parallel computation across sub-boxes
- `PI()`: Kernel density estimate of persistence diagram
- `V()`: 1D persistence image (death axis)
- `H()`: 1D persistence image (birth axis)
- `VH()`: 1D persistence image (persistence axis: death - birth)
- `Betti()`: Betti curve computation
- `calcSavePI()`: Calculate and save persistence image

**Output**: 2D arrays (images) that can be used as ML features
---
## 5. Workflow & Usage

### 5.1 Parameter Configuration

Edit `param.py` to configure:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cat` | 1 | Catalog type (1=Rockstar, 2=Minerva, 3=Flagship) |
| `box` | 2000 | Box size in Mpc/h |
| `smallbox` | 1000 | Sub-box size in Mpc/h |
| `z` | 1 | Redshift |
| `mmin` | 9.189e12 | Minimum halo mass in M☉/h |
| `mmax` | 1e20 | Maximum halo mass in M☉/h |
| `method` | 3 | Filtration method (0-3) |
| `nhalo` | 252321 | Halos per sub-box |
| `thresh` | 1.3 | Outlier threshold |

### 5.2 Command-Line Arguments

```bash
python main.py \
  --hfile <halo_catalog_path> \
  --pd True \
  --pdfile <output_pd_path> \
  --pi True \
  --npool <num_processors> \
  --cut <persistence_cut> \
  --boundsName <bounds_file>
```

**Arguments**:
- `--hfile`: Path to halo catalog (required)
- `--pd`: Compute persistence diagrams (True/False)
- `--pdfile`: Output/input persistence diagram file (required)
- `--pi`: Compute persistence images (True/False)
- `--npool`: Number of parallel processors (default: 2)
- `--cut`: Minimum persistence threshold (default: 0)
- `--boundsName`: Bounds file for statistics (default: '0')
- `--arrayjob`: Enable SLURM array job mode (for HPC)
- `--onedim`: Compute 1D functions (advanced)
---
### 5.3 Example Execution

#### Basic Example: Test Dataset

```bash
cd Code
python main.py \
  --hfile ../TestData/out_z2_nosub_M200b_20p.list \
  --pd True \
  --pdfile ../Output/Subsampling/PD/test_alphaDTM_ \
  --pi True \
  --npool 4
```

#### Full-Scale Simulation

```bash
python main.py \
  --hfile /path/to/full/halo/catalog.list \
  --pd True \
  --pdfile Output/PD/full_sim_run1_ \
  --pi True \
  --npool 16 \
  --boundsName Output/bounds_run1.txt
```

#### HPC Submission (SLURM)

```bash
sbatch Scripts/python_script_cartesius.sh
```
---
## 6. Code Components Analysis

### 6.1 Main Entry Point (`main.py`)

**Flow**:
1. Parse command-line arguments via `user.parse_arguments()`
2. Load halo catalog with `readcat.readcats()`
3. Filter halos by mass
4. If `--arrayjob`, modify filenames for SLURM array tasks
5. If `--pd True`, compute persistence diagrams via `DIAG` class
6. Compute bounds file for statistical analysis
7. If `--pi True`, compute persistence images via `IMAG` class

### 6.2 Catalog Reading (`DataReading/readcatalogs.py`)

**Supported formats**:
- **Type 1 (Rockstar)**: Plain text, columns = [Mass, X, Y, Z]
- **Type 2 (Minerva)**: Binary format with structured array
- **Type 3 (Flagship)**: Parquet format with log-mass extraction

**Output**: Numpy array shape (N_halos, 4) with [Mass, X, Y, Z]

### 6.3 Persistence Diagram Computation

**Key process**:
1. Split box into 8×8×8 = 512 sub-boxes (for 2000 Mpc/h with 1000 Mpc/h sub-boxes)
2. For each sub-box:
   - Extract halos within boundaries
   - Subsample to `nhalo` = 252,321 halos per sub-box
   - Build simplicial complex using GUDHI
   - Compute persistence (algorithm: Zomorodian-Carlsson)
   - Save results as CSV
3. Parallelize across processors using `parmap`
---
### 6.4 Distance-To-Measure (DTM) Implementation

```python
def DTM(X, query_pts, m):
    '''Compute distance-to-measure function'''
    N_tot = X.shape[0]
    k = int(np.floor(m * N_tot)) + 1
    
    # Build KD-tree for efficient nearest neighbor search
    kdt = KDTree(X, leaf_size=30, metric='euclidean')
    
    # Query k-nearest neighbors
    NN_Dist, NN = kdt.query(query_pts, k, return_distance=True)
    
    # Compute DTM
    DTM_result = np.sqrt(np.sum(NN_Dist * NN_Dist, axis=1) / (k-1))
    
    return DTM_result
```

**Intuition**: DTM estimates local density. Dense regions = small DTM, sparse regions = large DTM.

### 6.5 Persistence Image Computation

**Persistence Image (PI)**: Convert persistence diagram to 2D image

$$\text{PI}(x,y) = \int \rho(b, p) \cdot K(x,y,b,p) \, db \, dp$$

where:
- $\rho(b, p)$: Point weights (persistence values)
- $K$: Kernel (epanechnikov, gaussian, etc.)
- $x,y$: Grid coordinates

**Implementation in code**:
- Use kernel density estimation (KDE) with bandwidth $\sigma$
- Grid resolution: 30×30 default
- Weight by persistence value
---
### 6.6 Auxiliary Functions (`fclaux.py`)

**`loadfile(fname)`**
- Parse text-based halo catalog
- Extract mass and 3D position

**`savePD(diag, filename)`**
- Convert GUDHI persistence output to CSV
- Format: [dimension, birth, death]
- Skip features with birth = death (trivial features)

**`cleanPD(pd, cut, power)`**
- Transform coordinates: (birth, death) → (birth, persistence)
- Apply power-law rescaling: $\text{value} \to \text{value}^p$
- Filter by persistence threshold
- Separate by dimension (H₀, H₁, H₂)

**Purpose of rescaling**: Emphasize different scales depending on cosmic analysis
---
### 8.2 Cosmological Interpretation

**Why persistence homology is valuable for cosmology**:

1. **Primordial non-Gaussianity (PNG)**
   - Traditional power spectrum only captures Gaussian statistics
   - PNG leaves imprints on topological features
   - Persistence homology captures these signatures

2. **Scale-dependent information**
   - Birth scale: When structure forms
   - Death scale: When it merges
   - Different scales → different physics

3. **Complementary to other methods**
   - Power spectrum: 2-point statistics
   - Bispectrum: 3-point statistics  
   - Persistent homology: **Topology** - independent information

4. **Mass selection effects**
   - Mass threshold changes halo population
   - Affects topological features
   - Can constrain halo mass function
---
## Summary & Best Practices

### Key Takeaways

1. **PHLSS** applies topological data analysis to cosmology
2. **Workflow**: Halo catalog → Persistence diagrams → Persistence images → ML
3. **Four methods** with increasing noise robustness (Method 3 best)
4. **Parallel processing** essential for large catalogs (512 sub-boxes)
5. **Parameter tuning** affects physics extraction

### Tuning Recommendations

```python
# param.py settings
method = 3                    # Use most robust method
nhalo = 252321               # Good balance computational/statistics
mmin = 9.189e12              # Massive halos (more robust)
thresh = 1.3                 # Outlier removal
smallbox = 1000              # Balance parallel efficiency

# Command line
npool = min(32, n_cores)     # Max parallelization
cut = 0.01                   # Minimal persistence threshold
```

### Output Interpretation

- **Persistence diagrams**: Raw topological features (birth, death) by dimension
- **Persistence images**: Machine-learning-ready 2D feature images
- **Betti curves**: 1D summaries of topological evolution
- **Statistics files**: Bounds and outlier information for further analysis
---
## References & Further Reading

### Primary Literature
- Biagetti, M. & Cole, A. (2020). *The Persistence of Large Scale Structure I: Primordial non-Gaussianity*. [arXiv:2009.04819](https://arxiv.org/abs/2009.04819)
- N-body simulations: [Cosmos simulation suite](https://mbiagetti.gitlab.io/cosmos/nbody/eos/)

### Topological Data Analysis
- Edelsbrunner, H. & Harer, J. (2010). *Computational Topology: An Introduction*. AMS
- Adams, H., et al. (2017). *Persistence Images: A Stable Vector Representation of Persistent Homology*. [arXiv:1505.04087](https://arxiv.org/abs/1505.04087)
- Chazal, F., et al. (2016). *Statistical topological data analysis using persistence landscapes*. JMLR

### Software & Libraries
- **GUDHI**: [Geometry Understanding in Higher Dimensions](http://gudhi.gforge.inria.fr/) - Simplex tree & persistence computation
- **scikit-learn**: KD-trees, KDE, machine learning
- **parmap**: Parallel mapping for easy parallelization
- **NumPy/SciPy**: Scientific computing foundations