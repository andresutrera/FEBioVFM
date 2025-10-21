# FEBio VFM

## Architecture Overview

```
src/
├── task/             # FEBio FECoreTask implementations (VFMTask)
├── optimization/     # VFMProblem builder, LM solver, internal/external work assemblers
├── state/            # Strongly-typed storage for measured/virtual histories and tensors
├── domain/           # Generic nodal/element field containers with ragged support
├── services/         # Kinematic reconstruction and stress evaluation utilities
├── FE/               # FEBio integration helpers (shape, material provider, parameter applier)
├── io/               # XML reader and XPLT exporter
├── vfm/              # Algebraic assemblers for internal and external work
├── levmar/           # Levenberg–Marquardt routine (dlevmar_dif)
└── plugin/           # DLL entry glue for FEBio
```

The code is organised around a small number of data-carrier classes (`state/`, `domain/`) and stateless services (`services/`, `optimization/`, `vfm/`) that act on them. FEBio-specific hooks live in `task/` and `FE/`, while IO concerns are isolated in `io/`.

## Task Pipeline

1. **Task initialisation (`task/vfm_task.cpp`)**
   - `VFMTask::Init` loads the optimisation XML through `VFMXmlReader`.
   - `prepare_vfm_problem` (from `optimization/vfm_problem.cpp`) builds a `VFMProblem` containing:
     - FEM mesh metadata (`MeshDims`, `MeshConn`, `MeshQuad`).
     - Experiment data in a `VFMState` instance (measured nodal histories, virtual fields, loads).
     - Pre-sized tensor containers (`Deformations`, `VirtualDeformations`, `Stresses`).
     - Precomputed virtual-kinematic fields (`services/kinematics.cpp`).
     - Initial FEBio stresses via `StressEval::cauchy` and `StressEval::first_piola`.
   - FEBio parameters are applied through `FEBioParameterApplier`; materials are reinitialised to reflect the starting guess.

2. **Optimisation (`optimization/vfm_solver.cpp`)**
   - `solve_vfm_problem` orchestrates:
     - Parameter updates (lambda closures) and stress recomputation.
     - Internal work assembly (`vfm/InternalWork.hpp`).
     - External work assembly (`vfm/ExternalVirtualWork.hpp`).
     - Non-linear least squares solve via `run_vfm_levmar` (wrapper over `dlevmar_dif`).
   - Upon convergence the physical parameters inside `VFMState` are updated and a final stress refresh ensures tensors match the optimal solution.

3. **Export (`io/exporter.cpp`)**
   - `export_vfm_results` writes an XPLT file using FEBio’s `FEBioPlotFile`:
     - Measured displacement field at nodes for every frame.
     - Measured deformation gradient per element (Gauss-point average).
     - Cauchy stress and First Piola stress per element.
     - For each virtual field: displacement and deformation gradient histories.

## Data Containers

- **Measured data (`MeasuredData`)** – ragged time series of nodal `vec3d`. Each frame is a `MeasuredFrame` containing a `NodalField<vec3d>`.
- **Virtual fields (`VirtualFields`)** – collection of time series indexed by virtual-field id; storage mirrors `MeasuredData`.
- **Loads (`MeasuredLoad`)** – time steps with named surface forces. Each entry contains `{surface, vec3d}` pairs aligned with `MeshDims`.
- **Deformations (`Deformations`)** – ragged element×Gauss storage (`RaggedElemField<mat3d>`) for measured deformation gradients per frame.
- **Virtual deformations (`VirtualDeformations`)** – same structure per virtual field.
- **Stresses (`Stresses`)** – per-frame RaggedElemField for both Cauchy (`sigma`) and First Piola (`P`) stresses.

These containers expose `setElemShape`, `addTime`, and direct `ref`/`cref` accessors so assemblers can operate without copying.

## FEBio Integration Helpers

- `FE/shape_provider_febio` supplies shape-function gradients for each Gauss point to the kinematics service.
- `FE/material_provider_febio` clones FEBio material points, injects trial deformation gradients, and evaluates stress consistently for both uncoupled and full materials.
- `FE/params_febio` caches raw parameter pointers and synchronises them with `VFMState`.

## Assembly Mechanics

### Kinematics

`services/kinematics.cpp` transforms nodal displacements into deformation gradients:

$$
F(\mathbf{x}) = \mathbf{I} + \sum_{a} \mathbf{u}_a \otimes \nabla N_a
$$

with gradients computed in the reference configuration via the shape provider. Virtual fields use the same routine but operate on prescribed virtual displacements.

### Stress Evaluation

`StressEval::cauchy` and `StressEval::first_piola` run in two phases:

1. For each Gauss point, clone FEBio’s material point, assign $F$, zero stateful variables, evaluate either $\boldsymbol{\sigma}$ or its deviatoric part (for uncoupled materials).
2. Convert to First Piola through:

$$
\mathbf{P} = J \, \boldsymbol{\sigma} \, \mathbf{F}^{-T}
$$

The tensors are stored in `Stresses` so downstream assemblers have random access.

### External Virtual Work

For each virtual field $v$ and frame $t$ the external work is

$$
W_{\text{ext}}^{(v,t)} = \sum_{k \in \text{surfaces}} \mathbf{F}_k^{(t)} \cdot \mathbf{u}_k^{* (v,t)}
$$

where $\mathbf{F}_k$ are resultant forces supplied in the XML per named surface, and $\mathbf{u}_k^*$ are the virtual displacements mapped to those surface nodes. `ExternalVirtualWorkAssembler` pre-resolves surface names to node sets using `build_surface_info` so evaluation is $O(\text{#nodes on surface})$.

### Internal Virtual Work

Internal work uses the Gauss-point form:

$$
W_{\text{int}}^{(v,t)} = \sum_{e} \sum_{g} \left( \mathbf{P}_{e,g}^{(t)} : \mathbf{G}_{e,g}^{(v,t)} \right) \, w_{e,g}
$$

with $\mathbf{G} = F^{*(v)} - \mathbf{I}$, and weights $w_{e,g} = \det J_0 \, w_g$ from `MeshQuad`. `InternalWorkAssembler` computes the double contraction using `mat3d::dotdot` and accumulates one scalar per $(v,t)$ pair.

### Non-linear Solve

The unknown parameter vector $\boldsymbol{\theta}$ (one entry per optimised FEBio material parameter) is passed to `dlevmar_dif`. The residual for virtual field $v$ and frame $t$ is

$$
r_{(v,t)}(\boldsymbol{\theta}) = W_{\text{ext}}^{(v,t)} - W_{\text{int}}^{(v,t)}(\boldsymbol{\theta})
$$

The solver minimises $\| \mathbf{r}(\boldsymbol{\theta}) \|_2^2$ with Levenberg–Marquardt trust-region steps. Parameter updates are applied directly to the FEBio model via `FEBioParameterApplier`, followed by a stress refresh so the residual reflects the current iterate.

## Exported Data Sets

The XPLT exporter registers the following plot variables:

- `measured displacement` – nodal measured displacements.
- `measured deformation gradient` – element-averaged deformation gradient.
- `cauchy stress` – element-averaged symmetric stress.
- `first piola stress` – element-averaged first Piola stress.
- `virtual displacement i` / `virtual deformation gradient i` – one pair per virtual field index $i$ (name suffixed when multiple fields exist).

Frames are emitted in time-index order; each variable is set to the correct pointer-backed storage immediately before `plot.Write`.

## Mathematical Snapshot

Combining the pieces, the discretised VFM statement solved by the code is:

$$
\min_{\boldsymbol{\theta}} \; \left\| \mathbf{W}_{\text{ext}} - \mathbf{W}_{\text{int}}(\boldsymbol{\theta}) \right\|_2^2
$$

subject to FEBio material evaluations providing $\mathbf{P}(\boldsymbol{\theta})$. The equality $W_{\text{ext}}^{(v,t)} = W_{\text{int}}^{(v,t)}$ encodes momentum balance for each admissible virtual field. The assemblers and solver ensure that Jacobian-free residual evaluations remain consistent with FEBio’s constitutive response, allowing the solver to refine parameters until the virtual work balance is satisfied over all prescribed frames and virtual fields.
