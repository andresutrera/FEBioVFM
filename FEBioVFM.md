# FEBio Virtual Fields Method Plug-in – Implementation Notes

This document summarises the current state of the FEBio VFM plug-in, the data
structures that have been implemented, and the mathematical background that
underpins each component.  It is intended as living documentation while the
full Virtual Fields Method workflow is assembled.

---

## 1. High-Level Flow

At the moment the plug-in focuses on bootstrapping a Virtual Fields analysis by:

1. Parsing FEBio optimisation input files (`<febio_optimize>`), including scalar
   model parameters and externally supplied displacement fields.
2. Validating that the underlying FE model satisfies the assumptions of the
   VFM scaffold (solid-only domains and consistent displacement data).
3. Reconstructing the deformation gradient field \(\mathbf{F}\) at every Gauss
   point using the measured nodal displacements, without mutating the FEBio
   mesh state.

The entry point is `FEVFMTask::Init`, which orchestrates the steps above after
the optimisation data container (`FEOptimizeDataVFM`) has been populated.

---

## 2. Input Parameters and Scalar Linking

Scalar optimisation variables are exposed via two helper classes:

- `FEInputParameterVFM` – abstract base describing optimisation variables,
  storing metadata such as bounds, scaling, and initial guesses.
- `FEModelParameterVFM` – concrete adapter that resolves a specified FEBio
  model parameter, caches the pointer to its underlying double value, and keeps
  the optimiser in sync with the FE model.

`FEOptimizeDataVFM` owns the collection of input parameters and exposes
high-level operations (`Input`, `Init`, `Solve`, `FESolve`).  Although the solve
loop is still a placeholder, the infrastructure for registering parameters and
running initialisation has been completed.

---

## 3. Displacement Containers

Measured and virtual displacements are represented by the header-only
`DisplacementContainer`.  Each entry is a `NodeDisplacement` struct containing
a FEBio node identifier and a triple \((u_x, u_y, u_z)\).  The container offers:

- `Add` and `Clear` helpers,
- constant-time lookups via `TryGet(nodeId, disp)`,
- iteration over the insertion-ordered sample vector.

This structure is used both for measured data (coming from, for instance,
Digital Image Correlation) and for the prescribed virtual fields required by
the VFM.

---

## 4. XML Parsing

`FEVFMInput` now recognises the following sections inside the optimisation file:

```xml
<Parameters>
    <param name="...">init, min, max, scale</param>
</Parameters>
<MeasuredDisplacements>
    <node id="...">ux, uy, uz</node>
    ...
</MeasuredDisplacements>
<VirtualDisplacements>
    <node id="...">ux, uy, uz</node>
    ...
</VirtualDisplacements>
```

Both `<node>` and legacy `<elem>` tags are accepted for backwards compatibility,
but the data are always treated as nodal quantities.  After parsing, the two
containers live on a time history; access individual steps with
`FEOptimizeDataVFM::MeasuredHistory()[i]` and `::VirtualHistory()[i]`, each of
which exposes the samples collected for that time point.

Measured surface loads follow a similar pattern:

```xml
<MeasuredLoads>
    <time t="0.0">
        <surface id="left_grip">Fx, Fy, Fz</surface>
        <surface id="right_grip">Fx, Fy, Fz</surface>
    </time>
</MeasuredLoads>
```

Loads are available through `FEOptimizeDataVFM::MeasuredLoads()[i].loads`.
Individual surfaces can be queried by name with `SurfaceLoadSet::Find` or
`::TryGet`, both of which return `vec3d` force vectors.

For debugging, a `LogDebugSummary()` method lists every parsed parameter and
displacement sample using `feLogDebugEx`, so the information is only displayed
when FEBio runs in debug mode (e.g. with the `-g` flag).

---

## 5. Model Validation

Prior to any heavy computation we check that the FEBio model satisfies two
conditions, implemented in `VFMValidation`:

1. **Solid-only domains**: `ValidateSolidDomains` walks every FE domain,
   ensuring each is an instance of `FESolidDomain`.
2. **Displacement coverage**: `ValidateDisplacementCounts` checks that both
   measured and virtual displacement containers contain exactly as many samples
   as the FE mesh has nodes.  If the counts differ, an informative message is
   generated.

These validations are invoked in `FEVFMTask::Init`.  Any failure is logged via
`feLogErrorEx`, aborting the task initialisation early.

---

## 6. Deformation Gradient Reconstruction

The current focus of the implementation is to reconstruct the deformation
gradient field from the measured displacements without altering FEBio’s mesh.
We introduce two helpers:

- `DeformationGradientField`: a container storing, for each element, the
  Gauss-point tensors \(\mathbf{F}\) (as `mat3d` instances).
- `VFMKinematics`: a utility namespace that computes \(\mathbf{F}\) by applying
  the Total Lagrangian relation.

### 6.1 Mathematical Background

Let \(\mathbf{u}_i\) denote the displacement vector at node \(i\) and
\(\nabla_X N_i\) the gradient of the element shape function with respect to the
reference configuration.  The deformation gradient at a Gauss point is given by

\[
  \mathbf{F} = \mathbf{I} + \sum_{i=1}^{n_{\mathrm{eln}}}
    \mathbf{u}_i \otimes \nabla_X N_i.
\]

The gradient \(\nabla_X N_i\) is obtained by transforming the natural shape
function derivatives \(\partial N_i / \partial (r, s, t)\) with the inverse of
the reference Jacobian \(\mathbf{J}_0^{-1}\):

\[
  \nabla_X N_i = \mathbf{J}_0^{-T}
    \begin{bmatrix}
      \partial N_i / \partial r \\
      \partial N_i / \partial s \\
      \partial N_i / \partial t
    \end{bmatrix}.
\]

FEBio provides `FESolidDomain::invjac0` which computes \(\mathbf{J}_0^{-1}\) at
the requested integration point.  In code we therefore:

1. Call `invjac0` to fetch `Ji`.
2. Loop over the element’s nodes; for each node \(i\), evaluate \(\nabla_X N_i\)
   by multiplying `Ji` with the natural derivatives stored in
   `el.Gr(n)`, `el.Gs(n)`, `el.Gt(n)`.
3. Accumulate the contribution \(\mathbf{u}_i \otimes \nabla_X N_i\) into \(\mathbf{F}\).
4. Add the identity matrix to recover the total deformation gradient.

When \(\mathbf{u}_i = \mathbf{0}\) for all nodes, the summation vanishes and
\(\mathbf{F} = \mathbf{I}\), which is an important consistency check that is
also verified in the debugging logs.

For diagnostic purposes `VFMKinematics` prints both the inverse Jacobian and
the resulting deformation gradient (with five decimal places) when FEBio runs
in debug mode:

```
Ji = [[0.50000 0.00000 0.00000] ...]
F(elem 1, gp 0) = [[1.00000 0.00000 0.00000] ...], det(F) = 1.00000
```

### 6.2 Implementation Details

- The nodal displacement lookup uses FEBio node IDs (`GetID()`), avoiding the
  common pitfall of assuming node IDs are strictly consecutive.
- Displacements are stored as `vec3d` instances, with one vector per element
  node for each element processed.
- If a displacement sample is missing or the determinant of \(\mathbf{F}\) is
  non-positive, the computation aborts and an explanatory message is returned.

The computed field is stored in `FEOptimizeDataVFM::DeformationGradients()` so
that subsequent VFM stages (stress evaluation, virtual work assembly, parameter
updates) can reuse it without recomputing.

---

## 7. Task Initialisation Summary

The steps performed in `FEVFMTask::Init` are now:

1. Parse optimisation parameters and displacement datasets (`FEVFMInput`).
2. Initialise parameter adapters (`FEOptimizeDataVFM::Init`).
3. Validate FE model topology and displacement coverage (`VFMValidation`).
4. Reconstruct Gauss-point deformation gradients once (`VFMKinematics`).

If all steps succeed the task is ready for the eventual Virtual Fields solve.

---

## 8. Next Steps

Future development will focus on:

- Incorporating the virtual displacement field into the VFM balance equations.
- Evaluating stresses \(\boldsymbol{\sigma}\) at Gauss points using the
  reconstructed \(\mathbf{F}\) and the material constitutive law.
- Assembling the virtual work equation
  \(\int_{\Omega} \boldsymbol{\sigma} : \delta \boldsymbol{\varepsilon}^* \, \mathrm{d}\Omega = \int_{\partial \Omega} \mathbf{t} \cdot \mathbf{v}^* \, \mathrm{d}\Gamma\)
  (or its inverse counterpart) to update material parameters.
- Completing the optimisation loop (`FEOptimizeDataVFM::Solve`) with appropriate
  objective functions and gradient computations.

These notes will be updated as new components are implemented.
