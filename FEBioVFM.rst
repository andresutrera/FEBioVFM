=============================================
FEBio Virtual Fields Method Plug-in Overview
=============================================

.. contents::
   :depth: 2
   :local:

Purpose and Scope
=================

The FEBio Virtual Fields Method (VFM) plug-in augments FEBio with an inverse
identification workflow that fits constitutive parameters against full-field
deformation measurements.  This document provides an implementation-oriented
reference covering:

* the data model exposed by the plug-in,
* the numerical formulation of the residual and cost functions,
* the decomposition of the virtual work equation implemented in the code base,
* the optimisation strategy based on the ``dlevmar_bc_dif`` solver,
* file export facilities and logging conventions for debugging.

Although stored as a ``.md`` file for compatibility with existing tooling, the
syntax below follows reStructuredText conventions so that it can be rendered
by Sphinx or other docutils-driven pipelines.

Architecture at a Glance
========================

The plug-in is structured around an FEBio task (``FEVFMTask``) that orchestrates
the workflow and a data container (``FEOptimizeDataVFM``) that owns optimisation
state.

``FEVFMTask`` sequence
----------------------

1. **Input parsing** via ``FEOptimizeDataVFM::Input``.
2. **State initialisation** using ``FEOptimizeDataVFM::Init`` to register
   parameters and cache their initial values.
3. **Model validation** (solid domains only and consistent displacement/load
   coverage).
4. **Kinematic reconstruction** – builds deformation gradient histories for
   both measured and virtual fields.
5. **Stress history recovery** – computes Cauchy and first Piola-Kirchhoff
   stresses for every Gauss point and time step.
6. **External virtual work assembly** – integrates the measured tractions over
   named surface sets.
7. **Residual assembly** and **optimisation** via the Levenberg–Marquardt (LM)
   solver with bound constraints.
8. **Diagnostics and export** – logs the parameter evolution, residual cost,
   and writes an ``.xplt`` snapshot containing all histories evaluated at the
   converged parameters.

The optimisation container exposes convenience methods that the task calls:

* ``SetParameterVector`` and ``GetParameterVector`` for synchronising the FEBio
  model with the optimiser,
* ``RebuildStressHistories`` to recompute stresses whenever parameters change,
* ``AssembleResidual`` to evaluate the virtual work imbalance,
* ``MinimizeResidualWithLevmar`` to launch the LM solver using the residual as
  a zero-target least-squares system.

Data Model
==========

Parameters
----------

Scalar optimisation variables descend from ``FEInputParameterVFM``.  Each
instance stores:

* ``InitValue`` – the starting value written to the FEBio model,
* ``MinValue`` and ``MaxValue`` – bounds fed to ``dlevmar_bc_dif``,
* ``ScaleFactor`` – retained for future normalisation options,
* a string ``name`` to resolve the underlying FEBio parameter.

``FEModelParameterVFM`` is the predominant implementation, resolving a
``FEParamDouble`` from the model and caching a mutable pointer to its raw value.

Histories and Fields
--------------------

Displacement data and their derived quantities are organised into time-indexed
histories.  The following containers are used extensively:

* ``DisplacementHistory`` – a sequence of ``NodeDisplacement`` samples
  (per-node vectors) keyed by time.
* ``DeformationGradientHistory`` – holds ``DeformationGradientField`` instances,
  each mapping element identifiers to per-Gauss-point gradients ``F``.
* ``StressHistory`` – analogous container for Cauchy stress tensors ``σ``.
* ``FirstPiolaHistory`` – per-Gauss-point first Piola–Kirchhoff tensors ``P``.
* ``VirtualDisplacementCollection`` and
  ``VirtualDeformationGradientCollection`` – indexed arrays of user-supplied
  virtual fields, each with an optional identifier for logging/export.
* ``VirtualExternalWorkHistory`` – stores the scalar virtual external work
  sampled at the displacement time grid.

Measured load histories are keyed by surface names that correspond to FEBio
surface or facet sets.  During assembly, surfaces are resolved to node lists
through the FEBio mesh API.

Mathematical Formulation
========================

Virtual Fields Equation
-----------------------

For each virtual displacement field :math:`\hat{\mathbf{u}}^{(k)}` and time
step :math:`t_j` the VFM imposes equilibrium of internal and external virtual
work:

.. math::

   \mathcal{R}_{k,j} =
     \int_{\Omega_0} \mathbf{P}(t_j) : \nabla_X \hat{\mathbf{u}}^{(k)} \, \mathrm{d}V
     - \int_{\partial \Omega_0} \hat{\mathbf{u}}^{(k)} \cdot \bar{\mathbf{T}}(t_j)\, \mathrm{d}S
     = 0.

Here :math:`\mathbf{P}` is the first Piola–Kirchhoff stress derived from the
measured deformation gradient and current material parameters, and
:math:`\bar{\mathbf{T}}` represents the applied traction vector reconstructed
from measured surface loads.

The discrete residual implemented in ``FEOptimizeDataVFM::AssembleResidual``
evaluates

.. math::

   \mathcal{R}_{k,j}
   = \sum_{e \in \mathcal{T}} \sum_{q=1}^{n_q}
       w_{e,q} \left[ \mathbf{P}_{e,q}(t_j) : \mathbf{G}_{e,q}^{(k)} \right]
     - \mathcal{W}_{k}^{\text{ext}}(t_j),

where

* :math:`w_{e,q}` are the precomputed element integration weights (`detJ0 *
  gaussWeight`),
* :math:`\mathbf{G}_{e,q}^{(k)} = \nabla_X \hat{\mathbf{u}}^{(k)}` is recovered
  by removing the identity from the stored virtual deformation gradient
  (``VirtualGradientFromDeformation``),
* :math:`\mathcal{W}_{k}^{\text{ext}}(t_j)` is the virtual work of measured
  tractions, calculated by distributing surface forces to the corresponding
  nodes and projecting them through the virtual displacements.

Cost Functional
---------------

The optimisation cost is the standard least-squares metric:

.. math::

   J(\boldsymbol{p}) = \frac{1}{2}\, \mathbf{r}(\boldsymbol{p})^\mathsf{T}
   \mathbf{r}(\boldsymbol{p}),

with :math:`\mathbf{p}` denoting the stacked model parameters and
:math:`\mathbf{r}` the residual vector obtained by concatenating
:math:`\mathcal{R}_{k,j}` for all virtual fields and time steps.

Implementation details:

* ``FEVFMTask::Run`` logs ``J`` both before and after the LM solve.
* The residual is always reassembled with the current parameter set before
  evaluating the cost.

Stress Reconstruction
---------------------

The stress pipeline resides in ``VFMStress`` and ``FEOptimizeDataVFM``:

1. **Deformation Gradient** – ``VFMKinematics::ComputeDeformationGradients``
   evaluates the total Lagrangian expression

   .. math::

      \mathbf{F}_{e,q} = \mathbf{I} + \sum_{a=1}^{n_{\text{node}}}
        \mathbf{u}_a \otimes \nabla_X N_a(\xi_q),

   extracting nodal displacements from the measured history.

2. **Cauchy Stress** – ``VFMStress::ComputeCauchyStress`` calls into FEBio to
   update material point data and retrieves :math:`\boldsymbol{\sigma}`.

3. **First Piola–Kirchhoff Stress** – computed analytically from the Cauchy
   tensor and deformation gradient via

   .. math::

      \mathbf{P} = J\, \boldsymbol{\sigma} \mathbf{F}^{-\mathsf{T}},

   where :math:`J = \det{\mathbf{F}}`.

4. **Virtual Deformation Gradients** – when virtual fields are supplied as
   displacement histories, their gradients are computed using the same
   kinematic routine as the measured data.  The tensor
   :math:`\mathbf{G} = \mathbf{F}_\text{virtual} - \mathbf{I}` is used in the
   internal virtual work calculation.

Optimisation Strategy
=====================

Residual Callback
-----------------

``FEOptimizeDataVFM::MinimizeResidualWithLevmar`` wraps the residual assembly
into a callback compatible with ``dlevmar_bc_dif``:

* Parameters are copied from the LM-provided vector and applied to the FEBio
  model without touching the load step state.
* ``AssembleResidual(parameters, false, residual, error)`` evaluates the
  residual directly (no temporary stress rebuild) because the current
  parameter set is already synchronised.
* The residual is written to the buffer supplied by LM.  Any failure in the
  process (invalid bounds, stress reconstruction error, etc.) sets a flag that
  aborts the optimisation and restores the original parameter vector and stress
  histories.

Solver Configuration
--------------------

Key options passed to ``dlevmar_bc_dif``:

* **Bounds** – extracted from ``FEInputParameterVFM::MinValue`` and
  ``::MaxValue`` for every parameter; validation ensures ``min ≤ max``.
* **Initial guess** – uses the current parameter vector as returned by the
  FEBio model (after initialisation or after a prior solve).
* **Iteration limit** – defaults to 100 when the caller supplies a non-positive
  limit; the task currently passes zero to select this default.
* **Options array** – ``[ LM_INIT_MU, 1e-12, 1e-12, 1e-12, LM_DIFF_DELTA ]``,
  providing tight tolerance thresholds and the default finite-difference step
  for the Jacobian approximation.
* **Workspace** – allocated dynamically using ``LM_BC_DIF_WORKSZ`` for the
  requested problem size.

The solver returns the number of LM iterations (non-negative on success), which
is stored in ``FEOptimizeDataVFM::m_niter``.  The raw ``info`` vector is passed
back to callers on demand for advanced diagnostics.

Logging and Diagnostics
=======================

* ``LogParameterValues`` prints current parameter values along with their
  identifiers in debug mode and (after optimisation) in release mode.
* ``FEVFMTask::Run`` reports the initial and final costs, the LM termination
  metrics, and the final parameter vector.
* ``LogStressDiagnostics`` is invoked only after successful optimisation,
  guaranteeing that stresses correspond to the converged parameters.
* External virtual work integrations provide detailed messages (surface names,
  nodal forces, virtual displacements) when debug logging is enabled.

Exporting Results
=================

``FEVFMTask::ExportState`` writes an FEBio ``.xplt`` file reflecting the
current model state:

* measured displacements and deformation gradients,
* virtual displacement and deformation gradient histories,
* recovered Cauchy and first Piola–Kirchhoff stresses,
* virtual external work histories.

During ``Init`` an initial snapshot is exported using the file specified on the
command line (``-task="VFM" path/to/VFMData.feb``).  After optimisation,
``Run`` reuses the same base path to overwrite the ``.xplt`` file with the
stresses and histories that correspond to the final parameter values.

Extending the Plug-in
=====================

Several extension points remain intentionally lightweight:

* **Forward solver integration** – ``FEOptimizeDataVFM::FESolve`` is currently
  a stub.  Hooking into FEBio's non-linear solver will enable iterative forward
  updates between LM iterations, which is necessary for strongly non-linear
  constitutive laws.
* **Analytical Jacobians** – ``dlevmar_bc_dif`` uses finite differences.
  Implementing ``dlevmar_bc_der`` with custom Jacobian code could improve
  convergence in challenging problems.
* **Regularisation** – ``MinimizeResidualWithLevmar`` focuses on pure least
  squares.  Optional Tikhonov or Bayesian priors can be incorporated by
  augmenting the residual vector with synthetic measurements.
* **Adaptive tolerances** – the current tolerance triplet ``1e-12`` is chosen
  to favour accuracy.  Problem-specific tuning hooks can be added to the
  ``FEVFMTask`` options once user-facing configuration is defined.

Testing and Validation Strategy
===============================

Unit-style validation is split across deterministic checks:

* **Parser regression tests** – sample ``VFMData.feb`` files exercise the XML
  reader and populate all histories.
* **Kinematics sanity checks** – deformation gradients reconstructed from
  identity deformation fields match ``F = I`` to machine precision.
* **Residual invariance** – for a linear elastic benchmark with analytical
  solution, the residual norms match reference values.
* **Optimization smoke tests** – the LM wrapper is executed on synthetic
  problems to verify bound enforcement and stress history rebuild logic.

When integrating with FEBio, set ``-g`` for additional log output, and monitor
the exported ``.xplt`` file in PostView to confirm that virtual fields and
stresses evolve as expected.

Appendix: Symbol Reference
==========================

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Symbol
     - Meaning
   * - :math:`\mathbf{u}`
     - Measured displacement vector in the reference configuration.
   * - :math:`\hat{\mathbf{u}}`
     - Virtual displacement field supplied by the user.
   * - :math:`\mathbf{F}`
     - Deformation gradient, :math:`\nabla_X \mathbf{x}`.
   * - :math:`\mathbf{P}`
     - First Piola–Kirchhoff stress, :math:`J\, \boldsymbol{\sigma} \mathbf{F}^{-\mathsf{T}}`.
   * - :math:`\boldsymbol{\sigma}`
     - Cauchy stress tensor (true stress).
   * - :math:`\mathbf{r}`
     - Stacked virtual work residual vector.
   * - :math:`J(\boldsymbol{p})`
     - Least-squares cost function, :math:`\frac{1}{2} \mathbf{r}^\mathsf{T} \mathbf{r}`.
   * - :math:`w_{e,q}`
     - Quadrature weight (Gauss weight times reference Jacobian determinant).
   * - :math:`\mathcal{W}^{\text{ext}}`
     - External virtual work accumulated from measured tractions.

This document should serve as the canonical reference while the plug-in evolves
toward a full-fledged VFM identification tool integrated within FEBio.
