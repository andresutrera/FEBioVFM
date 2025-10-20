# FEBioVFM Refactor TODO

## Pain Points Observed
- `src/FEData.h:25` and `src/FEData.cpp:200` bundle parameter abstractions, runtime state, solver glue, and heavy numerical routines in a single class (`FEOptimizeDataVFM`), making it hard to reason about ownership and reuse.
- `src/VFM.cpp:435` keeps the task orchestration monolithic: input, validation, kinematics, stress evaluation, external work assembly, diagnostics, and exporting all live in one class (`FEVFMTask`).
- The XML reader in `src/FEVFMInput.cpp:24` writes directly into `FEOptimizeDataVFM`, so there is no decoupled representation of parsed input data.
- Multiple history containers (`DisplacementHistory`, `MeasuredLoadHistory`, `DeformationGradientHistory`, `StressHistory`, `FirstPiolaHistory`) repeat the same bookkeeping patterns across different headers (e.g., `src/DisplacementContainer.h:70`, `src/MeasuredLoadContainer.h:55`, `src/DeformationGradientField.h:55`), increasing maintenance cost.
- Utility concerns such as diagnostics (`src/VFM.cpp:350`), levmar hooks (`src/FEData.cpp:427`), and external work assembly (`src/VFM.cpp:667`) are interwoven with domain logic instead of being isolated services.

## Target Architecture
- **Domain layer**: immutable data descriptions (`MeasuredDataset`, `VirtualFieldSet`, `GaussPointField`) plus reusable timeline helpers.
- **State layer**: `VFMState` aggregates domain data (measured, virtual, calculated stresses) and parameter slots without embedding solver logic.
- **Services layer**: dedicated calculators (`KinematicsSolver`, `StressEvaluator`, `ExternalWorkAssembler`, `ResidualAssembler`) with clear inputs/outputs.
- **Optimization layer**: `VFMSolver` wraps levmar and only depends on service interfaces, not on FEBio task internals.
- **I/O layer**: `VFMInputParser` returns a plain `VFMConfig` struct; `VFMStateBuilder` maps that config into `VFMState`; `VFMExportSession` remains focused on writing.
- **Task layer**: `FEVFMTask` becomes a thin orchestrator that wires parser → builder → services → solver → export, delegating subtasks to the layers above.
- **Diagnostics**: `VFMDiagnostics` owns all logging helpers currently embedded in `src/VFM.cpp`.

## Workflow Orchestration Blueprint
1. `FEVFMTask::Init` creates an `VFMInputParser`, feeds the configuration file, and receives an immutable `VFMConfig`.
2. `FEVFMTask` hands the config to `VFMStateBuilder`, which produces a `VFMState` seeded with parameters, measured/virtual data, and empty result slots.
3. During setup, the task drives each service in sequence, passing references to the shared `VFMState`:
   - `KinematicsSolver` consumes measured displacement snapshots from the state and writes deformation gradients back into it.
   - `StressEvaluator` reads deformation history from the state and fills the stress timelines.
   - `ExternalWorkAssembler` operates on virtual displacements plus measured loads to populate the virtual work buffers.
   - `ResidualAssembler` binds the current state and exposes functors required by `VFMSolver`.
4. `FEVFMTask::Run` requests a `ParameterView` from the state, forwards it to `VFMSolver`, and applies the optimized values through the view once the solve finishes.
5. After optimization, `FEVFMTask` calls `VFMExportSession`, supplying read-only references to the data stored in `VFMState`, and `VFMDiagnostics` to log high-level summaries.
6. All cross-service communication goes through the state or explicit return values—services never invoke one another directly, keeping dependencies acyclic.

## Target Class Dependencies
- `FEVFMTask` depends on `VFMInputParser`, `VFMStateBuilder`, `VFMState`, `KinematicsSolver`, `StressEvaluator`, `ExternalWorkAssembler`, `ResidualAssembler`, `VFMSolver`, `VFMExportSession`, and `VFMDiagnostics`.
- `VFMStateBuilder` depends on `VFMConfig`, domain DTOs, and constructs `VFMState`; it does not see solver/services.
- `KinematicsSolver`, `StressEvaluator`, and `ExternalWorkAssembler` depend only on `VFMState` read/write views plus FEBio mesh/material interfaces; they never reference the task or solver.
- `ResidualAssembler` depends on `VFMState` views and exposes callbacks used by `VFMSolver`; no FEBio calls live inside the solver.
- `VFMSolver` depends on `ParameterSet`/`ResidualAssembler` abstractions and the levmar backend; it returns updated parameter values without touching the state directly.
- `VFMDiagnostics` accepts const references to `VFMState` and optional solver reports, keeping logging side effects isolated.
- `VFMExportSession` consumes the finalized `VFMState` but does not mutate it; the session remains decoupled from solver logic.

## Work Breakdown

### 1. Restructure Project Layout
- [ ] Introduce `src/domain/`, `src/services/`, `src/solver/`, `src/io/`, `src/task/` folders and move existing files accordingly.
- [ ] Split `src/FEData.h` into `src/domain/Parameter.h`, `src/state/VFMState.h`, and `src/solver/VFMSolver.h`; keep namespace-appropriate includes minimal.
- [ ] Update `CMakeLists.txt` to reflect the new file structure.

### 2. Normalize Data Containers
- [ ] Extract a generic `TimeHistory<T>` helper (`src/domain/TimeHistory.h`) and refactor `DisplacementHistory`, `MeasuredLoadHistory`, `DeformationGradientHistory`, `StressHistory`, and `FirstPiolaHistory` to build on it.
- [ ] Collapse virtual field wrappers into `VirtualField.h` and `VirtualFieldDefGrad.h` that share the same timeline primitives.
- [ ] Add lightweight value types (`NodeDisplacement`, `SurfaceLoadSample`, `GaussPointDeformation`, etc.) to a common header so services consume consistent DTOs.

### 3. Isolate Parameter Management
- [ ] Move `FEInputParameterVFM` and `FEModelParameterVFM` into `src/domain/Parameter.h/.cpp`.
- [ ] Replace the raw pointer vector (`m_Var`) with `std::unique_ptr<Parameter>` plus non-owning spans for FEBio APIs.
- [ ] Provide a `ParameterSet` facade that exposes `GetValues()`, `Assign(values)`, and metadata, removing direct access from `FEVFMTask::LogParameterTable`.

### 4. Build Explicit State + Builder
- [ ] Create `VFMState` that owns measured/virtual histories, deformation/stress results, and parameter set accessors.
- [ ] Implement `VFMStateBuilder` under `src/state/` to translate `VFMConfig` into a populated `VFMState` (inject dependencies instead of mutating `FEOptimizeDataVFM` during parse).
- [ ] Ensure state exposes const views for read-only services and limited mutating APIs for calculators.

### 5. Refactor Input Pipeline
- [ ] Rework `FEVFMInput` into `VFMInputParser` (in `src/io/`) that returns a `VFMConfig` (parameters, measured data, virtual fields, loads).
- [ ] Add unit tests for parser edge cases (missing nodes, duplicate surfaces, invalid parameter types).
- [ ] Make the builder responsible for reporting validation errors instead of mixing them into task logging.

### 6. Extract Computational Services
- [ ] Move `ComputeDeformationGradients` logic (`src/VFMKinematics.cpp:15`) into `KinematicsSolver` with a stateless API `bool run(const DisplacementField&, DeformationGradientField&, std::string&)`.
- [ ] Convert `VFMStress` into a `StressEvaluator` class that provides both Cauchy and first Piola evaluations; inject FEBio material dependencies explicitly.
- [ ] Introduce `ExternalWorkAssembler` for the code in `src/VFM.cpp:650`, allowing it to be reused by solver stages and tested independently.
- [ ] Add `ResidualAssembler` to wrap the logic currently in `FEOptimizeDataVFM::ComputeInternalWork`/`ComputeResidual`, returning structured results (flattened vectors + mapping metadata).

### 7. Wrap Optimization Logic
- [ ] Encapsulate levmar callbacks (`src/FEData.cpp:427`) inside `VFMSolver`, exposing `std::vector<double> solve(const ParameterSet&, const ResidualAssembler&)`.
- [ ] Allow solver options (max iterations, tolerances) to be configured via `VFMConfig`.
- [ ] Provide hooks to plug alternative optimizers in the future by defining a `IVFMOptimizer` interface.

### 8. Slim the FEBio Task
- [ ] Replace direct method implementations in `FEVFMTask` with calls into the new services; the class should mainly translate FEBio lifecycle events into service invocations.
- [ ] Move all diagnostic/logging helpers to `VFMDiagnostics` and inject it where needed.
- [ ] Ensure `FEVFMTask::Run` reads final parameters and exports through well-defined APIs (`ParameterSet`, `VFMExportSession`), without peeking into internal vectors.

### 9. Improve Testing & Documentation
- [ ] Add unit tests per service (kinematics, stress, external work, residual, solver) covering the refactored APIs.
- [ ] Update `README.md` and `FEBioVFM.md` to describe the new module layout and responsibilities.
- [ ] Ship Doxygen groups aligned with the new folder structure so generated docs reflect the cleaned class division.

### 10. Migration & Cleanup
- [ ] Remove legacy fields/methods rendered obsolete after the split (e.g., `FEOptimizeDataVFM::Solve` placeholder).
- [ ] Audit includes for minimal exposure, forward declare FEBio types where possible.
- [ ] Run clang-format / clang-tidy once the reorganization is complete to standardize style (keep tooling config in repo).
