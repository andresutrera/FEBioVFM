src/
  domain/                # pure data + tiny helpers (no FEBio)
    vfm_types.hpp        # vec3d, mat3d, indices
    vfm_containers.hpp   # NodalField, RaggedElemField, frames, TimeSeries
    timeline.hpp         # small time utilities (indexing, maps)
    measured_dataset.hpp # typedefs over containers (MeasuredData, MeasuredLoad)
    virtual_set.hpp      # VirtualFields
    stresses.hpp         # Stresses

  state/                 # mutable aggregate, no algorithms
    vfm_state.hpp        # VFMState { measured, virtuals, stresses, params }

  services/              # pure functions or small classes
    kinematics.hpp       # interfaces only
    kinematics.cpp       # F from u using mesh + shapes
    stress_eval.hpp
    stress_eval.cpp
    ext_work.hpp
    ext_work.cpp
    residual.hpp
    residual.cpp

  optimization/
    vfm_solver.hpp       # VFMSolver interface (run, options)
    vfm_solver_levmar.cpp# levmar glue; depends on services/* only

  io/
    xml_reader.hpp       # VFMXmlReader → XMLInput (sparse DTOs)
    xml_reader.cpp
    loader.hpp           # VFMLoader: DTO → containers in VFMState (nodal u, loads)
    loader.cpp
    export.hpp           # VFMExportSession
    export.cpp

  build/                 # mesh facts, IDs, gp shapes (no algorithms)
    mesh_maps.hpp        # nodeId2idx, elemNodes, nGPperElem
    mesh_maps.cpp
    state_builder.hpp    # VFMStateBuilder: takes XMLInput+MeshMaps → VFMState
    state_builder.cpp

  task/                  # FEBio-facing thin layer
    fevfm_task.hpp
    fevfm_task.cpp       # wires parser → builder → services → solver → export

  diag/
    diagnostics.hpp      # logging API
    diagnostics.cpp

  plugin/                # FEBio plugin entry points only
    register.cpp




Dependency rules

domain ← nothing.

state ← domain.

services ← domain, state (read-only), build (for mesh shapes).

optimization ← services, state, diag.

io:

xml_reader ← nothing from project.

loader ← domain, state.

export ← state.

build ← domain, io::xml_reader (reads DTO types only).

task ← io, build, services, optimization, diag, state.

plugin ← task only.

Targets (CMake)

vfm_domain → src/domain/*

vfm_state → src/state/*

vfm_services → src/services/*

vfm_optimization → src/optimization/*

vfm_io → src/io/*

vfm_build → src/build/*

vfm_task → src/task/*

vfm_plugin (SHARED) → src/plugin/*
Link order: vfm_plugin → vfm_task → (vfm_optimization vfm_services vfm_build vfm_io vfm_state vfm_domain diag).