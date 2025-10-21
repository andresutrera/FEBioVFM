// #pragma once
// #include <string>
// #include "domain/vfm_types.hpp"
// #include "domain/vfm_containers.hpp"

// // small helpers
// std::string to_string(const vec3d& v);
// std::string to_string(const mat3d& m);

// // summaries
// std::string summary_measured(const MeasuredData& md);     // counts per time
// std::string summary_virtuals(const VirtualFields& vf);
// std::string summary_loads(const MeasuredLoad& ml);

// // optional heavy prints (guard with level)
// std::string dump_nodes(const NodalField<vec3d>& u, size_t max_items=20);
// std::string dump_elemgp(const RaggedElemField<mat3d>& F, size_t max_items=10);