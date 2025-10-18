/**
 * @file DisplacementContainer.h
 * @brief Container for nodal displacement definitions.
 *
 * The Virtual Fields Method compares simulated responses against externally
 * measured data and leverages admissible virtual fields. This helper class
 * stores the mapping between a mesh node and an associated displacement vector
 * (ux, uy, uz) so that downstream stages can access the values without
 * re-parsing the XML input.
 */
#pragma once

#include <array>
#include <vector>
#include <cstddef>
#include <algorithm>

/**
 * @brief Represents displacement data associated with a single mesh node.
 */
struct NodeDisplacement
{
	int id = -1;                                  ///< Node identifier (1-based as in FEBio input).
	std::array<double, 3> displacement{{0, 0, 0}}; ///< Measured ux, uy, uz tuple.
};

/**
 * @brief Stores measured displacement vectors indexed by node ID.
 *
 * The class provides simple insert and query helpers that are used by the XML
 * parser as well as future inverse-solver stages. Data is kept in insertion
 * order which matches the sequence found in the XML file.
 */
class DisplacementContainer
{
public:
	/**
	 * @brief Remove all stored displacement entries.
	 */
	inline void Clear() { m_data.clear(); }

	/**
	 * @brief Append a new nodal displacement measurement.
	 * @param nodeId Node identifier from the XML input.
	 * @param disp Measured ux, uy, uz tuple.
	 */
	inline void Add(int nodeId, const std::array<double, 3>& disp)
	{
		NodeDisplacement entry;
		entry.id = nodeId;
		entry.displacement = disp;
		m_data.push_back(entry);
	}

	/**
	 * @brief Number of stored displacement samples.
	 */
	inline size_t Size() const { return m_data.size(); }

	/**
	 * @brief Provide read-only access to all stored displacement samples.
	 */
	const std::vector<NodeDisplacement>& Samples() const { return m_data; }

	/**
	 * @brief Find the measurement associated with a given element ID.
	 * @param elemId Element identifier to look up.
	 * @return Pointer to the measurement or nullptr when not present.
	 */
	inline const NodeDisplacement* Find(int nodeId) const
	{
		auto it = std::find_if(m_data.begin(), m_data.end(), [nodeId](const NodeDisplacement& e) {
			return e.id == nodeId;
		});
		return (it != m_data.end() ? &(*it) : nullptr);
	}

private:
	std::vector<NodeDisplacement> m_data; ///< Insertion-ordered displacement samples.
};
