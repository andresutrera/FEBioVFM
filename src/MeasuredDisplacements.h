/**
 * @file MeasuredDisplacements.h
 * @brief Container for nodal or elemental displacement definitions.
 *
 * The Virtual Fields Method compares simulated responses against externally
 * measured data and leverages admissible virtual fields. This helper class
 * stores the mapping between an element ID and an associated displacement
 * vector (ux, uy, uz) so that downstream stages can access the values without
 * re-parsing the XML input.
 */
#pragma once

#include <array>
#include <vector>

/**
 * @brief Represents displacement data associated with a single finite element.
 */
struct ElementDisplacement
{
	int id = -1;                                  ///< Element identifier.
	std::array<double, 3> displacement{{0, 0, 0}}; ///< Measured ux, uy, uz tuple.
};

/**
 * @brief Stores measured displacement vectors indexed by element ID.
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
	void Clear();

	/**
	 * @brief Append a new element displacement measurement.
	 * @param elemId Element identifier from the XML input.
	 * @param disp Measured ux, uy, uz tuple.
	 */
	void Add(int elemId, const std::array<double, 3>& disp);

	/**
	 * @brief Number of stored displacement samples.
	 */
	size_t Size() const;

	/**
	 * @brief Provide read-only access to all stored displacement samples.
	 */
	const std::vector<ElementDisplacement>& Samples() const { return m_data; }

	/**
	 * @brief Find the measurement associated with a given element ID.
	 * @param elemId Element identifier to look up.
	 * @return Pointer to the measurement or nullptr when not present.
	 */
	const ElementDisplacement* Find(int elemId) const;

private:
	std::vector<ElementDisplacement> m_data; ///< Insertion-ordered displacement samples.
};
