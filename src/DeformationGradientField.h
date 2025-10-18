/**
 * @file DeformationGradientField.h
 * @brief Storage for gauss-point deformation gradient tensors.
 */
#pragma once

#include <vector>
#include <algorithm>
#include <utility>

#include <FECore/mat3d.h>

/**
 * @brief Holds deformation gradients for a single finite element.
 */
struct GaussPointDeformation
{
	int elementId = -1;              ///< Element identifier (`GetID()` from FEBio).
	std::vector<mat3d> gradients;    ///< Deformation gradient tensors per Gauss point.
};

/**
 * @brief Aggregates deformation gradients for all elements in the analysed mesh.
 */
class DeformationGradientField
{
public:
	/// Remove all stored deformation gradients.
	inline void Clear() { m_data.clear(); }

	/// Append gradients for a new element.
	inline void Add(GaussPointDeformation entry) { m_data.push_back(std::move(entry)); }

	/// Number of elements with stored gradients.
	inline size_t Size() const { return m_data.size(); }

	/// Read/write access to the stored data.
	inline std::vector<GaussPointDeformation>& Data() { return m_data; }

	/// Read-only access to the stored data.
	inline const std::vector<GaussPointDeformation>& Data() const { return m_data; }

	/// Find gradients for a particular element ID.
	inline const GaussPointDeformation* Find(int elementId) const
	{
		auto it = std::find_if(m_data.begin(), m_data.end(), [elementId](const GaussPointDeformation& e) {
			return e.elementId == elementId;
		});
		return (it != m_data.end() ? &(*it) : nullptr);
	}

private:
	std::vector<GaussPointDeformation> m_data;
};
