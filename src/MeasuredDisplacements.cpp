/**
 * @file MeasuredDisplacements.cpp
 * @brief Implementation of the displacement container.
 */

#include "MeasuredDisplacements.h"

#include <algorithm>

void DisplacementContainer::Clear()
{
	m_data.clear();
}

void DisplacementContainer::Add(int elemId, const std::array<double, 3>& disp)
{
	ElementDisplacement entry;
	entry.id = elemId;
	entry.displacement = disp;
	m_data.push_back(entry);
}

size_t DisplacementContainer::Size() const
{
	return m_data.size();
}

const ElementDisplacement* DisplacementContainer::Find(int elemId) const
{
	auto it = std::find_if(m_data.begin(), m_data.end(), [elemId](const ElementDisplacement& e) {
		return e.id == elemId;
	});
	return (it != m_data.end() ? &(*it) : nullptr);
}
