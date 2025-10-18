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

/**
 * @brief Timeline wrapper around deformation gradients per time step.
 */
class DeformationGradientHistory
{
public:
	struct TimeStep
	{
		double time = 0.0;
		DeformationGradientField field;
	};

	void Clear()
	{
		m_steps.clear();
		m_active = 0;
	}

	TimeStep& AddStep(double time)
	{
		TimeStep step;
		step.time = time;
		m_steps.push_back(std::move(step));
		return m_steps.back();
	}

	size_t Steps() const { return m_steps.size(); }
	bool Empty() const { return m_steps.empty(); }

	TimeStep& StepAt(size_t index) { return m_steps.at(index); }
	const TimeStep& StepAt(size_t index) const { return m_steps.at(index); }

	void SetActiveStepByIndex(size_t index)
	{
		if (m_steps.empty())
		{
			m_active = 0;
		}
		else
		{
			m_active = (index < m_steps.size()) ? index : (m_steps.size() - 1);
		}
	}

	TimeStep& ActiveStep()
	{
		if (m_steps.empty())
		{
			m_steps.push_back(TimeStep{});
		}
		if (m_active >= m_steps.size()) m_active = 0;
		return m_steps[m_active];
	}

	const TimeStep& ActiveStep() const
	{
		if (m_steps.empty())
		{
			static const TimeStep empty;
			return empty;
		}
		size_t idx = (m_active < m_steps.size()) ? m_active : 0;
		return m_steps[idx];
	}

	std::vector<TimeStep>& StepsRef() { return m_steps; }
	const std::vector<TimeStep>& StepsRef() const { return m_steps; }

private:
	std::vector<TimeStep> m_steps;
	size_t m_active = 0;
};
