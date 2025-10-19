/**
 * @file FirstPiolaField.h
 * @brief Storage helpers for Gauss-point first Piola-Kirchhoff stresses across time steps.
 */
#pragma once

#include <vector>
#include <algorithm>
#include <utility>

#include <FECore/mat3d.h>

/**
 * @brief Holds first Piola-Kirchhoff stress tensors for a single finite element.
 */
struct GaussPointFirstPiola
{
	int elementId = -1;              ///< Element identifier (`GetID()` from FEBio).
	std::vector<mat3d> stresses;     ///< First Piola-Kirchhoff stress tensors per Gauss point.
};

/**
 * @brief Aggregates first Piola stresses for all elements in the analysed mesh.
 */
class FirstPiolaField
{
public:
	/// Remove all stored stresses.
	inline void Clear() { m_data.clear(); }

	/// Append stresses for a new element.
	inline void Add(GaussPointFirstPiola entry) { m_data.push_back(std::move(entry)); }

	/// Number of elements with stored stresses.
	inline size_t Size() const { return m_data.size(); }

	/// Read/write access to the stored data.
	inline std::vector<GaussPointFirstPiola>& Data() { return m_data; }

	/// Read-only access to the stored data.
	inline const std::vector<GaussPointFirstPiola>& Data() const { return m_data; }

	/// Find stresses for a particular element ID.
	inline const GaussPointFirstPiola* Find(int elementId) const
	{
		auto it = std::find_if(m_data.begin(), m_data.end(), [elementId](const GaussPointFirstPiola& e) {
			return e.elementId == elementId;
		});
		return (it != m_data.end() ? &(*it) : nullptr);
	}

private:
	std::vector<GaussPointFirstPiola> m_data;
};

/**
 * @brief Timeline wrapper around first Piola stresses per time step.
 */
class FirstPiolaHistory
{
public:
	struct TimeStep
	{
		double time = 0.0;
		FirstPiolaField field;
	};

	void Clear()
	{
		m_steps.clear();
		m_active = 0;
	}

	void Reserve(size_t count) { m_steps.reserve(count); }

	TimeStep& AddStep(double time)
	{
		TimeStep step;
		step.time = time;
		m_steps.push_back(std::move(step));
		return m_steps.back();
	}

	size_t Steps() const { return m_steps.size(); }
	bool Empty() const { return m_steps.empty(); }

	TimeStep& operator[](size_t index) { return m_steps.at(index); }
	const TimeStep& operator[](size_t index) const { return m_steps.at(index); }

	std::vector<TimeStep>& StepsRef() { return m_steps; }
	const std::vector<TimeStep>& StepsRef() const { return m_steps; }

	using iterator = std::vector<TimeStep>::iterator;
	using const_iterator = std::vector<TimeStep>::const_iterator;

	iterator begin() { return m_steps.begin(); }
	iterator end() { return m_steps.end(); }
	const_iterator begin() const { return m_steps.begin(); }
	const_iterator end() const { return m_steps.end(); }
	const_iterator cbegin() const { return m_steps.cbegin(); }
	const_iterator cend() const { return m_steps.cend(); }

	template <typename Functor>
	void ForEachTime(Functor fn) const
	{
		for (const auto& step : m_steps)
		{
			fn(step.time);
		}
	}

	template <typename Functor>
	void ForEachStress(Functor fn) const
	{
		for (const auto& step : m_steps)
		{
			for (const GaussPointFirstPiola& entry : step.field.Data())
			{
				fn(step.time, entry);
			}
		}
	}

private:
	std::vector<TimeStep> m_steps;
	size_t m_active = 0;
};
