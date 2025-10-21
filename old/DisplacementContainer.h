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
#include <cmath>

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
	 * @brief Find the measurement associated with a given node ID.
	 * @param nodeId Node identifier to look up.
	 * @return Pointer to the measurement or nullptr when not present.
	 */
	inline const NodeDisplacement* Find(int nodeId) const
	{
		auto it = std::find_if(m_data.begin(), m_data.end(), [nodeId](const NodeDisplacement& e) {
			return e.id == nodeId;
		});
		return (it != m_data.end() ? &(*it) : nullptr);
	}

	/**
	 * @brief Try to obtain the displacement for a node.
	 * @param nodeId Node identifier (1-based) to look up.
	 * @param disp Filled with the displacement when found.
	 * @return true when the displacement is available.
	 */
	inline bool TryGet(int nodeId, std::array<double, 3>& disp) const
	{
		const NodeDisplacement* entry = Find(nodeId);
		if (entry == nullptr) return false;
		disp = entry->displacement;
		return true;
	}

private:
	std::vector<NodeDisplacement> m_data; ///< Insertion-ordered displacement samples.
};

/**
 * @brief Time history wrapper that stores displacement containers for multiple steps.
 */
class DisplacementHistory
{
public:
	struct TimeStep
	{
		double time = 0.0;
		DisplacementContainer displacements;
	};

	/// Remove all steps.
	inline void Clear()
	{
		m_steps.clear();
		m_active = 0;
	}

	/// Append a new time step (steps are stored in insertion order).
	inline TimeStep& AddStep(double time)
	{
		TimeStep step;
		step.time = time;
		m_steps.push_back(std::move(step));
		return m_steps.back();
	}

	/// Remove all steps and reserve storage.
	inline void Reserve(size_t count) { m_steps.reserve(count); }

	/// Number of stored time steps.
	inline size_t Steps() const { return m_steps.size(); }

	/// True when no time steps are present.
	inline bool Empty() const { return m_steps.empty(); }

	/// Random access to a time step by index.
	inline TimeStep& operator[](size_t index) { return m_steps.at(index); }
	inline const TimeStep& operator[](size_t index) const { return m_steps.at(index); }

	/// Find an existing step by time value within a tolerance.
	inline TimeStep* FindStepByTime(double time, double tol = 1e-12)
	{
		for (auto& step : m_steps)
		{
			if (std::fabs(step.time - time) <= tol) return &step;
		}
		return nullptr;
	}

	inline const TimeStep* FindStepByTime(double time, double tol = 1e-12) const
	{
		for (const auto& step : m_steps)
		{
			if (std::fabs(step.time - time) <= tol) return &step;
		}
		return nullptr;
	}

	/// Direct access to the underlying storage.
	inline std::vector<TimeStep>& StepsRef() { return m_steps; }
	inline const std::vector<TimeStep>& StepsRef() const { return m_steps; }

	using iterator = std::vector<TimeStep>::iterator;
	using const_iterator = std::vector<TimeStep>::const_iterator;

	inline iterator begin() { return m_steps.begin(); }
	inline iterator end() { return m_steps.end(); }
	inline const_iterator begin() const { return m_steps.begin(); }
	inline const_iterator end() const { return m_steps.end(); }
	inline const_iterator cbegin() const { return m_steps.cbegin(); }
	inline const_iterator cend() const { return m_steps.cend(); }

	template <typename Functor>
	void ForEachTime(Functor fn) const
	{
		for (const auto& step : m_steps)
		{
			fn(step.time);
		}
	}

	template <typename Functor>
	void ForEachMeasurement(Functor fn) const
	{
		for (const auto& step : m_steps)
		{
			for (const NodeDisplacement& entry : step.displacements.Samples())
			{
				fn(step.time, entry);
			}
		}
	}

private:
	std::vector<TimeStep> m_steps;
	size_t m_active = 0;
};
