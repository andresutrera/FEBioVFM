/**
 * @file MeasuredLoadContainer.h
 * @brief Storage helpers for measured surface loads across time steps.
 */
#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <cstddef>
#include <cmath>
#include <utility>

#include <FECore/vec3d.h>

/**
 * @brief Load sample applied to a named surface.
 */
struct SurfaceLoadSample
{
	std::string id; ///< User-supplied surface identifier (e.g. "left_grip").
	vec3d load = vec3d(0, 0, 0); ///< Load components Fx, Fy, Fz.

	SurfaceLoadSample() = default;
	SurfaceLoadSample(std::string surfaceId, vec3d value)
		: id(std::move(surfaceId)), load(value)
	{
	}
};

/**
 * @brief Container that stores surface loads for a single time step.
 */
class SurfaceLoadSet
{
public:
	/// Remove all stored loads.
	void Clear()
	{
		m_samples.clear();
		m_index.clear();
	}

	/// Add or overwrite a load entry for a surface.
	void Add(const std::string& surfaceId, const vec3d& value)
	{
		auto it = m_index.find(surfaceId);
		if (it != m_index.end())
		{
			m_samples[it->second].load = value;
			return;
		}

		const size_t idx = m_samples.size();
		m_samples.push_back(SurfaceLoadSample{ surfaceId, value });
		m_index.emplace(surfaceId, idx);
	}

	/// Number of stored load entries.
	size_t Size() const { return m_samples.size(); }

	/// Read-only access to the stored samples.
	const std::vector<SurfaceLoadSample>& Samples() const { return m_samples; }

	/// Attempt to access a particular surface load.
	const SurfaceLoadSample* Find(const std::string& surfaceId) const
	{
		auto it = m_index.find(surfaceId);
		return (it != m_index.end() ? &m_samples[it->second] : nullptr);
	}

	/// Try to get a load vector for a named surface.
	bool TryGet(const std::string& surfaceId, vec3d& value) const
	{
		const SurfaceLoadSample* sample = Find(surfaceId);
		if (sample == nullptr) return false;
		value = sample->load;
		return true;
	}

private:
	std::vector<SurfaceLoadSample> m_samples;
	std::unordered_map<std::string, size_t> m_index;
};

/**
 * @brief Timeline of measured loads keyed by time.
 */
class MeasuredLoadHistory
{
public:
	struct TimeStep
	{
		double time = 0.0;
		SurfaceLoadSet loads;
	};

	void Clear()
	{
		m_steps.clear();
	}

	TimeStep& AddStep(double time)
	{
		TimeStep step;
		step.time = time;
		m_steps.push_back(std::move(step));
		return m_steps.back();
	}

	void Reserve(size_t count) { m_steps.reserve(count); }

	size_t Steps() const { return m_steps.size(); }
	bool Empty() const { return m_steps.empty(); }

	TimeStep& operator[](size_t index) { return m_steps.at(index); }
	const TimeStep& operator[](size_t index) const { return m_steps.at(index); }

	TimeStep& StepAt(size_t index) { return m_steps.at(index); }
	const TimeStep& StepAt(size_t index) const { return m_steps.at(index); }

	TimeStep* FindStepByTime(double time, double tol = 1e-12)
	{
		for (auto& step : m_steps)
		{
			if (std::fabs(step.time - time) <= tol) return &step;
		}
		return nullptr;
	}

	const TimeStep* FindStepByTime(double time, double tol = 1e-12) const
	{
		for (const auto& step : m_steps)
		{
			if (std::fabs(step.time - time) <= tol) return &step;
		}
		return nullptr;
	}

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
	void ForEachLoad(Functor fn) const
	{
		for (const auto& step : m_steps)
		{
			for (const SurfaceLoadSample& entry : step.loads.Samples())
			{
				fn(step.time, entry);
			}
		}
	}

private:
	std::vector<TimeStep> m_steps;
};
