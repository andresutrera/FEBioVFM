/**
 * @file VFMValidation.h
 * @brief Validation helpers for FEBio models used by the VFM plugin.
 */
#pragma once

#include <string>
#include <algorithm>
#include <cmath>
#include <unordered_set>

#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/FESolidDomain.h>

#include "FEData.h"

/**
 * @brief Encapsulates reusable checks that must pass before running the VFM task.
 *
 * The validator is intentionally lightweight so additional rules can be added
 * without cluttering the task or optimisation code.
 */
class VFMValidation
{
public:
	/**
	 * @brief Validate that the supplied FEModel only contains solid domains.
	 * @param fem FEBio model to inspect.
	 * @param errorMessage Populated with a description when the check fails.
	 * @return true if all domains are solid.
	 */
	static inline bool ValidateSolidDomains(FEModel& fem, std::string& errorMessage)
	{
		FEMesh& mesh = fem.GetMesh();
		const int domainCount = mesh.Domains();
		for (int i = 0; i < domainCount; ++i)
		{
			FEDomain& dom = mesh.Domain(i);
			if (dynamic_cast<FESolidDomain*>(&dom) == nullptr)
			{
				errorMessage = "The Virtual Fields Method only supports solid domains; found non-solid domain \"" + dom.GetName() + "\".";
				return false;
			}
		}
		return true;
	}

	/**
	 * @brief Validate that displacement counts align with the mesh size.
	 * @param fem FEBio model whose mesh provides the node count.
	 * @param data Optimisation container with measured and virtual displacement sets.
	 * @param errorMessage Populated with a description when the check fails.
	 * @return true when both displacement sets contain one entry per mesh node.
	 */
	static inline bool ValidateDisplacementCounts(FEModel& fem, const FEOptimizeDataVFM& data, std::string& errorMessage)
	{
		const int nodeCount = fem.GetMesh().Nodes();

		const auto& measuredHistory = data.MeasuredHistory();
		const auto& virtualHistory = data.VirtualHistory();
		

		if (measuredHistory.Empty())
		{
			errorMessage = "Measured displacement history is empty.";
			return false;
		}

		if (virtualHistory.Empty())
		{
			errorMessage = "Virtual displacement history is empty.";
			return false;
		}

		if (measuredHistory.Steps() != virtualHistory.Steps())
		{
			errorMessage = "Measured and virtual displacement histories contain a different number of time steps.";
			return false;
		}

		if (measuredHistory.Steps() != data.DeformationHistory().Steps())
		{
			errorMessage = "Deformation gradient history does not match displacement histories.";
			return false;
		}


		for (size_t i = 0; i < measuredHistory.Steps(); ++i)
		{
			const auto& measStep = measuredHistory[i];
			const auto& virtStep = virtualHistory[i];
			const size_t countMeasured = measStep.displacements.Size();
			const size_t countVirtual = virtStep.displacements.Size();

			if ((int)countMeasured != nodeCount)
			{
				errorMessage = "Measured displacement count at time " + std::to_string(measStep.time) + " (" + std::to_string(countMeasured) + ") does not match mesh node count (" + std::to_string(nodeCount) + ").";
				return false;
			}

			if ((int)countVirtual != nodeCount)
			{
				errorMessage = "Virtual displacement count at time " + std::to_string(virtStep.time) + " (" + std::to_string(countVirtual) + ") does not match mesh node count (" + std::to_string(nodeCount) + ").";
				return false;
			}

		}

		return true;
	}

	/**
	 * @brief Ensure measured load history has the same time coverage as measured displacements.
	 * @param data Optimisation container with displacement and load histories.
	 * @param errorMessage Populated when the validation fails.
	 * @return true when each measured displacement time has corresponding loads for all surfaces.
	 */
	static inline bool ValidateMeasuredLoads(const FEOptimizeDataVFM& data, std::string& errorMessage)
	{
		const auto& measuredHistory = data.MeasuredHistory();
		const auto& loadHistory = data.MeasuredLoads();

		if (measuredHistory.Empty())
		{
			errorMessage = "Measured displacement history is empty.";
			return false;
		}

		if (loadHistory.Empty())
		{
			errorMessage = "Measured load history is empty.";
			return false;
		}

		if (loadHistory.Steps() != measuredHistory.Steps())
		{
			errorMessage = "Measured load history contains " + std::to_string(loadHistory.Steps()) +
				" time steps; expected " + std::to_string(measuredHistory.Steps()) + " like measured displacements.";
			return false;
		}

		// Gather the union of surfaces defined across the load history.
		std::unordered_set<std::string> referenceSurfaces;
		for (const auto& loadStep : loadHistory)
		{
			for (const SurfaceLoadSample& sample : loadStep.loads.Samples())
			{
				if (sample.id.empty())
				{
					errorMessage = "Measured load history contains a surface entry with an empty id.";
					return false;
				}
				referenceSurfaces.insert(sample.id);
			}
		}

		if (referenceSurfaces.empty())
		{
			errorMessage = "Measured load history does not define any surfaces.";
			return false;
		}

		for (const auto& measStep : measuredHistory)
		{
			const auto* loadStep = loadHistory.FindStepByTime(measStep.time);
			if (loadStep == nullptr)
			{
				errorMessage = "Measured load history missing timestep for t = " + std::to_string(measStep.time) + ".";
				return false;
			}

			const double loadTime = loadStep->time;
			const SurfaceLoadSet& loads = loadStep->loads;
			if (loads.Size() != referenceSurfaces.size())
			{
				errorMessage = "Measured load history at t = " + std::to_string(loadTime) +
					" defines " + std::to_string(loads.Size()) + " surfaces; expected " +
					std::to_string(referenceSurfaces.size()) + ".";
				return false;
			}

			std::unordered_set<std::string> encountered;
			for (const SurfaceLoadSample& sample : loads.Samples())
			{
				if (referenceSurfaces.find(sample.id) == referenceSurfaces.end())
				{
					errorMessage = "Unexpected surface \"" + sample.id + "\" in measured loads at t = " + std::to_string(loadTime) + ".";
					return false;
				}
				encountered.insert(sample.id);
			}

			if (encountered.size() != referenceSurfaces.size())
			{
				errorMessage = "Duplicate or missing surface entries detected in measured loads at t = " + std::to_string(loadTime) + ".";
				return false;
			}
		}

		return true;
	}
};
