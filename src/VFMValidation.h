/**
 * @file VFMValidation.h
 * @brief Validation helpers for FEBio models used by the VFM plugin.
 */
#pragma once

#include <string>
#include <algorithm>
#include <cmath>

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
			const double tm = measuredHistory.StepAt(i).time;
			const double tv = virtualHistory.StepAt(i).time;
			const double tf = data.DeformationHistory().StepAt(i).time;
			if (std::fabs(tm - tv) > 1e-12 || std::fabs(tm - tf) > 1e-12)
			{
				errorMessage = "Measured, virtual, and deformation gradient histories do not share the same time sequence.";
				return false;
			}
		}

		for (size_t i = 0; i < measuredHistory.Steps(); ++i)
		{
			const auto& measStep = measuredHistory.StepAt(i);
			const auto& virtStep = virtualHistory.StepAt(i);
			const auto& defStep = data.DeformationHistory().StepAt(i);

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

			// Deformation gradients are stored per element, so no node-count check here.
		}

		return true;
	}
};
