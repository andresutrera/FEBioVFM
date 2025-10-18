/**
 * @file VFMValidation.h
 * @brief Validation helpers for FEBio models used by the VFM plugin.
 */
#pragma once

#include <string>

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

		const size_t measuredCount = data.MeasuredData().Size();
		const size_t virtualCount = data.VirtualData().Size();

		if ((int)measuredCount != nodeCount)
		{
			errorMessage = "Measured displacement count (" + std::to_string(measuredCount) + ") does not match mesh node count (" + std::to_string(nodeCount) + ").";
			return false;
		}

		if ((int)virtualCount != nodeCount)
		{
			errorMessage = "Virtual displacement count (" + std::to_string(virtualCount) + ") does not match mesh node count (" + std::to_string(nodeCount) + ").";
			return false;
		}

		return true;
	}
};
