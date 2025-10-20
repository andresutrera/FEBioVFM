/**
 * @file VFMKinematics.h
 * @brief Helpers for computing deformation gradients from measured displacements.
 */
#pragma once

#include <string>

class FEModel;
class DisplacementContainer;
class DeformationGradientField;

/**
 * @brief Computes deformation gradients without altering the FE model state.
 */
class VFMKinematics
{
public:
	/**
	 * @brief Populate deformation gradients at all Gauss points using the supplied displacements.
	 *
	 * @param fem FEBio model providing mesh and element information.
	 * @param displacements Nodal displacement definitions (measured or virtual).
	 * @param outField Destination container for computed gradients; cleared before filling.
	 * @param errorMessage Receives a human readable message when the computation fails.
	 * @return true on success, false otherwise.
	 */
	static bool ComputeDeformationGradients(FEModel& fem,
		const DisplacementContainer& displacements,
		DeformationGradientField& outField,
		std::string& errorMessage);
};
