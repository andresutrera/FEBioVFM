/**
 * @file VFMStress.h
 * @brief Helpers for evaluating stresses from precomputed deformation gradients.
 */
#pragma once

#include <string>

class FEModel;
class DeformationGradientField;
class StressField;
class FirstPiolaField;

/**
 * @brief Constitutive utilities that evaluate Gauss-point stresses without advancing the FE model.
 */
class VFMStress
{
public:
	/**
	 * @brief Compute Cauchy stresses for all Gauss points using cached deformation gradients.
	 *
	 * @param fem FEBio model that supplies mesh, material, and element definitions.
	 * @param defField Deformation gradients evaluated at each Gauss point.
	 * @param outField Destination container for the computed stresses; cleared before filling.
	 * @param errorMessage Populated with diagnostic text when the evaluation fails.
	 * @return true on success, false otherwise.
	 */
	static bool ComputeCauchyStress(FEModel& fem,
		const DeformationGradientField& defField,
		StressField& outField,
		std::string& errorMessage);

	/**
	 * @brief Transform Cauchy stresses into first Piola-Kirchhoff stresses.
	 *
	 * @param defField Deformation gradients evaluated at each Gauss point.
	 * @param cauchyField Cauchy stresses previously computed for the same Gauss points.
	 * @param outField Destination container for the transformed stresses; cleared before filling.
	 * @param errorMessage Populated with diagnostic text when the evaluation fails.
	 * @return true on success, false otherwise.
	 */
	static bool ComputeFirstPiolaStress(const DeformationGradientField& defField,
		const StressField& cauchyField,
		FirstPiolaField& outField,
		std::string& errorMessage);
};
