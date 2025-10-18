/**
 * @file VFMExport.h
 * @brief Helpers for exporting VFM kinematic quantities to FEBio plot files.
 */
#pragma once

#include <string>

class FEModel;
class DisplacementContainer;
class DeformationGradientField;

/**
 * @brief Export measured/virtual displacements and deformation gradients to an XPLT file.
 *
 * @param filePath Destination path for the generated plot file.
 * @param fem      Reference to the FEModel that owns the mesh being exported.
 * @param measured Container with measured nodal displacements.
 * @param virt     Container with virtual nodal displacements.
 * @param defGrad  Gauss-point deformation gradients reconstructed by the VFM kinematics pass.
 * @param error    Filled with a descriptive message when the export fails.
 * @return true on success, false otherwise.
 */
bool ExportVFMKinematics(const std::string& filePath,
	FEModel& fem,
	const DisplacementContainer& measured,
	const DisplacementContainer& virt,
	const DeformationGradientField& defGrad,
	std::string& error);
