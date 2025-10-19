/**
 * @file VFMExport.h
 * @brief Helpers for exporting VFM kinematic quantities to FEBio plot files.
 */
#pragma once

#include <string>

#include "VirtualDisplacementContainer.h"
#include "VirtualDeformationGradientContainer.h"

class FEModel;
class DisplacementHistory;
class DeformationGradientHistory;
class StressHistory;

/**
 * @brief Export measured/virtual displacements and deformation gradients to an XPLT file.
 *
 * @param filePath Destination path for the generated plot file.
 * @param fem      Reference to the FEModel that owns the mesh being exported.
 * @param measuredHist Timeline of measured nodal displacements.
 * @param virtualHist  Timeline of virtual nodal displacements.
 * @param defHist      Timeline of Gauss-point deformation gradients reconstructed by the VFM kinematics pass.
 * @param stressHist   Timeline of Gauss-point stresses reconstructed from deformation gradients.
 * @param error    Filled with a descriptive message when the export fails.
 * @return true on success, false otherwise.
 */
bool ExportVFMKinematics(const std::string& filePath,
	FEModel& fem,
	const DisplacementHistory& measuredHist,
	const VirtualDisplacementCollection& virtualFields,
	const VirtualDeformationGradientCollection& virtualGradients,
	const DeformationGradientHistory& defHist,
	const StressHistory& stressHist,
	std::string& error);
