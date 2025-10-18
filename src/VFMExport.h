/**
 * @file VFMExport.h
 * @brief Helpers for exporting VFM kinematic quantities to FEBio plot files.
 */
#pragma once

#include <string>

class FEModel;
class DisplacementHistory;
class DeformationGradientHistory;

/**
 * @brief Export measured/virtual displacements and deformation gradients to an XPLT file.
 *
 * @param filePath Destination path for the generated plot file.
 * @param fem      Reference to the FEModel that owns the mesh being exported.
 * @param measuredHist Timeline of measured nodal displacements.
 * @param virtualHist  Timeline of virtual nodal displacements.
 * @param defHist      Timeline of Gauss-point deformation gradients reconstructed by the VFM kinematics pass.
 * @param error    Filled with a descriptive message when the export fails.
 * @return true on success, false otherwise.
 */
bool ExportVFMKinematics(const std::string& filePath,
	FEModel& fem,
	const DisplacementHistory& measuredHist,
	const DisplacementHistory& virtualHist,
	const DeformationGradientHistory& defHist,
	std::string& error);
