/**
 * @file VFMExport.h
 * @brief Helpers for exporting VFM kinematic quantities to FEBio plot files.
 */
#pragma once

#include <memory>
#include <string>

#include "VirtualDisplacementContainer.h"
#include "VirtualDeformationGradientContainer.h"

class FEModel;
class DisplacementHistory;
class DeformationGradientHistory;
class StressHistory;

/**
 * @brief Helper class that stages VFM export data and writes it to an XPLT file.
 *
 * The session registers plot variables lazily as the caller adds data sets and
 * postpones the actual file write until @c Finalize is invoked.
 */
class VFMExportSession
{
public:
	VFMExportSession(const std::string& filePath, FEModel& fem);
	~VFMExportSession();

	bool AddMeasuredDisplacements(const DisplacementHistory& hist, std::string& error);
	bool AddVirtualDisplacements(const VirtualDisplacementCollection& fields, std::string& error);
	bool AddVirtualDeformationGradients(const VirtualDeformationGradientCollection& fields, std::string& error);
	bool AddMeasuredDeformationGradients(const DeformationGradientHistory& hist, std::string& error);
	bool AddMeasuredStress(const StressHistory& hist, std::string& error);
	bool Finalize(std::string& error);

private:
	struct Impl;
	std::unique_ptr<Impl> m_impl;
};
