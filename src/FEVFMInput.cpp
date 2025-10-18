#include <FEBioMech/stdafx.h>
#include <FECore/log.h>
#include <FEBioXML/xmltool.h>
#include <FECore/FECoreKernel.h>
#include <cstdlib>
#include <iostream>
#include "FEData.h"
//=============================================================================
// FEVFMInput
//=============================================================================

/**
 * @brief Read the data from the XML optimization input file.
 *
 * The function validates the root tag (<febio_optimize version="2.0">) before
 * delegating parsing to finer grained helpers. The parser currently recognizes
 * <Parameters>, <MeasuredDisplacements>, and <VirtualDisplacements>; unsupported
 * tags trigger an exception that bubbles up to the caller as a failure.
 *
 * @note Error handling mirrors the classic FEBio loader by printing messages to
 * stderr so the user receives immediate feedback even when the FEBio log is muted.
 */
bool FEVFMInput::Input(const char* szfile, FEOptimizeDataVFM* pOpt)
{
	// try to open the file
	XMLReader xml;
	if (xml.Open(szfile) == false)
	{
		fprintf(stderr, "\nFATAL ERROR: Failed to load file %s\n", szfile);
		return false;
	}

	// find the root element
	XMLTag tag;
	if (xml.FindTag("febio_optimize", tag) == false) return false;

	// check for the version attribute
	const char* szversion = tag.AttributeValue("version", true);
	if ((szversion == 0) || (strcmp(szversion, "2.0") != 0))
	{
		fprintf(stderr, "\nFATAL ERROR: Invalid version number for febio_optimize!\n\n");
		return false;
	}

	m_opt = pOpt;

	// process the file
	bool ret = true;
	try {
	++tag;
	while (!tag.isend())
	{
		std::cout << "Processing tag: " << tag.Name() << std::endl;

		if (tag.Name()[0] == '#')
		{
			tag.skip();
			++tag;
			continue;
		}

		if (tag == "Parameters")
		{
			ParseParameters(tag);
			tag.skip();
		}
		else if (tag == "MeasuredDisplacements")
		{
			ParseMeasuredDisplacements(tag);
			tag.skip();
		}
		else if (tag == "VirtualDisplacements")
		{
			ParseVirtualDisplacements(tag);
			tag.skip();
		}
		else
		{
			throw XMLReader::InvalidTag(tag);
		}

		++tag;
	}
	}
	catch (...)
	{
		fprintf(stderr, "Fatal exception while reading optimization input file.\n");
		ret = false;
	}

	// all done
	xml.Close();

	if (ret) LogDebugSummary();

	return ret;
}


/**
 * @brief Parse the <Parameters> block of the optimization input file.
 *
 * The implementation instantiates FEModelParameterVFM objects for each
 * <param name="..."> entry and copies the numeric tuple (initial value, lower
 * bound, upper bound, scale). Unrecognized tags result in a thrown exception to
 * make configuration errors visible early.
 *
 * @note Only scalar double parameters are supported at the moment because the
 * downstream FEModelParameterVFM class is limited to FE_PARAM_DOUBLE.
 */
void FEVFMInput::ParseParameters(XMLTag& tag)
{
	FEModel& fem = *m_opt->GetFEModel();

	// read the parameters
	++tag;
	do
	{
		if (tag == "param")
		{
			FEModelParameterVFM* var = new FEModelParameterVFM(&fem);

			// get the variable name
			const char* sz = tag.AttributeValue("name");
			var->SetName(sz);

			// set initial values and bounds
			double d[4] = { 0, 0, 0, 1 };
			tag.value(d, 4);
			var->InitValue() = d[0];
			var->MinValue() = d[1];
			var->MaxValue() = d[2];
			var->ScaleFactor() = d[3];

			// add the variable
			m_opt->AddInputParameter(var);
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	} while (!tag.isend());
}

namespace {

// Helper that parses either measured or virtual displacement blocks.
void ParseDisplacementBlock(XMLTag& tag, DisplacementHistory& history)
{
	std::cout << "VFM parse: clearing displacement history" << std::endl;
	history.Clear();

	auto parseNode = [](XMLTag& nodeTag, DisplacementContainer& container)
	{
		const char* szId = nodeTag.AttributeValue("id");
		int nodeId = szId ? atoi(szId) : -1;
		double disp[3] = { 0, 0, 0 };
		nodeTag.value(disp, 3);
		container.Add(nodeId, { disp[0], disp[1], disp[2] });
		std::cout << "  node " << nodeId << " -> (" << disp[0] << ", " << disp[1] << ", " << disp[2] << ")" << std::endl;
	};

	XMLTag timeTag(tag);
	++timeTag;
	bool found = false;
	int timeIndex = 0;
	while (!timeTag.isend())
	{
		if (timeTag == "time")
		{
			found = true;
			double timeValue = 0.0;
			timeTag.AttributeValue("t", timeValue, true);
			auto& step = history.AddStep(timeValue);
			std::cout << "VFM parse: created time step " << timeIndex++ << " at t = " << timeValue << std::endl;

			XMLTag nodeTag = timeTag;
			++nodeTag;
			int nodeIndex = 0;
			while (!nodeTag.isend())
			{
				if ((nodeTag == "node") || (nodeTag == "elem"))
				{
					std::cout << "  parsing node index " << nodeIndex++ << std::endl;
					parseNode(nodeTag, step.displacements);
				}
				else
				{
					std::cout << "  encountered unexpected tag <" << nodeTag.Name() << ">" << std::endl;
					throw XMLReader::InvalidTag(nodeTag);
				}

				nodeTag.skip();
				++nodeTag;
			}
		}
		else
		{
			std::cout << "VFM parse: skipping non-time tag <" << timeTag.Name() << ">" << std::endl;
		}

		timeTag.skip();
		++timeTag;
	}

	if (!found)
	{
		std::cout << "VFM parse: no <time> entries detected, aborting" << std::endl;
		throw XMLReader::InvalidTag(tag);
	}

	history.SetActiveStepByIndex(0);
	std::cout << "VFM parse: active step set to index 0 (t = " << history.ActiveStep().time << ")" << std::endl;
}

} // namespace

/**
 * @brief Parse experimental displacement measurements from the XML input.
 *
 * @param tag XML tag positioned at <MeasuredDisplacements>.
 *
 * Each child <node id="..."> (or legacy <elem id="...">) tag must provide three floating-point numbers that
 * correspond to ux, uy, and uz. Values may be comma-separated or whitespace
 * separated per FEBio's XML parser conventions. The container is cleared before
 * new samples are appended so repeated definitions overwrite previous entries.
 */
void FEVFMInput::ParseMeasuredDisplacements(XMLTag& tag)
{
    std::cout << "ParseMeasuredDisplacements begin" << std::endl;
 	ParseDisplacementBlock(tag, m_opt->MeasuredHistory());
    std::cout << "ParseMeasuredDisplacements end" << std::endl;
}

/**
 * @brief Parse user-defined virtual displacement fields from the XML input.
 *
 * @param tag XML tag positioned at <VirtualDisplacements>.
 *
 * Virtual displacements share the same syntax as measured data but correspond to
 * the admissible virtual fields required by the Virtual Fields Method.
 */
void FEVFMInput::ParseVirtualDisplacements(XMLTag& tag)
{
    std::cout << "ParseVirtualDisplacements begin" << std::endl;
	ParseDisplacementBlock(tag, m_opt->VirtualHistory());
    std::cout << "ParseVirtualDisplacements end" << std::endl;
}

/**
 * @brief Emit a formatted summary of the parsed configuration using debug logging.
 */
void FEVFMInput::LogDebugSummary() const
{
	if (m_opt == nullptr) return;

	FEModel* fem = m_opt->GetFEModel();
	if (fem == nullptr) return;

	feLogDebugEx(fem, "---- VFM Input Summary --------------------------------");

	const int paramCount = m_opt->InputParameters();
	feLogDebugEx(fem, "  Parameters to optimise: %d", paramCount);
	for (int i = 0; i < paramCount; ++i)
	{
		FEInputParameterVFM* param = m_opt->GetInputParameter(i);
		if (param == nullptr) continue;

		const std::string name = param->GetName();
		const double initVal = param->InitValue();
		const double minVal = param->MinValue();
		const double maxVal = param->MaxValue();

		feLogDebugEx(fem, "    %-20s init=%-12g min=%-12g max=%-12g", name.c_str(), initVal, minVal, maxVal);
	}

	const auto& measuredHistory = m_opt->MeasuredHistory();
	feLogDebugEx(fem, "  Measured displacement steps: %zu", measuredHistory.Steps());
	for (const auto& step : measuredHistory.StepsRef())
	{
		feLogDebugEx(fem, "    t = %-12g (%zu entries)", step.time, step.displacements.Size());
		for (const NodeDisplacement& entry : step.displacements.Samples())
		{
			feLogDebugEx(fem, "      node %6d : ux=%-12g uy=%-12g uz=%-12g", entry.id, entry.displacement[0], entry.displacement[1], entry.displacement[2]);
		}
	}

	const auto& virtualHistory = m_opt->VirtualHistory();
	feLogDebugEx(fem, "  Virtual displacement steps: %zu", virtualHistory.Steps());
	for (const auto& step : virtualHistory.StepsRef())
	{
		feLogDebugEx(fem, "    t = %-12g (%zu entries)", step.time, step.displacements.Size());
		for (const NodeDisplacement& entry : step.displacements.Samples())
		{
			feLogDebugEx(fem, "      node %6d : ux=%-12g uy=%-12g uz=%-12g", entry.id, entry.displacement[0], entry.displacement[1], entry.displacement[2]);
		}
	}
}
