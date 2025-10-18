#include <FEBioMech/stdafx.h>
#include <FECore/log.h>
#include <FEBioXML/xmltool.h>
#include <FECore/FECoreKernel.h>
#include <cstdlib>
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
		do
		{
			if (tag == "Parameters" ) ParseParameters(tag);
			else if (tag == "MeasuredDisplacements") ParseMeasuredDisplacements(tag);
			else if (tag == "VirtualDisplacements") ParseVirtualDisplacements(tag);
			else throw XMLReader::InvalidTag(tag);
			++tag;
		} while (!tag.isend());
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
void ParseDisplacementBlock(XMLTag& tag, DisplacementContainer& container)
{
	container.Clear();

	++tag;
	do
	{
		if (tag == "elem")
		{
			const char* szId = tag.AttributeValue("id");
			int elemId = szId ? atoi(szId) : -1;
			double disp[3] = { 0, 0, 0 };
			tag.value(disp, 3);

			container.Add(elemId, { disp[0], disp[1], disp[2] });
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	} while (!tag.isend());
}

} // namespace

/**
 * @brief Parse experimental displacement measurements from the XML input.
 *
 * @param tag XML tag positioned at <MeasuredDisplacements>.
 *
 * Each child <elem id="..."> tag must provide three floating-point numbers that
 * correspond to ux, uy, and uz. Values may be comma-separated or whitespace
 * separated per FEBio's XML parser conventions. The container is cleared before
 * new samples are appended so repeated definitions overwrite previous entries.
 */
void FEVFMInput::ParseMeasuredDisplacements(XMLTag& tag)
{
	ParseDisplacementBlock(tag, m_opt->MeasuredData());
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
	ParseDisplacementBlock(tag, m_opt->VirtualData());
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

	const auto& measuredSamples = m_opt->MeasuredData().Samples();
	if (measuredSamples.empty())
	{
		feLogDebugEx(fem, "  Measured displacements: (none)");
	}
	else
	{
		feLogDebugEx(fem, "  Measured displacements (%zu entries)", measuredSamples.size());
		for (const ElementDisplacement& entry : measuredSamples)
		{
			feLogDebugEx(fem, "    elem %6d : ux=%-12g uy=%-12g uz=%-12g", entry.id, entry.displacement[0], entry.displacement[1], entry.displacement[2]);
		}
	}

	const auto& virtualSamples = m_opt->VirtualData().Samples();
	if (virtualSamples.empty())
	{
		feLogDebugEx(fem, "  Virtual displacements: (none)");
	}
	else
	{
		feLogDebugEx(fem, "  Virtual displacements (%zu entries)", virtualSamples.size());
		for (const ElementDisplacement& entry : virtualSamples)
		{
			feLogDebugEx(fem, "    elem %6d : ux=%-12g uy=%-12g uz=%-12g", entry.id, entry.displacement[0], entry.displacement[1], entry.displacement[2]);
		}
	}
}
