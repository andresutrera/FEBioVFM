#include <FEBioMech/stdafx.h>
#include <FECore/log.h>
#include <FEBioXML/xmltool.h>
#include <FECore/FECoreKernel.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include "FEData.h"
//=============================================================================
// FEVFMInput
//=============================================================================

/**
 * @brief Read the data from the XML optimization input file.
 *
 * The function validates the root tag (<febio_optimize version="2.0">) before
 * delegating parsing to finer grained helpers. The parser currently recognizes
 * <Parameters>, <MeasuredDisplacements>, <VirtualDisplacements>, and <MeasuredLoads>; unsupported
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
		else if (tag == "MeasuredLoads")
		{
			ParseMeasuredLoads(tag);
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
	history.Clear();

	auto parseNode = [](XMLTag& nodeTag, DisplacementContainer& container)
	{
		const char* szId = nodeTag.AttributeValue("id");
		int nodeId = szId ? atoi(szId) : -1;
		double disp[3] = { 0, 0, 0 };
		nodeTag.value(disp, 3);
		container.Add(nodeId, { disp[0], disp[1], disp[2] });
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

			XMLTag nodeTag = timeTag;
			++nodeTag;
			int nodeIndex = 0;
			while (!nodeTag.isend())
			{
				if ((nodeTag == "node") || (nodeTag == "elem"))
				{
					parseNode(nodeTag, step.displacements);
				}
				else
				{
					throw XMLReader::InvalidTag(nodeTag);
				}

				nodeTag.skip();
				++nodeTag;
			}
		}
		else
		{
		}

		timeTag.skip();
		++timeTag;
	}

	if (!found)
	{
		throw XMLReader::InvalidTag(tag);
	}


}

void ParseMeasuredLoadsBlock(XMLTag& tag, MeasuredLoadHistory& history)
{
	history.Clear();

	XMLTag timeTag(tag);
	++timeTag;
	bool found = false;
	while (!timeTag.isend())
	{
		if (timeTag == "time")
		{
			found = true;
			double timeValue = 0.0;
			timeTag.AttributeValue("t", timeValue, true);
			auto& step = history.AddStep(timeValue);

			XMLTag surfaceTag = timeTag;
			++surfaceTag;
			while (!surfaceTag.isend())
			{
				if (surfaceTag == "surface")
				{
					const char* szId = surfaceTag.AttributeValue("id");
					std::string surfaceId = (szId != nullptr) ? szId : "";

					double loadValues[3] = { 0, 0, 0 };
					surfaceTag.value(loadValues, 3);
					step.loads.Add(surfaceId, vec3d(loadValues[0], loadValues[1], loadValues[2]));
				}
				else
				{
					throw XMLReader::InvalidTag(surfaceTag);
				}

				surfaceTag.skip();
				++surfaceTag;
			}
		}

		timeTag.skip();
		++timeTag;
	}

	if (!found)
	{
		throw XMLReader::InvalidTag(tag);
	}
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
 	ParseDisplacementBlock(tag, m_opt->MeasuredHistory());
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
	m_opt->ClearVirtualFields();

	bool found = false;
	XMLTag fieldTag(tag);
	++fieldTag;
	while (!fieldTag.isend())
	{
			if (fieldTag == "virtualdisplacement")
			{
				found = true;
				const char* szId = fieldTag.AttributeValue("id", false);
				std::string fieldId = szId ? szId : "";
				auto& field = m_opt->AddVirtualField(fieldId);
				ParseDisplacementBlock(fieldTag, field.history);
			}
		else if (fieldTag == "time")
		{
				// Legacy format: times directly under <VirtualDisplacements>.
				if (!found)
				{
					auto& field = m_opt->AddVirtualField("");
				ParseDisplacementBlock(tag, field.history);
				found = true;
				break;
			}
			else
			{
				throw XMLReader::InvalidTag(fieldTag);
			}
		}
		else
		{
			throw XMLReader::InvalidTag(fieldTag);
		}

		fieldTag.skip();
		++fieldTag;
	}

	if (!found)
	{
		throw XMLReader::InvalidTag(tag);
	}
}

/**
 * @brief Parse measured surface loads from the XML input.
 *
 * @param tag XML tag positioned at <MeasuredLoads>.
 */
void FEVFMInput::ParseMeasuredLoads(XMLTag& tag)
{
	ParseMeasuredLoadsBlock(tag, m_opt->MeasuredLoads());
}
