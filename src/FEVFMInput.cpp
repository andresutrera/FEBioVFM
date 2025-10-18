#include <FEBioMech/stdafx.h>
#include <FECore/log.h>
#include <FEBioXML/xmltool.h>
#include <FECore/FECoreKernel.h>
#include "FEData.h"
//=============================================================================
// FEVFMInput
//=============================================================================

/**
 * @brief Read the data from the XML optimization input file.
 *
 * The function validates the root tag (<febio_optimize version="2.0">) before
 * delegating parsing to finer grained helpers. Only the <Parameters> section is
 * recognized today; unsupported tags trigger an exception that bubbles up to the
 * caller as a failure.
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
