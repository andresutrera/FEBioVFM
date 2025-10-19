/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#pragma once
#include <XML/XMLReader.h>

class FEOptimizeDataVFM;

/**
 * @brief Reader for Virtual Fields Method optimization input files.
 *
	 * This helper extracts parameter definitions from the XML format used by the
	 * original FEBio optimization module and forwards the results to
	 * FEOptimizeDataVFM. It currently understands <Parameters>, <MeasuredDisplacements>,
	 * and <VirtualDisplacements>, mirroring the limited data required by the experimental
	 * scaffold.
 *
 * @note Additional sections (objective functions, constraints, etc.) can be
 * introduced later without changing the ownership model; the parser already
 * holds a pointer back to FEOptimizeDataVFM for dispatch.
 */
class FEVFMInput
{
public:
	/**
	 * @brief Read the XML optimization file and populate FEOptimizeDataVFM.
	 * @param szfile Path to the XML file passed on the FEBio command line.
	 * @param pOpt Destination container that receives parsed parameters.
	 * @return true on success, false when the file could not be parsed.
	 */
	bool Input(const char* szfile, FEOptimizeDataVFM* pOpt);

private:
		/**
		 * @brief Parse the <Parameters> subsection of the optimization input.
		 * @param tag XML tag positioned at the <Parameters> node.
		 *
		 * @note The routine instantiates FEModelParameterVFM objects dynamically and
		 * leaves ownership with FEOptimizeDataVFM, mimicking the behaviour of the
		 * legacy plugin for compatibility.
		 */
		void ParseParameters(XMLTag& tag);

		/**
		 * @brief Parse the <MeasuredDisplacements> subsection of the optimization input.
		 * @param tag XML tag positioned at the <MeasuredDisplacements> node.
		 *
		 * @note Each <node id="..."> entry (or legacy <elem id="...">) is expected to
		 * contain three comma- or space-separated floating-point values representing
		 * ux, uy, and uz.
		 */
		void ParseMeasuredDisplacements(XMLTag& tag);

		/**
		 * @brief Parse the <VirtualDisplacements> subsection of the optimization input.
		 * @param tag XML tag positioned at the <VirtualDisplacements> node.
		 *
		 * @note Virtual fields follow the same format as measured displacements and
		 * drive the inverse analysis portion of the VFM workflow.
		 */
			void ParseVirtualDisplacements(XMLTag& tag);

private:
			FEOptimizeDataVFM*	m_opt; ///< Destination optimization container populated during parsing.
	};
