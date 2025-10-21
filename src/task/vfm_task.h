// src/task/fevfm_task.hpp
#pragma once
#include <string>
#include <FECore/FECoreTask.h>

#include "optimization/vfm_problem.hpp"
#include "io/xml_reader.hpp"   // VFMXmlReader, XMLInput

class VFMTask final : public FECoreTask {
public:
    explicit VFMTask(FEModel* fem);

    bool Init(const char* xmlPath) override; // parse XML only
    bool Run() override;                     // no-op for now

    const XMLInput& input() const { return m_input; }

private:
    std::string m_inputPath;
    XMLInput    m_input;     // plain DTO parsed from XML
    VFMProblem m_problem;
};
