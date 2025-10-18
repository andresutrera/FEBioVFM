#pragma once
#include <FECore/FEModel.h>
#include <vector>
#include <string>
#include "FEVFMInput.h"
#include "DisplacementContainer.h"
#include "DeformationGradientField.h"

/**
 * @brief Abstract base for a scalar optimization variable used by the VFM plugin.
 *
 * Derived classes bridge between the optimizer and FEBio by presenting a scalar
 * value that can be queried and overwritten, while internally mapping those reads
 * and writes to model state. Metadata such as bounds, scaling, and an initial guess
 * are cached here so Doxygen exposes the parameterization in the generated docs.
 *
 * @note The class only stores the FEModel pointer; it never takes ownership and
 * assumes the model outlives the parameter representation.
 */
class FEInputParameterVFM
{
public:
	/**
	 * @brief Construct an input parameter bound to a FEBio model context.
	 * @param fem Owning FEBio model that provides parameter metadata and storage.
	 *
	 * @note Lookup of concrete model data is deferred to derived classes because
	 * some parameters are not available until after the input file has been parsed.
	 */
	FEInputParameterVFM(FEModel* fem) : m_fem(fem) { m_min = -1e99; m_max = 1e99; m_scale = 1.0; }
	virtual ~FEInputParameterVFM() {}

	/**
	 * @brief Optional hook that gives derived classes a chance to prepare state.
	 * @return true if the parameter is ready for use by the optimizer.
	 *
	 * @note The base implementation simply returns true; override in subclasses
	 * when the parameter depends on FEBio structures that need to be resolved.
	 */
	virtual bool Init() { return true; }

	/**
	 * @brief Retrieve the current parameter value from the underlying data source.
	 */
	virtual double GetValue() = 0;

	/**
	 * @brief Push a new value to the underlying data source.
	 * @param newValue Proposed value in physical units.
	 * @return true if the update was accepted and applied.
	 *
	 * @note Implementations should reject invalid values (e.g. outside allowable
	 * bounds) by returning false and leaving the underlying model unchanged.
	 */
	virtual bool SetValue(double newValue) = 0;

	/**
	 * @brief Get or set the initial value used to seed the optimization run.
	 */
	double& InitValue() { return m_initVal; }

	/**
	 * @brief Get or set the lower admissible bound for this variable.
	 */
	double& MinValue() { return m_min; }

	/**
	 * @brief Get or set the upper admissible bound for this variable.
	 */
	double& MaxValue() { return m_max; }

	/**
	 * @brief Get or set the scaling factor that normalizes sensitivity during search.
	 *
	 * @note Scaling is currently only stored; no automatic normalization is applied
	 * elsewhere in the code base yet.
	 */
	double& ScaleFactor() { return m_scale; }

	/**
	 * @brief Set the human readable name of the parameter.
	 * @param name Identifier used to resolve the underlying FEBio parameter.
	 */
	void SetName(const string& name) { m_name = name; }

	/**
	 * @brief Retrieve the human readable name of the parameter.
	 */
	string GetName() { return m_name; }

	/**
	 * @brief Access the FEBio model associated with the parameter.
	 */
	FEModel* GetFEModel() { return m_fem; }

private:
	string		m_name;			//!< name of input parameter
	double		m_initVal;		//!< initial value
	double		m_min, m_max;	//!< min, max values for parameter
	double		m_scale;		//!< scale factor
	FEModel*	m_fem;			//!< pointer to model data
};
/**
 * @brief Adapter that exposes a FEBio model parameter as an optimization variable.
 *
 * The class resolves an FEModel parameter by name, caches a pointer to the raw
 * value, and mirrors optimizer updates back to the model. Doxygen consumers can
 * see that this class is meant for direct parameter tuning rather than custom
 * derived quantities.
 *
 * @note The adapter expects the target parameter to be of type FE_PARAM_DOUBLE.
 * Other types will trigger a descriptive log error during initialization.
 */
class FEModelParameterVFM : public FEInputParameterVFM
{
public:
	/**
	 * @brief Construct an adapter bound to the supplied FEBio model.
	 * @param fem Owning FEBio model.
	 */
	FEModelParameterVFM(FEModel* fem);

	/**
	 * @brief Locate the targeted model parameter and cache a pointer to it.
	 * @return true when the parameter exists and is accessible.
	 *
	 * @note Initialization must be called after the parent FEModel parsed its
	 * input file, otherwise the lookup performed here will fail.
	 */
	bool Init();

	/**
	 * @brief Return the current value of the model parameter.
	 *
	 * @note Falls back to zero when initialization never resolved the pointer,
	 * mirroring the defensive behaviour in the original plugin.
	 */
	double GetValue();

	/**
	 * @brief Write a new value into the underlying FEModel parameter.
	 * @param newValue Proposed value in physical units.
	 * @return true if the pointer was valid and the assignment succeeded.
	 *
	 * @note The current implementation rejects updates when @c Init() did not
	 * resolve a pointer, allowing the optimizer to detect misconfigured inputs.
	 */
	bool SetValue(double newValue);

private:
	string	m_name;		//!< name of the FEModel parameter this adapter targets
	double*	m_pd;		//!< pointer to the raw model data once resolved in Init
	double	m_val;		//!< cached value used when the model data pointer is absent
};


/**
 * @brief Container for optimization state managed by the VFM plugin.
 *
 * The class owns the collection of input parameters, coordinates their
 * initialization, and is responsible for driving the eventual forward solve
 * and inverse update loop. Several methods remain placeholders to highlight
 * the intended flow while more of the algorithm is being implemented.
 */
class FEOptimizeDataVFM
{
public:
	/**
	 * @brief Construct the optimization state wrapper.
	 * @param fem FEBio model that will be solved during optimization.
	 */
	FEOptimizeDataVFM(FEModel* fem);
	~FEOptimizeDataVFM(void);

	/**
	 * @brief Parse the VFM-specific input file.
	 * @param sz Path to the XML file that describes optimization parameters.
	 */
	bool Input(const char* sz);

	/**
	 * @brief Initialize the optimization problem before the first solve.
	 *
	 * @note This routine currently disables FEBio plot output and seeds each
	 * registered parameter with its initial value. Solver and objective wiring
	 * are left as future work and logged in comments for traceability.
	 */
	bool Init();

	/**
	 * @brief Execute the optimization loop.
	 *
	 * @return Always returns false until the loop is implemented.
	 * @note The function is not implemented yet; the body currently returns
	 * false so callers and generated documentation surface the missing feature.
	 */
	bool Solve();

	/**
	 * @brief Access the FEModel instance used during optimization.
	 */
	FEModel* GetFEModel() { return m_fem; }

	/**
	 * @brief Mutable access to the measured displacement data set.
	 */
	DisplacementContainer& MeasuredData() { return m_measured; }

	/**
	 * @brief Read-only access to the measured displacement data set.
	 */
	const DisplacementContainer& MeasuredData() const { return m_measured; }

	/**
	 * @brief Mutable access to the prescribed virtual displacement fields.
	 */
	DisplacementContainer& VirtualData() { return m_virtual; }

	/**
	 * @brief Read-only access to the prescribed virtual displacement fields.
	 */
	const DisplacementContainer& VirtualData() const { return m_virtual; }

	/**
	 * @brief Mutable access to the reconstructed deformation gradients.
	 */
	DeformationGradientField& DeformationGradients() { return m_defGrad; }

	/**
	 * @brief Read-only access to the reconstructed deformation gradients.
	 */
	const DeformationGradientField& DeformationGradients() const { return m_defGrad; }

	/**
	 * @brief Solve the forward FE problem with a proposed parameter vector.
	 * @param a Ordered list of parameter values sourced from the optimizer.
	 *
	 * @return Always returns false until the solver hook is implemented.
	 * @note This is currently a stub that will eventually call into FEBio's
	 * task infrastructure once the inversion loop is assembled.
	 */
	bool FESolve(const std::vector<double>& a);

public:
	/**
	 * @brief Return the number of registered input parameters.
	 */
	int InputParameters() { return (int)m_Var.size(); }

	/**
	 * @brief Register a new input parameter with the optimization problem.
	 * @param var Parameter instance; ownership stays with the caller.
	 *
	 * @note Parameters are stored as raw pointers to mirror the original plugin.
	 * The container does not delete them, so lifetime must be managed externally.
	 */
	void AddInputParameter(FEInputParameterVFM* var) { m_Var.push_back(var); }

	/**
	 * @brief Retrieve a parameter by position.
	 */
	FEInputParameterVFM* GetInputParameter(int n) { return m_Var[n]; }

public:



public:
	int	m_niter;	///< Number of minor iterations (i.e. FE solves) executed so far.


protected:
	FEModel*	m_fem;   ///< FEBio model associated with the optimization run.
	std::vector<FEInputParameterVFM*>	    m_Var; ///< Registered optimization variables (non-owning).
	DisplacementContainer m_measured; ///< Experimentally measured displacements.
	DisplacementContainer m_virtual;  ///< User-specified virtual displacement fields.
	DeformationGradientField m_defGrad; ///< Cached deformation gradients per Gauss point.
};
