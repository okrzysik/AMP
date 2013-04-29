#ifndef included_AMP_TensorProperty
#define included_AMP_TensorProperty

#include "Property.h"

namespace AMP
{
namespace Materials
{

template<class Number>
class TensorProperty: public Property<Number>
{
public:
	/**
	 * Constructor
	 * \param name name of property
	 * \param source literature reference for model and notes
	 * \param params default parameter values
	 * \param nparams number of parameter values
	 * \param args names of arguments
	 * \param nargs number of arguments
	 * \param ranges ranges of arguments
	 * \param dimensions the first and second dimensions of return value tensor
	 */
	TensorProperty(
			const std::string& name = std::string("NotDefined"),
			const std::string& source = std::string("None"),
			const Number *params = NULL,
			const unsigned int nparams = 0,
			const std::string* args = NULL,
			const unsigned int nargs = 0,
			const Number ranges[][2] = NULL,
			const std::vector<size_t>& dimensions = std::vector<size_t>(2,1)) :
				Property<Number>(name,source,params,nparams,args,nargs,ranges),
				d_dimensions(dimensions), d_variableDimensions(false)
	{
		AMP_INSIST(d_dimensions.size()==2, "there must be two dimensions");
		AMP_INSIST(d_dimensions[0]>0, "must have first return tensor dimension > 0");
		AMP_INSIST(d_dimensions[1]>0, "must have second return tensor dimension > 0");
	}

	/**
	 * Destructor
	 */
	virtual ~TensorProperty()
	{
	}

	/// get dimensions of evalv return value tensor
	std::vector<size_t> get_dimensions(){return d_dimensions;}

	/// ok to change dimensions
	bool variable_dimensions(){return d_variableDimensions;}

	/// set dimensions of evalv return value tensor
	void set_dimensions(std::vector<size_t> dimensions)
	{
		AMP_INSIST(d_variableDimensions, "can not change dimensions for this property");
		d_dimensions = dimensions;
		AMP_INSIST(d_dimensions.size()==2, "there must be two dimensions");
		AMP_INSIST(d_dimensions[0]>0, "must have first return tensor dimension > 0");
		AMP_INSIST(d_dimensions[1]>0, "must have second return tensor dimension > 0");
	}

	/// indicator for scalar evaluator
	virtual bool isScalar(){return false;}

	/// indicator for tensor evaluator
	virtual bool isTensor(){return true;}

protected:
	std::vector<size_t> d_dimensions; ///< dimensions of return value tensor
	bool d_variableDimensions; ///< true if ok to change dimensions

///////////////////// Evaluators /////////////////////

private:
	/* Loops through input vectors, calling the child eval function, returning tensor results */
	template<class INPUT_VTYPE, class RETURN_VTYPE>
	void
	evalvActual(
		std::vector<std::vector<boost::shared_ptr<RETURN_VTYPE> > >& r,
		const std::map<std::string, boost::shared_ptr<INPUT_VTYPE> >& args);
public:

	/**
	 * tensor evaluation function for a single argument set
	 * \param args list of argument values, in correct order, given by  get_arguments()
	 * \return tensor of property values with dimensions get_dimensions()
	 */
	virtual std::vector<std::vector<Number> >
	evalTensor(std::vector<Number>& args)=0;

	/** Wrapper function that calls evalvActual for each argument set
	 *  \param r tensor vector of return values
	 *  \param args map of vectors of arguments, indexed by strings which are members of  get_arguments()
	 *
	 *  The  \a r parameter must have dimensions get_dimensions().
	 *  The  \a args  parameter need not contain all the members of  get_arguments()  as indices.
	 *  Arguments left out will have values supplied by the entries in  get_defaults() .
	 *  Sizes of  *r[i][j]  and \a args["name"] must match. Members of
	 *  \a args  indexed by names other than those in  get_arguments()  are ignored.
	 *  The list {args["name-1"][i], ..., args["name-n"][i]} will be passed to eval() and the j-th result
	 *  returned in (*r[j])[i].
	 */
	virtual void
	evalv(
		std::vector<std::vector<boost::shared_ptr<std::vector<Number> > > >& r,
		const std::map<std::string, boost::shared_ptr<std::vector<Number> > >& args);

	/** Wrapper function that calls evalvActual for each argument set
	 *  \param r tensor of AMP vectors of return values
	 *  \param args map of AMP vectors of arguments, indexed by strings which are members of  get_arguments()
	 *
	 *  The  \a args  parameter need not contain all the members of  get_arguments()  as indices.
	 *  Arguments left out will have values supplied by the entries in  get_defaults() .
	 *  Sizes of  \a *r[i][j]  and \a args["name"] must match. Members of
	 *  \a args  indexed by names other than those in  get_arguments()  are ignored.
	 *  The list {args["name-1"][i], ..., args["name-n"][i]} will be passed to eval() and the k-j-th result
	 *  returned in (*r[k][j])[i].
	 */
	virtual void
	evalv(
		std::vector<std::vector<boost::shared_ptr<AMP::LinearAlgebra::Vector> > >& r,
		const std::map<std::string, boost::shared_ptr<AMP::LinearAlgebra::Vector> >& args);

	/** Wrapper function that calls evalvActualVector for each argument set
	 *  \param r tensor of AMP vectors of return values
	 *  \param args AMP multivector of arguments
	 *
	 *  Before this function is used, a translation table must be assigned by means of the set_translator() function
	 *  which gives the correspondence between entries in get_arguments() and the \a args multivector.
	 *  Upon invocation, the \a args parameter is converted to a map of AMP vectors via make_map() and passed to another version of evalv.
	 */
	virtual void
	evalv(
		std::vector<std::vector<boost::shared_ptr<AMP::LinearAlgebra::Vector> > >& r,
		const boost::shared_ptr<AMP::LinearAlgebra::MultiVector>& args);

	// disable scalar evaluator
	virtual Number
	eval(std::vector<Number>& args)
	{AMP_INSIST(false, "cannot use scalar evaluator from tensor property"); return 0;}

	// disable scalar evaluator
	virtual void
	evalv(std::vector<Number>& r,
		const std::map<std::string, boost::shared_ptr<std::vector<Number> > >& args)
	{AMP_INSIST(false, "cannot use scalar evaluator from tensor property");}

	// disable scalar evaluator
	virtual void
	evalv(boost::shared_ptr<AMP::LinearAlgebra::Vector>& r,
		const std::map<std::string, boost::shared_ptr<AMP::LinearAlgebra::Vector> >& args)
	{AMP_INSIST(false, "cannot use scalar evaluator from tensor property");}

	// disable scalar evaluator
	virtual void
	evalv(boost::shared_ptr<AMP::LinearAlgebra::Vector>& r,
		const boost::shared_ptr<AMP::LinearAlgebra::MultiVector>& args)
	{AMP_INSIST(false, "cannot use scalar evaluator from tensor property");}
};

template<>
void
TensorProperty<double>::evalv(
	std::vector<std::vector<boost::shared_ptr<AMP::LinearAlgebra::Vector> > >& r,
	const std::map<std::string, boost::shared_ptr<AMP::LinearAlgebra::Vector> >& args);

template<>
void
TensorProperty<double>::evalv(
	std::vector<std::vector<boost::shared_ptr<AMP::LinearAlgebra::Vector> > >& r,
	const boost::shared_ptr<AMP::LinearAlgebra::MultiVector>& args);

} // namespace Materials
} // namespace AMP

#include "TensorProperty.i.h"

#endif /* included_AMP_TensorProperty */
