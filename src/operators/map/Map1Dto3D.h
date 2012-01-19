#ifndef included_AMP_Map1Dto3D
#define included_AMP_Map1Dto3D


#include "boost/shared_ptr.hpp"
#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "operators/map/MapOperator.h"
#include "operators/map/MapOperatorParameters.h"
#include "vectors/Vector.h"
#include "vectors/Variable.h"

#include <string>

namespace AMP {
namespace Operator {

  /**
    Class Map1Dto3D is used to map the 1D solution from simple vector 
    onto the 2D Surface by injunction as the 1D location are constructed 
    by the nodal points on the 2D surface. The 2D surface has to be on the
    d_MapAdapater from the MapOperator class. 
    */
    class Map1Dto3D : public MapOperator
	{

		public :

			Map1Dto3D(const boost::shared_ptr<OperatorParameters>& params)
				: MapOperator (params)
			{
				boost::shared_ptr<MapOperatorParameters> myparams = 
					boost::dynamic_pointer_cast<MapOperatorParameters>(params);

				d_MapMesh = myparams->d_MapMesh;
				reset(myparams);
			}

			virtual ~Map1Dto3D() { }

			/**
			  This function reads the entries of the database for the operator
			  and can also be used to change the parameters if required.
			  */
			void reset(const boost::shared_ptr<OperatorParameters>& );

			/**
			  For this operator the apply function would map the solution by injunction from 
			  Simple Vector in u to NodalScalar Vector in r. 
			  @param [in]  f auxillary vector. 
			  @param [in]  u input vector. 
			  @param [out] r output vector. 
			  @param [in]  a first constant used in the expression: r = a*A(u) + b*f. The default value is -1.
			  @param [in]  b second constant used in the expression: r = a*A(u) + b*f. The default value is 1.
			  */
			void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
					AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

			AMP::LinearAlgebra::Variable::shared_ptr createInputVariable (const std::string & name, int  = -1)
			{
				return d_inpVariable->cloneVariable(name);
			}

			AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable (const std::string & name, int  = -1) 
			{
				return d_outVariable->cloneVariable(name);
			}

			AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1) {
				return d_inpVariable;
			}

			AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
				return d_outVariable;
			}

			/**
			  This member function is used to compute 1D locations.
			  */
			void computeZLocations();

			/**
			  This member function returns the number of 1D locations.
			  */
			size_t  getNumZlocations(){return d_zLocations.size();}

			/**
			  This member function returns the 1D locations stl vector.
			  */
			std::vector<double> getZLocations()
			{
				return d_zLocations ;
			}

			void setVector (AMP::LinearAlgebra::Vector::shared_ptr vec)
			{
				outputVec = vec->subsetVectorForVariable(d_outVariable);
			}


		protected :

			AMP::LinearAlgebra::Vector::shared_ptr outputVec;

			boost::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;

			boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable; 

			std::vector<double> d_zLocations; /**< std vector to store 1D z locations. */


		private :


	};

}
}

#endif
