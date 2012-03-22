
#ifndef included_AMP_NeutronicsRhsExtras
#define included_AMP_NeutronicsRhsExtras

/* AMP files */
#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "NeutronicsRhsExtrasParameters.h"
#include "vectors/Variable.h"
#include "utils/Utilities.h"

#include "utils/InputDatabase.h"

/*Boost files */
#include "boost/shared_ptr.hpp"

#include <vector>

namespace AMP {
namespace Operator {

//===========================================================================//
/*!
 * \class NeutronicsRhsExtras 
 * \brief A class for representing the neutronics source operator.
 */
//===========================================================================//

  class NeutronicsRhsExtras : public  Operator {

    public:
      //typedef AMP::LinearAlgebra::VectorVariable<AMP::Mesh::IntegrationPointVariable, 8>      HexGaussPointVariable;
      //typedef boost::shared_ptr<HexGaussPointVariable>      SP_HexGaussPointVariable;
      typedef boost::shared_ptr<NeutronicsRhsExtrasParameters>               SP_Parameters;
      typedef boost::shared_ptr<OperatorParameters>            SP_OperatorParameters;
      typedef boost::shared_ptr<AMP::LinearAlgebra::Vector>                SP_Vector; 
      typedef std::vector<double>                                            Vec_Dbl;
      typedef boost::shared_ptr<Vec_Dbl>                                  SP_Vec_Dbl; 
      typedef boost::shared_ptr<AMP::Database>                           SP_Database;
      typedef std::vector<double>                                           Vec_Dbl1;
      typedef std::vector<Vec_Dbl1>                                         Vec_Dbl2;
      typedef std::vector<Vec_Dbl2>                                         Vec_Dbl3;
    
      //! Neutronics Input Types
      enum SourceType{ Isotopes, Elements, NUM_SOURCE_TYPES };
      
    private:
      
    public:

      NeutronicsRhsExtras(SP_Parameters parameters);

      /**
       * Empty destructor for NeutronicsRhsExtras
       */
      virtual ~NeutronicsRhsExtras();

      /**
       * Print out all members of integrator instance to given output stream.
       */
      void printClassData(std::ostream& os) const;

      /**
       * Write out state of object to given database.
       *
       * When assertion checking is active, the database pointer must be non-null.
       */
      void putToDatabase(SP_Database db);

      /**
        The function that computes the residual.
       * @param f: rhs vector for A(u)=f, this may be a null pointer if f=0. 
       * @param u: multivector of the state.
       * @param r: specific power in Watts per gram 
       The result of apply is
       * r = b*f+a*A(u)
       */
      void apply(const  SP_Vector & f, 
                 const  SP_Vector & u, 
                        SP_Vector & r,
                 const  double      a = 1.0,
                 const  double      b = 0.0);

      /**
        A function to reinitialize this object.
        */
      void reset(const SP_OperatorParameters & parameters);

      AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {return d_outputVariable;}

      void   setTimeStep ( int tStep ) { d_timeStep = tStep; }
      void   setExtrasId ( int extrasId ) { d_extrasId = extrasId; }
      void   setTimeInSeconds ( double seconds );

    protected:

      /*
       * Read input data from specified database and initialize class members.
       * If run is from restart, a subset of the restart values may be replaced
       * with those read from input.
       *
       * When assertion checking is active, the database pointer must be non-null.
       */
      void getFromInput(SP_Database db);

      SP_Database                               d_db;
      bool                                      d_useFixedValue;
      int                                       d_numTimeSteps;
      Vec_Dbl                                   d_timeStepsInDays;
      SourceType                                d_type;      
      Vec_Dbl                                   d_fixedValues;
      int                                       d_timeStep;
      int                                       d_extrasId;
      int                                       d_numExtras;
      std::vector<std::string>                  d_extrasName;
      double                                    d_timeStepInSeconds;
      AMP::LinearAlgebra::Variable::shared_ptr  d_outputVariable;
      Vec_Dbl3                                  d_values;
      AMP::Mesh::Mesh::shared_ptr               d_Mesh;
      double                                    d_secondsPerDay;               
      SourceType str2id(std::string str);
  };

}
}

#include "NeutronicsRhsExtras.i.h"

#endif



