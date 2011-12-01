
#ifndef included_AMP_PelletStackOperator
#define included_AMP_PelletStackOperator

#include "operators/Operator.h"

namespace AMP {
  namespace Operator {

    class PelletStackOperator : public Operator 
    {
      public :
        PelletStackOperator(const boost::shared_ptr<OperatorParameters> & params);

        ~PelletStackOperator() { }

        bool hasPellet(unsigned int pellId); 

        void setCurrentPellet(unsigned int pellId);

        unsigned int getTotalNumberOfPellets();

        bool useSerial();

        bool onlyZcorrection();

        bool useScaling();

        void applyScaling(AMP::LinearAlgebra::Vector::shared_ptr f);

        void applyUnscaling(AMP::LinearAlgebra::Vector::shared_ptr f);

        void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
            AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1);

      protected:

        bool d_useSerial;
        bool d_onlyZcorrection;
        bool d_useScaling;
        double d_scalingFactor;

    };

  }
}

#endif


