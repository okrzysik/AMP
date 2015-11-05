#include "RowOperator.h"

namespace AMP {
namespace Operator {


void RowOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u, 
			 AMP::LinearAlgebra::Vector::shared_ptr r )
{
  
     AMP::LinearAlgebra::Vector::shared_ptr fNull;

     d_OutputVariable = (d_Operators[0])->getOutputVariable();

     std::vector<AMP::LinearAlgebra::Vector::shared_ptr> rInternal(d_Operators.size());
     AMP::LinearAlgebra::Vector::shared_ptr rOriginal = r->subsetVectorForVariable(d_OutputVariable);

     rOriginal->zero();

     for(unsigned int i = 0; i < d_Operators.size(); i++)
     {
       rInternal[i] = rOriginal->cloneVector(); 
       rInternal[i]->zero();
       d_Operators[i]->apply(u, rInternal[i]);
       rInternal[i]->scale(scalea[i]);
     }

     for(unsigned int i = 0; i < d_Operators.size(); i++)
     {
       rOriginal->add(rOriginal, (rInternal[i]) );
     }
     
}

}
}

