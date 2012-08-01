#include "RowOperator.h"

namespace AMP {
namespace Operator {


void RowOperator::apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
	const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
	const double a, const double b)
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
       d_Operators[i]->apply(fNull, u, rInternal[i], scalea[i], 0);
     }

     for(unsigned int i = 0; i < d_Operators.size(); i++)
     {
       rOriginal->add(rOriginal, (rInternal[i]) );
     }
     
     if(f.get()==NULL)
     {
       rOriginal->scale(a);
     }
     else
     {
       r->axpby(b, a, *f);
     }

}

}
}

