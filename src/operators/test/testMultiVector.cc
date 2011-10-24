#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "materials/Material.h"

#include "vectors/MultiVariable.h"
#include "vectors/Vector.h"

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshVariable.h"
#include "ampmesh/SiloIO.h"

#include "libmesh.h"


extern "C"{
#ifdef MPICH_SKIP_MPICXX
#define _FIX_FOR_PETSC_MPICH_CXX
#undef MPICH_SKIP_MPICXX
#endif

#ifdef OMPI_SKIP_MPICXX
#define _FIX_FOR_PETSC_OMPI_CXX
#undef OMPI_SKIP_MPICXX
#endif

#include "petsc.h"

#ifdef _FIX_FOR_PETSC_MPICH_CXX
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#endif

#ifdef _FIX_FOR_PETSC_OMPI_CXX
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif
#endif
}

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

void myTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "brick" );

  //--------------------------------------------------------------------------------
  //Construct Variables
  boost::shared_ptr<AMP::LinearAlgebra::Variable> Variable1 (new AMP::Mesh::NodalScalarVariable ("Var1", meshAdapter));
  boost::shared_ptr<AMP::LinearAlgebra::Variable> Variable2 (new AMP::Mesh::NodalScalarVariable ("Var2", meshAdapter));
  boost::shared_ptr<AMP::LinearAlgebra::Variable> Variable3 (new AMP::Mesh::NodalScalarVariable ("Var3", meshAdapter));
  boost::shared_ptr<AMP::LinearAlgebra::Variable> Variable4 (new AMP::Mesh::NodalScalarVariable ("Var4", meshAdapter));

  //-------------------------------------------------------

  boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> subVariable (new AMP::LinearAlgebra::MultiVariable("subVar"));
  subVariable->add( Variable1   );
  subVariable->add( Variable2   );

  //-
  boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> fullVariable(new AMP::LinearAlgebra::MultiVariable("fullVariable"));
  fullVariable->add( Variable1   );
  fullVariable->add( Variable2   );
  fullVariable->add( Variable3   );
  fullVariable->add( Variable4   );

  AMP::LinearAlgebra::Vector::shared_ptr FullVector = manager->createVector ( fullVariable );

  AMP::LinearAlgebra::Vector::shared_ptr SubVector  = FullVector->subsetVectorForVariable( subVariable );

  AMP_INSIST( SubVector.get()!=NULL , "Sub Vector is NULL");


  ut->passes(exeName);

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back("testMultiVector");
    for(size_t i = 0; i < exeNames.size(); i++) {
        try {
            myTest(&ut, exeNames[i]);
        } catch (std::exception &err) {
            AMP::pout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
            ut.failure("ERROR");
        } catch( ... ) {
            AMP::pout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
            ut.failure("ERROR");
        }
    } //end for i

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   


