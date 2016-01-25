/*
 * testManufacturedSolution.cc
 *
 *  Created on: Jul 27, 2010
 *      Author: gad
 */

#include <exception>
#include <iostream>
#include <valarray>

#include "../Database.h"
#include "../ManufacturedSolution.h"
#include "../MemoryDatabase.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/shared_ptr.h"


void testit( AMP::UnitTest *ut,
             std::string geom,
             std::string order,
             std::string bc,
             double x,
             double y,
             double z )
{
    AMP::shared_ptr<AMP::Database> db(
        AMP::dynamic_pointer_cast<AMP::Database>( AMP::shared_ptr<AMP::MemoryDatabase>(
            new AMP::MemoryDatabase( "ManufacturedSolution" ) ) ) );

    db->putString( "Geometry", geom );
    db->putString( "Order", order );
    db->putString( "BoundaryType", bc );
    db->putDouble( "MinX", 4. );
    db->putDouble( "MaxX", 14. );
    db->putDouble( "MinY", 40. );
    db->putDouble( "MaxY", 80. );
    db->putDouble( "MinZ", 400. );
    db->putDouble( "MaxZ", 1000. );
    db->putDouble( "MinR", 30. );
    db->putDouble( "MaxR", 100. );
    db->putDouble( "MinTh", 0. );
    db->putDouble( "MaxTh", 6.284 );

    AMP::ManufacturedSolution ms( db );
    size_t nc = ms.getNumberOfInputs();
    size_t na = ms.getNumberOfParameters();

    std::valarray<double> out( 10 ), in( nc ), param( na );

    for ( size_t i = 0; i < na; i++ )
        param[i]   = i;
    for ( size_t i = 0; i < nc; i++ )
        in[i]      = i + 1;
    ms.setTricubicParams( in, param );
    ms.evaluate( out, x, y, z );

    std::string msg = geom + " " + order + " " + bc + " basics";
    ut->passes( msg );
}


int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    AMP::AMP_MPI globalComm = AMP::AMP_MPI( AMP_COMM_WORLD );

    if ( globalComm.getRank() == 0 ) {
        testit( &ut, "Brick", "Quadratic", "Neumann", 5., 60., 700. );
        testit( &ut, "Brick", "Quadratic", "Dirichlet-1", 5., 60., 700. );
        testit( &ut, "Brick", "Quadratic", "Dirichlet-2", 5., 60., 700. );
        testit( &ut, "Brick", "Cubic", "Neumann", 5., 60., 700. );
        testit( &ut, "Brick", "Cubic", "Dirichlet-1", 5., 60., 700. );
        testit( &ut, "Brick", "Cubic", "Dirichlet-2", 5., 60., 700. );
        testit( &ut, "CylindricalRod", "Cubic", "None", 55., 4.2, 700. );
        testit( &ut, "CylindricalRod", "Cubic", "Dirichlet-2-z", 55., 4.2, 700. );
        testit( &ut, "CylindricalShell", "Quadratic", "Neumann", 55., 4.2, 700. );
    } else {
        ut.expected_failure( "Manufactured Solutions only apply to scalar tests." );
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
