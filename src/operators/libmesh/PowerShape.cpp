
//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   operators/PowerShape.cc
 * \author Kevin Clarno and Gokhan Yesilyurt
 * \brief  Power shape
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: PowerShape.cc,v 1.3 2010/06/15 12:00:00 h51 Exp $
//---------------------------------------------------------------------------//

#include "PowerShape.h"
#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/Operator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/Database.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "PowerShapeParameters.h"
#include "VolumeIntegralOperator.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/string_to_enum.h"
#include <memory>

#include <cmath>
#include <vector>

namespace AMP {
namespace Operator {

/*!
 *************************************************************************
 * \brief Constructor for PowerShape.  The constructor initializes the  *
 * values from the parameters.                                          *
 *************************************************************************
 */
PowerShape::PowerShape( std::shared_ptr<PowerShapeParameters> parameters ) : Operator( parameters )
{
    AMP_ASSERT( parameters );
    d_Mesh = parameters->d_Mesh;
    d_db   = parameters->d_db;
    getFromDatabase( d_db );
}

/*!
 *************************************************************************
 * \brief Destructor.                                                    *
 *************************************************************************
 */
PowerShape::~PowerShape() = default;

/*!
 *************************************************************************
 * \brief Resets all of the available coordinate system variables to the *
 *  default values.                                                      *
 *************************************************************************
 */
void PowerShape::reset( const std::shared_ptr<OperatorParameters> &parameters )
{
    AMP_ASSERT( parameters.get() != nullptr );
    d_db = parameters->d_db;

    if ( d_coordinateSystem == "cartesian" ) {

        if ( d_type == "legendre" ) {
            d_useFixedPower = true;
            d_numXmoments   = 0;
            d_numYmoments   = 0;
            d_numZmoments   = 0;
        } else if ( d_type == "gaussian" ) {
            d_muX    = 0.;
            d_muY    = 0.;
            d_sigmaX = 3;
            d_sigmaY = 3;
        } else {
            AMP_INSIST(
                0, "The power shape type used is not valid for cartesian coordinate systems." );
        }
    } else if ( d_coordinateSystem == "cylindrical" ) {

        if ( d_type == "frapcon" ) {
            //      AMP_ASSERT(!(d_type == "frapcon"))
            d_numZmoments     = 0;
            d_frapconConstant = 3.45;
            d_angularConstant = 0.0;
        } else if ( d_type == "zernikeRadial" ) {
            d_numMoments  = 0;
            d_numZmoments = 0;
        } else if ( d_type == "zernike" ) {
            d_numZmoments = 0;
            d_numMoments  = 0;
        } else if ( d_type == "diffusion" ) {
            d_numZmoments = 0;
        } else {
            AMP_INSIST(
                0, "The power shape type used is not valid for cylindrical coordinate systems." );
        }
    } else {
        AMP_INSIST( 0, "The coordinate system is not valid." );
    }
}

/*!
 ****************************************************************************
 * \brief If simulation is not from restart, read data from input database. *
 * Otherwise, override restart values for a subset of the data members      *
 * with those found in input.                                               *
 ****************************************************************************
 */
void PowerShape::getFromDatabase( std::shared_ptr<AMP::Database> db )
{
    AMP_ASSERT( db );

    // d_Variable = createOutputVariable("RelativePower");

    // Coordinate System
    d_coordinateSystem = db->getWithDefault<std::string>( "coordinateSystem", "cartesian" );

    // Create the cylindrical bounding box
    auto min_max_pos = d_Mesh->getBoundingBox();
    double centerx   = 0.5 * ( min_max_pos[0] + min_max_pos[1] );
    double centery   = 0.5 * ( min_max_pos[2] + min_max_pos[3] );
    double minR = 1e100, maxR = -1e100;

    // Create the cylindrical bounding box
    d_radialBoundingBox = min_max_pos;
    auto iterator       = d_Mesh->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
    for ( size_t i = 0; i < iterator.size(); i++ ) {
        auto coord = iterator->coord();
        double rx  = ( coord[0] - centerx );
        double ry  = ( coord[1] - centery );
        minR       = std::min( minR, sqrt( rx * rx + ry * ry ) );
        maxR       = std::max( maxR, sqrt( rx * rx + ry * ry ) );
        ++iterator;
    }
    d_radialBoundingBox[0] = centerx;
    d_radialBoundingBox[1] = centery;
    d_radialBoundingBox[2] = ( d_Mesh->getComm() ).minReduce( minR );
    d_radialBoundingBox[3] = ( d_Mesh->getComm() ).maxReduce( maxR );
    d_radialBoundingBox[4] = min_max_pos[4];
    d_radialBoundingBox[5] = min_max_pos[5];

    if ( d_coordinateSystem == "cartesian" ) {

        d_type = db->getWithDefault<std::string>( "type", "legendre" );
        if ( d_type == "legendre" ) {

            // Number of moments in the X direction. Default = 0
            d_numXmoments = db->getWithDefault( "numXmoments", 0 );
            AMP_ASSERT( (int) d_numXmoments > -1 );
            if ( d_numXmoments > 0 ) {
                AMP_ASSERT( db->keyExists( "Xmoments" ) );
                d_Xmoments = db->getVector<double>( "Xmoments" );
                for ( unsigned int i = 0; i < d_numXmoments; i++ ) {
                    AMP_ASSERT( fabs( d_Xmoments[i] ) <= 1.0 );
                }
            }

            // Number of moments in the Y direction. Default = 0
            d_numYmoments = db->getWithDefault( "numYmoments", 0 );
            AMP_ASSERT( (int) d_numYmoments > -1 );
            if ( d_numYmoments > 0 ) {
                AMP_ASSERT( db->keyExists( "Ymoments" ) );
                d_Ymoments = db->getVector<double>( "Ymoments" );
                for ( unsigned int i = 0; i < d_numYmoments; i++ ) {
                    AMP_ASSERT( fabs( d_Ymoments[i] ) <= 1.0 );
                }
            }

            // Number of moments in the Z direction. Default = 0
            d_numZmoments = db->getWithDefault( "numZmoments", 0 );
            AMP_ASSERT( (int) d_numZmoments > -1 );
            if ( d_numZmoments > 0 ) {
                AMP_ASSERT( db->keyExists( "Zmoments" ) );
                d_Zmoments = db->getVector<double>( "Zmoments" );
                for ( unsigned int i = 0; i < d_numZmoments; i++ ) {
                    AMP_ASSERT( fabs( d_Zmoments[i] ) <= 1.0 );
                }
            }
        } else if ( d_type == "gaussian" ) {

            min_max_pos = d_Mesh->getBoundingBox();
            centerx     = 0.5 * ( min_max_pos[0] + min_max_pos[1] );
            centery     = 0.5 * ( min_max_pos[2] + min_max_pos[3] );

            // Read mu and sigma for gaussian distribution.
            d_muX    = db->getWithDefault<double>( "muX", centerx );
            d_muY    = db->getWithDefault<double>( "muY", centery );
            d_sigmaX = db->getWithDefault<double>( "sigmaX", 3.0 );
            d_sigmaY = db->getWithDefault<double>( "sigmaY", 3.0 );

            // Number of moments in the Z direction. Default = 0
            d_numZmoments = db->getWithDefault( "numZmoments", 0 );
            AMP_ASSERT( (int) d_numZmoments > -1 );
            if ( d_numZmoments > 0 ) {
                AMP_ASSERT( db->keyExists( "Zmoments" ) );
                d_Zmoments = db->getVector<double>( "Zmoments" );
                for ( unsigned int i = 0; i < d_numZmoments; i++ ) {
                    AMP_ASSERT( fabs( d_Zmoments[i] ) <= 1.0 );
                }
            }
        } else {
            AMP_INSIST(
                0, "The power shape type used is not valid for cartesian coordinate systems." );
        }
    } else if ( d_coordinateSystem == "cylindrical" ) {

        // Number of moments in the Z direction. Default =0
        d_numZmoments = db->getWithDefault( "numZmoments", 0 );
        AMP_ASSERT( (int) d_numZmoments > -1 );
        if ( d_numZmoments > 0 ) {
            AMP_ASSERT( db->keyExists( "Zmoments" ) );
            d_Zmoments = db->getVector<double>( "Zmoments" );
            for ( unsigned int i = 0; i < d_numZmoments; i++ ) {
                AMP_ASSERT( fabs( d_Zmoments[i] ) <= 1.0 );
            }
        }

        // Read the powershope model from input database.
        d_type = db->getWithDefault<std::string>( "type", "frapcon" );

        // frapcon power shape
        if ( d_type == "frapcon" ) {

            // Read the Frapcon constant from input database.
            d_frapconConstant = db->getWithDefault<double>( "frapconConstant", 3.45 );

            // Read the angular constant from input database.
            d_angularConstant = db->getWithDefault<double>( "angularConstant", 0.0 );
            AMP_ASSERT( fabs( d_angularConstant ) <= 1.0 );
        } else if ( d_type == "zernikeRadial" ) {

            // Number of moments in the Radial direction. Default =0
            d_numMoments = db->getWithDefault( "numMoments", 4 );
            d_Moments.resize( d_numMoments );

            // Establish default values
            if ( d_numMoments == 4 ) {
                d_Moments.resize( 4 );
                d_Moments[0] = 0.215;
                d_Moments[1] = 0.00879;
                d_Moments[2] = 0.0000335;
                d_Moments[3] = 0.0000450;
            }

            // these are only the m=even, n=zero moments.
            AMP_ASSERT( (int) d_numMoments > -1 );
            if ( d_numMoments > 0 ) {
                if ( d_numMoments != 4 )
                    AMP_INSIST( db->keyExists( "Moments" ),
                                "if numMoments are not zero, and "
                                "default is not used, you must define "
                                "the Moments to use." );
                if ( db->keyExists( "Moments" ) ) {
                    d_Moments = db->getVector<double>( "Moments" );
                    for ( unsigned int i = 0; i < d_numMoments; i++ ) {
                        AMP_ASSERT( fabs( d_Moments[i] ) <= 1.0 );
                    }
                }
            }

            // Read the angular constant from input database.
            d_angularConstant = db->getWithDefault<double>( "angularConstant", 0.0 );
            AMP_ASSERT( fabs( d_angularConstant ) <= 1.0 );
        } else if ( d_type == "zernike" ) {

            // Number of "m" moments in the zernike polynomial. Default =0
            d_numMoments = db->getWithDefault( "numMoments", 0 );
            // m=1; n=-1, 1
            // m=2; n=-2, 0, 2
            // m=3; n=-3, -1, 1, 3
            d_numMNmoments = ( d_numMoments + 2 ) * ( d_numMoments + 1 ) / 2 - 1;

            // these are only the m=even, n=zero moments.
            AMP_ASSERT( (int) d_numMoments > -1 );
            if ( d_numMoments > 0 ) {
                AMP_ASSERT( db->keyExists( "Moments" ) );
                d_Moments = db->getVector<double>( "Moments" );
                for ( unsigned int i = 0; i < d_numMNmoments; i++ ) {
                    AMP_ASSERT( fabs( d_Moments[i] ) <= 1.0 );
                }
            }
        } else if ( d_type == "diffusion" ) {
            d_numZmoments = 0;
        } else {
            AMP_INSIST(
                0, "The power shape type used is not valid for cylindrical coordinate systems!" );
        }
    } else {
        AMP_INSIST( 0, "The coordinate system is not valid." );
    }

    // create a shape function.
    {
        // std::string feTypeOrderName = d_db->getWithDefault<std::string>("FE_ORDER", "SECOND");
        std::string feTypeOrderName = d_db->getWithDefault<std::string>( "FE_ORDER", "FIRST" );
        auto feTypeOrder = libMesh::Utility::string_to_enum<libMeshEnums::Order>( feTypeOrderName );
        std::string feFamilyName = d_db->getWithDefault<std::string>( "FE_FAMILY", "LAGRANGE" );
        auto feFamily = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( feFamilyName );
        std::string qruleTypeName = d_db->getWithDefault<std::string>( "QRULE_TYPE", "QGAUSS" );
        auto qruleType =
            libMesh::Utility::string_to_enum<libMeshEnums::QuadratureType>( qruleTypeName );
        const unsigned int dimension = 3;
        d_feType.reset( new libMesh::FEType( feTypeOrder, feFamily ) );
        d_fe.reset( ( libMesh::FEBase::build( dimension, ( *d_feType ) ) ).release() );
        d_fe->get_xyz();
        std::string qruleOrderName = d_db->getWithDefault<std::string>( "QRULE_ORDER", "DEFAULT" );
        libMeshEnums::Order qruleOrder;
        if ( qruleOrderName == "DEFAULT" ) {
            qruleOrder = d_feType->default_quadrature_order();
        } else {
            qruleOrder = libMesh::Utility::string_to_enum<libMeshEnums::Order>( qruleOrderName );
        }

        d_qrule.reset( ( libMesh::QBase::build( qruleType, dimension, qruleOrder ) ).release() );
        d_fe->attach_quadrature_rule( d_qrule.get() );
    }
}

/*!
 *************************************************************************
 * \brief Provide the Apply function for the power shape. Power shape   *
 * is chosen based on the coordinate system and the model type.         *
 * Available power shapes:                                              *
 *                         Cartesian:   Legendre                        *
 *                                      Gaussian                        *
 *                         Cylindrical: Frapcon                         *
 *                                      Polynomial (not implemented yet)*
 *                                                                      *
 * Once the coordinate system and the model type is chosen, exact shape *
 * is calculated by using the model parameters that are input by the    *
 * user. Then, the resultant shape function is applied to a constant    *
 * power by keeping the total power the same over the problem domain.   *
 *************************************************************************
 */
void PowerShape::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr r )
{

    AMP_INSIST( ( ( u.get() ) != nullptr ), "NULL Power Vector" );
    AMP_INSIST( ( ( r.get() ) != nullptr ), "NULL PowerWithShape Vector" );

    constexpr double PI = 3.14159265359;
    double newval, val;
    int countGP = 0;

    // quick exit if the answer is obvious.
    if ( ( u->max() < 1e-14 ) && ( u->min() > -1e-14 ) ) {
        r->setToScalar( 0. );
        return;
    }

    auto min_max_pos = d_Mesh->getBoundingBox();

    double xmin    = min_max_pos[0];
    double xmax    = min_max_pos[1];
    double ymin    = min_max_pos[2];
    double ymax    = min_max_pos[3];
    double zmin    = d_radialBoundingBox[4];
    double zmax    = d_radialBoundingBox[5];
    double rmax    = d_radialBoundingBox[3];
    double centerx = d_radialBoundingBox[0];
    double centery = d_radialBoundingBox[1];

    // Create a DOF manager for a gauss point vector
    int DOFsPerElement = 8;
    int ghostWidth     = 0;
    bool split         = true;
    auto dof_map       = AMP::Discretization::simpleDOFManager::create(
        d_Mesh, AMP::Mesh::GeomType::Volume, ghostWidth, DOFsPerElement, split );

    // Create a shared pointer to a Variable - Power - Output because it will be used in the
    // "residual" location of
    // apply.
    r->setToScalar( 1. );

    auto elem      = d_Mesh->getIterator( AMP::Mesh::GeomType::Volume, ghostWidth );
    auto end_elems = elem.end();

    if ( d_coordinateSystem == "cartesian" ) {

        if ( d_type == "legendre" ) {

            if ( d_iDebugPrintInfoLevel > 3 )
                AMP::pout << "Starting Power Shape Loop over Gauss Points." << std::endl;

            // Loop over all elements on the mesh
            for ( ; elem != end_elems; ++elem ) {
                d_currNodes = elem->getElements( AMP::Mesh::GeomType::Vertex );
                createCurrentLibMeshElement();
                d_fe->reinit( d_currElemPtr );

                // Loop over all gauss-points on the element.
                for ( int i = 0; i < DOFsPerElement; i++ ) {
                    double x = d_fe->get_xyz()[i]( 0 );
                    double y = d_fe->get_xyz()[i]( 1 );
                    double z = d_fe->get_xyz()[i]( 2 );

                    // X moments
                    x = ( 2 * x - ( xmax + xmin ) ) /
                        ( xmax - xmin ); // map x coordinate to (-1 to 1)
                    newval = 1.;
                    for ( unsigned int m = 0; m < d_numXmoments; m++ ) {
                        newval += d_Xmoments[m] * evalLegendre( m + 1, x );
                    } // end for
                    val = newval;

                    // Y moments
                    y = ( 2 * y - ( ymax + ymin ) ) /
                        ( ymax - ymin ); // map y coordinate to (-1 to 1)
                    newval = 1.;
                    for ( unsigned int m = 0; m < d_numYmoments; m++ ) {
                        newval += d_Ymoments[m] * evalLegendre( m + 1, y );
                    } // end for m
                    val *= newval;

                    // Z moments
                    z = ( 2 * z - ( zmax + zmin ) ) /
                        ( zmax - zmin ); // map z coordinate to (-1 to 1)
                    newval = 1.;
                    for ( unsigned int m = 0; m < d_numZmoments; m++ ) {
                        newval += d_Zmoments[m] * evalLegendre( m + 1, z );
                    } // end for m
                    val *= newval;

                    countGP++;

                    // Set newval to this gauss point on this element of the vector r.
                    std::vector<size_t> ndx;
                    dof_map->getDOFs( elem->globalID(), ndx );
                    size_t offset = ndx[i];
                    val *= u->getValueByGlobalID( offset );
                    r->setValuesByGlobalID( 1, &offset, &val );
                    AMP_ASSERT(
                        AMP::Utilities::approx_equal( r->getValueByGlobalID( offset ), val ) );
                } // end for gauss-points
                destroyCurrentLibMeshElement();
            } // end for elements
            r->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

            if ( d_iDebugPrintInfoLevel > 3 )
                AMP::pout << "End Power Shape Loop over : " << countGP << " Gauss Points."
                          << std::endl;
        } else if ( d_type == "gaussian" ) {

            if ( d_iDebugPrintInfoLevel > 3 )
                AMP::pout << "Power Shape: Processing all Gauss-Points." << std::endl;
            // Loop over all elements on the mesh
            for ( ; elem != end_elems; ++elem ) {
                d_currNodes = elem->getElements( AMP::Mesh::GeomType::Vertex );
                createCurrentLibMeshElement();
                d_fe->reinit( d_currElemPtr );

                // Loop over all gauss-points on the element.
                for ( int i = 0; i < DOFsPerElement; i++ ) {
                    double x = d_fe->get_xyz()[i]( 0 );
                    double y = d_fe->get_xyz()[i]( 1 );
                    double z = d_fe->get_xyz()[i]( 2 );

                    // 2D Gaussian (Normal) distribution.
                    newval = ( xmax - xmin ) * ( ymax - ymin ) * getGaussianF( x, y );
                    newval =
                        newval / ( ( erf( -( xmax - d_muX ) / ( sqrt( 2. ) * d_sigmaX ) ) *
                                         erf( -( ymax - d_muY ) / ( sqrt( 2. ) * d_sigmaY ) ) +
                                     erf( -( xmin - d_muX ) / ( sqrt( 2. ) * d_sigmaX ) ) *
                                         erf( -( ymin - d_muY ) / ( sqrt( 2. ) * d_sigmaY ) ) -
                                     erf( -( xmax - d_muX ) / ( sqrt( 2. ) * d_sigmaX ) ) *
                                         erf( -( ymin - d_muY ) / ( sqrt( 2. ) * d_sigmaY ) ) -
                                     erf( -( xmin - d_muX ) / ( sqrt( 2. ) * d_sigmaX ) ) *
                                         erf( -( ymax - d_muY ) / ( sqrt( 2. ) * d_sigmaY ) ) ) *
                                   d_sigmaX * d_sigmaY * PI / 2.0 );
                    val = newval;

                    // Z moments
                    z = ( 2 * z - ( zmax + zmin ) ) /
                        ( zmax - zmin ); // Mapping z coordinate to (-1 to 1)
                    newval = 1.;
                    for ( unsigned int m = 0; m < d_numZmoments; m++ ) {
                        newval += d_Zmoments[m] * evalLegendre( m + 1, z );
                    } // end for m
                    val *= newval;

                    countGP++;

                    // Set newval to this gauss point on this element of the vector r.
                    std::vector<size_t> ndx;
                    dof_map->getDOFs( elem->globalID(), ndx );
                    size_t offset = ndx[i];
                    val *= u->getValueByGlobalID( offset );
                    r->setValuesByGlobalID( 1, &offset, &val );
                    AMP_ASSERT(
                        AMP::Utilities::approx_equal( r->getValueByGlobalID( offset ), val ) );
                } // end for gauss-points
                destroyCurrentLibMeshElement();
            } // end for elements
            r->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

            if ( d_iDebugPrintInfoLevel > 3 )
                AMP::pout << "Power Shape: Processing GP #: " << countGP << std::endl;
        } else {
            AMP_INSIST(
                0, "The power shape type used is not valid for cylindrical coordinate systems." );
        }
    } else if ( d_coordinateSystem == "cylindrical" ) {

        if ( d_type == "frapcon" ) {

            // Note: Dimensions are all in meter (m).

            // Choose the type of volume integral calculation.
            double volumeIntegral = getVolumeIntegralSum( rmax, centerx, centery );

            if ( d_iDebugPrintInfoLevel > 3 )
                AMP::pout << "Power Shape: Processing all Gauss-Points." << std::endl;
            // Loop over all elements on the mesh
            for ( ; elem != end_elems; ++elem ) {
                d_currNodes = elem->getElements( AMP::Mesh::GeomType::Vertex );
                createCurrentLibMeshElement();
                d_fe->reinit( d_currElemPtr );

                // Loop over all gauss-points on the element.
                for ( int i = 0; i < DOFsPerElement; i++ ) {
                    double x = d_fe->get_xyz()[i]( 0 ) - centerx;
                    double y = d_fe->get_xyz()[i]( 1 ) - centery;
                    double z = d_fe->get_xyz()[i]( 2 );

                    // r based on Frapcon.
                    double radius = sqrt( x * x + y * y );
                    double Fr     = getFrapconFr( radius, rmax );
                    newval        = Fr / volumeIntegral;
                    val           = newval;

                    // phi.
                    double theta = atan2( y, x );
                    newval       = 1 + d_angularConstant * sin( theta );
                    val *= newval;

                    // Z moments.
                    z = ( 2 * z - ( zmax + zmin ) ) /
                        ( zmax - zmin ); // mapping z coordinate to (-1 to 1)
                    newval = 1.;
                    for ( unsigned int m = 0; m < d_numZmoments; m++ ) {
                        newval += d_Zmoments[m] * evalLegendre( m, z );
                    } // end for m
                    val *= newval;

                    countGP++;

                    // Set newval to this gauss point on this element of the vector r.
                    std::vector<size_t> ndx;
                    dof_map->getDOFs( elem->globalID(), ndx );
                    size_t offset = ndx[i];
                    val *= u->getValueByGlobalID( offset );
                    r->setValuesByGlobalID( 1, &offset, &val );
                    AMP_ASSERT(
                        AMP::Utilities::approx_equal( r->getValueByGlobalID( offset ), val ) );
                } // end for gauss-points
                destroyCurrentLibMeshElement();
            } // end for elements
            r->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

            if ( d_iDebugPrintInfoLevel > 3 )
                AMP::pout << "Power Shape: Processing GP #: " << countGP << std::endl;
        } else if ( d_type == "diffusion" ) {

            // Infinite cylinder diffusion shape
            // Note: Dimensions are all in meter (m).
            if ( d_iDebugPrintInfoLevel > 3 )
                AMP::pout << "Power Shape: Processing all Gauss-Points." << std::endl;
            // Loop over all elements on the mesh
            for ( ; elem != end_elems; ++elem ) {
                d_currNodes = elem->getElements( AMP::Mesh::GeomType::Vertex );
                createCurrentLibMeshElement();
                d_fe->reinit( d_currElemPtr );

                // Loop over all gauss-points on the element.
                for ( int i = 0; i < DOFsPerElement; i++ ) {
                    double x = d_fe->get_xyz()[i]( 0 ) - centerx;
                    double y = d_fe->get_xyz()[i]( 1 ) - centery;
                    // r based on Frapcon.
                    double relativeRadius = sqrt( x * x + y * y ) / rmax;
                    double besArg         = 2.405 * relativeRadius;
                    // Taylor expansion of bessel function

                    val = 1 + ( besArg * besArg ) / 4 + ( besArg * besArg * besArg * besArg ) / 64 +
                          ( besArg * besArg * besArg * besArg * besArg * besArg ) / 2304;


                    // Set newval to this gauss point on this element of the vector r.
                    std::vector<size_t> ndx;
                    dof_map->getDOFs( elem->globalID(), ndx );
                    size_t offset = ndx[i];
                    val *= u->getValueByGlobalID( offset );
                    r->setValuesByGlobalID( 1, &offset, &val );
                    AMP_ASSERT(
                        AMP::Utilities::approx_equal( r->getValueByGlobalID( offset ), val ) );
                } // end for gauss-points
                destroyCurrentLibMeshElement();
            } // end for elements
            r->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

            AMP::LinearAlgebra::Vector::shared_ptr nullVec;

            // this is the volume integral operator that will be used to make the integral of the
            // power density the
            // value the user specified.
            d_db->putDatabase( "VolumeIntegral" );
            std::shared_ptr<AMP::Database> volume_db = d_db->getDatabase( "VolumeIntegral" );
            std::shared_ptr<AMP::Database> act_db;

            volume_db->putScalar( "name", "VolumeIntegralOperator" );

            volume_db->putScalar( "InputVariableType", "IntegrationPointScalar" );
            volume_db->putScalar( "Number_Active_Variables", 1 );
            volume_db->putScalar( "Number_Auxillary_Variables", 0 );
            volume_db->putScalar( "Constant_Source", true );
            volume_db->putScalar( "OutputVariable", "Temperature" );
            volume_db->putScalar( "print_info_level", 1 );
            volume_db->putDatabase( "ActiveInputVariables" );
            volume_db->putDatabase( "SourceElement" );

            std::shared_ptr<AMP::Database> source_db = volume_db->getDatabase( "SourceElement" );
            source_db->putScalar( "name", "SourceNonlinearElement" );

            act_db = volume_db->getDatabase( "ActiveInputVariables" );
            act_db->putScalar( "ActiveVariable_0", ( u->getVariable() )->getName() );
            std::shared_ptr<AMP::Operator::ElementPhysicsModel> emptyModel;
            std::shared_ptr<AMP::Operator::VolumeIntegralOperator> volumeIntegralOperator =
                std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
                    AMP::Operator::OperatorBuilder::createOperator(
                        d_Mesh, "VolumeIntegral", d_db, emptyModel ) );

            int DOFsPerNode     = 1;
            int nodalGhostWidth = 1;
            auto nodalDofMap    = AMP::Discretization::simpleDOFManager::create(
                d_Mesh, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
            auto nodalVariable = std::make_shared<AMP::LinearAlgebra::Variable>( "Temperature" );
            auto nodalVector =
                AMP::LinearAlgebra::createVector( nodalDofMap, nodalVariable, split );
            auto unodalPower = nodalVector->cloneVector();
            auto rnodalPower = nodalVector->cloneVector();

            volumeIntegralOperator->apply( u, unodalPower );
            volumeIntegralOperator->apply( r, rnodalPower );

            double denominator( rnodalPower->L1Norm() );
            AMP_INSIST( !AMP::Utilities::approx_equal( denominator, 0. ),
                        "The denominator is zero - not good." );
            r->scale( double( unodalPower->L1Norm() ) / double( rnodalPower->L1Norm() ) );

            r->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
            if ( d_iDebugPrintInfoLevel > 3 )
                AMP::pout << "Power Shape: Processed GP #: " << countGP << std::endl;
        }

        else if ( d_type == "zernikeRadial" ) {

            // Note: Dimensions are all in meter (m).

            if ( d_iDebugPrintInfoLevel > 3 )
                AMP::pout << "Power Shape: Processing all Gauss-Points." << std::endl;
            // Loop over all elements on the mesh
            for ( ; elem != end_elems; ++elem ) {
                d_currNodes = elem->getElements( AMP::Mesh::GeomType::Vertex );
                createCurrentLibMeshElement();
                d_fe->reinit( d_currElemPtr );

                // Loop over all gauss-points on the element.
                for ( int i = 0; i < DOFsPerElement; i++ ) {
                    double x = d_fe->get_xyz()[i]( 0 ) - centerx;
                    double y = d_fe->get_xyz()[i]( 1 ) - centery;
                    double z = d_fe->get_xyz()[i]( 2 );

                    // r based on Frapcon.
                    double relativeRadius = sqrt( x * x + y * y ) / rmax;
                    val                   = 1 + getZernikeRadial( relativeRadius );

                    // phi.
                    double theta = atan2( y, x );
                    newval       = 1 + d_angularConstant * sin( theta );
                    val *= newval;

                    // Z moments.
                    z = ( 2 * z - ( zmax + zmin ) ) /
                        ( zmax - zmin ); // mapping z coordinate to (-1 to 1)
                    newval = 1.;
                    for ( unsigned int m = 0; m < d_numZmoments; m++ ) {
                        newval += d_Zmoments[m] * evalLegendre( m + 1, z );
                    } // end for m
                    val *= newval;

                    countGP++;

                    // Set newval to this gauss point on this element of the vector r.
                    std::vector<size_t> ndx;
                    dof_map->getDOFs( elem->globalID(), ndx );
                    size_t offset = ndx[i];
                    val *= u->getValueByGlobalID( offset );
                    r->setValuesByGlobalID( 1, &offset, &val );
                    AMP_INSIST(
                        AMP::Utilities::approx_equal( r->getValueByGlobalID( offset ), val ),
                        "getting what was just put" );
                } // end for gauss-points
                destroyCurrentLibMeshElement();
            } // end for elements
            r->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

            if ( d_iDebugPrintInfoLevel > 3 )
                AMP::pout << "Power Shape: Processed GP #: " << countGP << std::endl;
        } else if ( d_type == "zernike" ) {

            // Note: Dimensions are all in meter (m).

            if ( d_iDebugPrintInfoLevel > 3 )
                AMP::pout << "Power Shape: Processing all Gauss-Points." << std::endl;
            // Loop over all elements on the mesh
            for ( ; elem != end_elems; ++elem ) {
                d_currNodes = elem->getElements( AMP::Mesh::GeomType::Vertex );
                createCurrentLibMeshElement();
                d_fe->reinit( d_currElemPtr );

                // Loop over all gauss-points on the element.
                for ( int i = 0; i < DOFsPerElement; i++ ) {
                    double x = d_fe->get_xyz()[i]( 0 ) - centerx;
                    double y = d_fe->get_xyz()[i]( 1 ) - centery;
                    double z = d_fe->get_xyz()[i]( 2 );

                    // r based on Frapcon.
                    double relativeRadius = sqrt( x * x + y * y ) / rmax;
                    double phi            = atan2( y, x );
                    val                   = 1 + getZernike( relativeRadius, phi );

                    // Z moments.
                    z = ( 2 * z - ( zmax + zmin ) ) /
                        ( zmax - zmin ); // mapping z coordinate to (-1 to 1)
                    newval = 1.;
                    for ( unsigned int m = 0; m < d_numZmoments; m++ ) {
                        newval += d_Zmoments[m] * evalLegendre( m + 1, z );
                    } // end for m
                    val *= newval;

                    countGP++;

                    // Set newval to this gauss point on this element of the vector r.
                    std::vector<size_t> ndx;
                    dof_map->getDOFs( elem->globalID(), ndx );
                    size_t offset = ndx[i];
                    val *= u->getValueByGlobalID( offset );
                    r->setValuesByGlobalID( 1, &offset, &val );
                    AMP_ASSERT(
                        AMP::Utilities::approx_equal( r->getValueByGlobalID( offset ), val ) );
                } // end for gauss-points
                destroyCurrentLibMeshElement();
            } // end for elements
            r->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

            if ( d_iDebugPrintInfoLevel > 3 )
                AMP::pout << "Power Shape: Processed GP #: " << countGP << std::endl;
        } else {
            AMP_INSIST(
                0, "The power shape type used is not valid for cylindrical coordinate systems" );
        }
    } else {
        AMP_INSIST( 0, "The coordinate system is not valid." );
    }
}

/*!
 *************************************************************************
 * \brief Print class data members to given output stream.               *
 *************************************************************************
 */
void PowerShape::printClassData( std::ostream &os ) const
{
    os << "\nPowerShape::printClassData..." << std::endl;
    os << "d_numXmoments = " << d_numXmoments << std::endl;
}

/*!
 *************************************************************************
 * \brief Write out class version number and data members to database.   *
 *************************************************************************
 */
void PowerShape::putToDatabase( std::shared_ptr<AMP::Database> db )
{
    AMP_ASSERT( !db.use_count() );
    db->putScalar( "numXmoments", d_numXmoments );
    db->putVector( "Xmoments", d_Xmoments );
}

/*!
 *************************************************************************
 * \brief Evaluates the value of the Legendre polynomial Pn at x         *
 *************************************************************************
 */
double PowerShape::evalLegendre( const int n, const double x )
{
    double result = 1;
    double a      = 1;
    double b      = x;
    if ( n == 0 ) {
        result = a;
        return result;
    }
    if ( n == 1 ) {
        result = b;
        return result;
    }
    for ( int i = 2; i <= n; i++ ) {
        result = ( ( 2 * i - 1 ) * x * b - ( i - 1 ) * a ) / i;
        a      = b;
        b      = result;
    }
    return result;
}

/*!
 ****************************************************************************
 * \brief Evaluates the volume integral by Sum ( Sum( f(r) )*elemVolume/8 )*
 ****************************************************************************
 */
double PowerShape::getVolumeIntegralSum( double rmax, double cx, double cy )
{
    double integralFr            = 0;
    double numerator             = 0;
    int ghostWidth               = 0;
    AMP::Mesh::MeshIterator elem = d_Mesh->getIterator( AMP::Mesh::GeomType::Volume, ghostWidth );
    AMP::Mesh::MeshIterator end_elems = elem.end();
    for ( ; elem != end_elems; ++elem ) {
        d_currNodes = elem->getElements( AMP::Mesh::GeomType::Vertex );
        createCurrentLibMeshElement();
        d_fe->reinit( d_currElemPtr );
        double elemVolume = elem->volume();

        double elemSum = 0;
        // Loop over all gauss-points on the element.
        int DOFsPerElement = 8;
        for ( int i = 0; i < DOFsPerElement; i++ ) {
            double x      = d_fe->get_xyz()[i]( 0 ) - cx;
            double y      = d_fe->get_xyz()[i]( 1 ) - cy;
            double radius = sqrt( x * x + y * y );
            elemSum += getFrapconFr( radius, rmax );
        } // end for gauss-points
        integralFr += ( elemSum / 8.0 ) * elemVolume;
        numerator += elemVolume;
        destroyCurrentLibMeshElement();
    } // end for elements

    double bot = ( d_Mesh->getComm() ).sumReduce( numerator );
    double top = ( d_Mesh->getComm() ).sumReduce( integralFr );
    // integralFr = integralFr/numerator;
    integralFr = top / bot;

    return integralFr;
}

void PowerShape::createCurrentLibMeshElement()
{
    d_currElemPtr = new libMesh::Hex8;
    for ( unsigned int j = 0; j < d_currNodes.size(); j++ ) {
        auto pt                      = d_currNodes[j].coord();
        d_currElemPtr->set_node( j ) = new libMesh::Node( pt[0], pt[1], pt[2], j );
    } // end for j
}

void PowerShape::destroyCurrentLibMeshElement()
{
    for ( unsigned int j = 0; j < d_currElemPtr->n_nodes(); j++ ) {
        delete ( d_currElemPtr->node_ptr( j ) );
        d_currElemPtr->set_node( j ) = nullptr;
    } // end for j
    delete d_currElemPtr;
    d_currElemPtr = nullptr;
}

/*!
 *************************************************************************
 * \brief Evaluates Frapcon power shape F(r).                           *
 *************************************************************************
 */
double PowerShape::getFrapconFr( double radius, double rmax )
{
    double fR = ( 1.0 + d_frapconConstant * exp( -3.0 * pow( 1000.0 * ( rmax - radius ), 0.45 ) ) );
    return fR;
}

/*!
 *************************************************************************
 * \brief Evaluates Zernike Radial power shape F(r).                     *
 * \param rhor Relative radius.                                          *
 * d_Moments are the coeffiecients for the (2*i=2, 0) moments.           *
 *************************************************************************
 */
double PowerShape::getZernikeRadial( double rhor )
{
    double fR = 0.0;
    for ( unsigned int j = 0; j < d_numMoments; j++ ) {
        unsigned int m = 2 * j + 2;
        fR += d_Moments[j] * evalZernike( 0, m, rhor, 0. );
    }
    return fR;
}

/*!
 *************************************************************************
 * \brief Evaluates Zernike power shape F(r,phi).                       *
 *************************************************************************
 */
double PowerShape::getZernike( const double rhor, const double phi )
{
    double rho = rhor;
    double fR  = 0.0;
    int j      = 0;
    for ( int m = 1; m <= (int) d_numMoments; m++ ) {
        for ( int n = -m; n <= m; n += 2 ) {
            fR += d_Moments[j] * evalZernike( n, m, rho, phi );
            j++;
        }
    }

    if ( d_iDebugPrintInfoLevel > 8 ) {
        if ( d_numMoments < 5 ) {
            double Zr = fR;
            fR        = 0;
            if ( d_numMoments > 0 ) {
                fR += d_Moments[0] * rho * sin( phi );
                fR += d_Moments[1] * rho * cos( phi );
            }
            if ( d_numMoments > 1 ) {
                double rho2 = rho * rho;
                double phi2 = 2 * phi;
                fR += d_Moments[2] * rho2 * sin( phi2 );
                fR += d_Moments[3] * ( 2. * rho2 - 1. );
                fR += d_Moments[4] * rho2 * cos( phi2 );
            }
            if ( d_numMoments > 2 ) {
                double rho2 = rho * rho;
                double rho3 = rho2 * rho;
                double phi3 = 3 * phi;
                fR += d_Moments[5] * rho3 * sin( phi3 );
                fR += d_Moments[6] * rho * ( 3. * rho2 - 2. ) * sin( phi );
                fR += d_Moments[7] * rho * ( 3. * rho2 - 2. ) * cos( phi );
                fR += d_Moments[8] * rho3 * cos( phi3 );
            }
            if ( d_numMoments > 3 ) {
                double rho2 = rho * rho;
                double rho4 = rho2 * rho2;
                double phi2 = 2 * phi;
                double phi4 = 2 * phi2;
                fR += d_Moments[9] * rho4 * sin( phi4 );
                fR += d_Moments[10] * rho2 * ( 4 * rho2 - 3 ) * sin( phi2 );
                fR += d_Moments[11] * ( 1 - 6 * rho2 * ( 1. - rho2 ) );
                fR += d_Moments[12] * rho2 * ( 4 * rho2 - 3 ) * cos( phi2 );
                fR += d_Moments[13] * rho4 * cos( phi4 );
            }
            if ( d_iDebugPrintInfoLevel > 9 ) {
                if ( !AMP::Utilities::approx_equal( Zr, fR, 1e-8 ) ) {
                    AMP::pout << "getZernike fails for rho = " << rho << " and phi = " << phi
                              << " for m = " << d_numMoments << std::endl;
                    AMP::pout << "fR is " << fR << " Zr is " << Zr << std::endl;
                    if ( d_numMoments == 3 ) {
                        AMP::pout << "Moments are " << d_Moments[5] << " " << d_Moments[6] << " "
                                  << d_Moments[7] << " " << d_Moments[8] << std::endl;
                    }
                }
            } // if d_iDebugPrintInfoLevel>9

            AMP_ASSERT( AMP::Utilities::approx_equal( Zr, fR, 1e-8 ) );

        } // if d_numMoments<5

    } // if d_iDebugPrintInfoLevel>8

    return fR;
}

/*!
 *************************************************************************
 * \brief Computes a factorial, which is used by the Zernike basis.     *
 *************************************************************************
 */
double PowerShape::evalFactorial( int n )
{
    double f = 1;
    while ( n > 1 )
        f *= n--;
    return f;
}

/*!
 *************************************************************************
 * \brief Computes the combinatorial choose function
 *************************************************************************
 */
double PowerShape::choose( int n, int k )
{
    if ( k > n )
        return 0;

    int r = 1;
    for ( int d = 1; d <= k; ++d ) {
        r *= n--;
        r /= d;
    }
    return double( r );
}

/*!
 *************************************************************************
 * \brief Computes the Zernike value for (m,n) of (rho, phi).           *
 *************************************************************************
 */
double PowerShape::evalZernike( int m, int n, const double rho, const double phi )
{
    // For Steven's sake, n and m are reversed from Denovo usage
    std::swap( n, m );

    double rhoFact = 0;
    bool even      = n >= 0;
    if ( n < 0 )
        n = -n;
    if ( ( m - n ) % 2 )
        return 0;

    // This should be slightly more stable numerically than
    //  multiplying/dividing factorials, probably doesn't matter
    //  but can't hurt
    for ( int k = 0; k <= ( m - n ) / 2; ++k ) {
        rhoFact += pow( -1, k ) * pow( rho, m - 2 * k ) * choose( m - k, k ) *
                   choose( m - 2 * k, ( n + m ) / 2 - k );
    }

    // Apply normalization (should be consistent with Denovo?)
    /*    double pi = 4.0*atan(1.0);
          if( m==0 )
          rhoFact *= sqrt( (n+1)/pi );
          else
          rhoFact *= sqrt( 2*(n+1)/pi );
          */
    /*    for(int k=0;k<=(m-n)/2;++k) {
          rhoFact += pow(rho,m-2*k) * ( pow(-1,k) * evalFactorial(m-k) )
          / (             evalFactorial(k)
     * evalFactorial((m+n)/2-k)
     * evalFactorial((m-n)/2-k) );
     } */

    if ( even )
        return rhoFact * cos( n * phi );
    else
        return rhoFact * sin( n * phi );
}

/*!
 *************************************************************************
 * \brief Evaluates 2D Gaussian (Normal) distribution.                  *
 *************************************************************************
 */
double PowerShape::getGaussianF( double x, double y )
{
    double gaussianF = exp( -( pow( x - d_muX, 2 ) / ( 2 * pow( d_sigmaX, 2.0 ) ) +
                               pow( y - d_muY, 2 ) / ( 2 * pow( d_sigmaY, 2.0 ) ) ) );
    return gaussianF;
}
} // namespace Operator
} // end namespace AMP

//---------------------------------------------------------------------------//
//                 end of PowerShape.cc
//---------------------------------------------------------------------------//
