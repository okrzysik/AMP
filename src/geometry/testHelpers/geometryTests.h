#ifndef included_AMP_GeometryTests
#define included_AMP_GeometryTests

#include "AMP/geometry/Geometry.h"
#include "AMP/utils/UnitTest.h"


namespace AMP::Geometry {


/**
 * \brief Run all geometry based tests
 * \details  This test runs all the geometry-based tests
 * \param[in] geom          Geometry to test
 * \param[in] ut            Unit test
 */
void testGeometry( const AMP::Geometry::Geometry &geom, AMP::UnitTest &ut );


} // namespace AMP::Geometry

#endif
