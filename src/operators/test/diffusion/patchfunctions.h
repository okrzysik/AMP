/*
 * patchfunctions.h
 *
 *  Created on: Jun 22, 2010
 *      Author: gad
 */

#ifndef PATCHFUNCTIONS_H_
#define PATCHFUNCTIONS_H_

inline double constant( const double, const double, const double ) { return 1.; }
inline double x_linear( const double x, const double, const double ) { return x; }
inline double y_linear( const double, const double y, const double ) { return y; }
inline double z_linear( const double, const double, const double z ) { return z; }
inline double x_squared( const double x, const double, const double ) { return x * x; }
inline double y_squared( const double, const double y, const double ) { return y * y; }
inline double z_squared( const double, const double, const double z ) { return z * z; }
inline double r_squared( const double x, const double y, const double z )
{
    return x * x + y * y + z * z;
}
inline double radius( const double x, const double y, const double z )
{
    return std::sqrt( r_squared( x, y, z ) );
}

#endif /* PATCHFUNCTIONS_H_ */
