#include "ArraySize.h"


static_assert( AMP::ArraySize( { 10 } ).ndim() == 1 );
static_assert( AMP::ArraySize( { 10 } ).length() == 10 );
static_assert( AMP::ArraySize( { 10, 20 } ).ndim() == 2 );
static_assert( AMP::ArraySize( { 10, -1 } ).ndim() == 1 );
static_assert( AMP::ArraySize( { 10, 20 } ).length() == 200 );
static_assert( AMP::ArraySize( { 10, -1 } ).length() == 10 );
