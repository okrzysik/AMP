// This file contains a list of enums for units until we have a Units class
#ifndef included_AMP_Units
#define included_AMP_Units

namespace AMP {
namespace Units {

enum class Time { unknown = 0, seconds, minutes, hours, days };
enum class Length { unknown = 0, meters };
enum class Temperature { unknown = 0, kelvin };
enum class Mass { unknown = 0, kilograms };
}
}

#endif
