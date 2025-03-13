#include "AMP/utils/copycast/CopyCastHelper.h"
#include "CopyCastHelper.hpp"

namespace AMP::Utilities {

template struct copyCast_<double,
                          float,
                          AMP::Utilities::PortabilityBackend::Kokkos,
                          AMP::HostAllocator<void>>;
template struct copyCast_<float,
                          double,
                          AMP::Utilities::PortabilityBackend::Kokkos,
                          AMP::HostAllocator<void>>;
template struct copyCast_<float,
                          float,
                          AMP::Utilities::PortabilityBackend::Kokkos,
                          AMP::HostAllocator<void>>;
template struct copyCast_<double,
                          double,
                          AMP::Utilities::PortabilityBackend::Kokkos,
                          AMP::HostAllocator<void>>;

#ifdef USE_DEVICE
template struct copyCast_<double,
                          float,
                          AMP::Utilities::PortabilityBackend::Kokkos,
                          AMP::ManagedAllocator<void>>;
template struct copyCast_<float,
                          double,
                          AMP::Utilities::PortabilityBackend::Kokkos,
                          AMP::ManagedAllocator<void>>;
template struct copyCast_<float,
                          float,
                          AMP::Utilities::PortabilityBackend::Kokkos,
                          AMP::ManagedAllocator<void>>;
template struct copyCast_<double,
                          double,
                          AMP::Utilities::PortabilityBackend::Kokkos,
                          AMP::ManagedAllocator<void>>;

template struct copyCast_<double,
                          float,
                          AMP::Utilities::PortabilityBackend::Kokkos,
                          AMP::DeviceAllocator<void>>;
template struct copyCast_<float,
                          double,
                          AMP::Utilities::PortabilityBackend::Kokkos,
                          AMP::DeviceAllocator<void>>;
template struct copyCast_<float,
                          float,
                          AMP::Utilities::PortabilityBackend::Kokkos,
                          AMP::DeviceAllocator<void>>;
template struct copyCast_<double,
                          double,
                          AMP::Utilities::PortabilityBackend::Kokkos,
                          AMP::DeviceAllocator<void>>;
#endif


} // namespace AMP::Utilities
