#pragma once

#include <hdf5.h>
#include <mpi.h>

/**
 * A template function that returns the appropriate types to allow for templated HDF5 functions
 * @tparam T
 * @return
 */
template <typename T>
hid_t HdfType();

template <>
inline hid_t HdfType<short>()
{
    return H5T_NATIVE_SHORT;
}
template <>
inline hid_t HdfType<unsigned short>()
{
    return H5T_NATIVE_USHORT;
}
template <>
inline hid_t HdfType<int>()
{
    return H5T_NATIVE_INT;
}
template <>
inline hid_t HdfType<unsigned int>()
{
    return H5T_NATIVE_UINT;
}
template <>
inline hid_t HdfType<long>()
{
    return H5T_NATIVE_LONG;
}
template <>
inline hid_t HdfType<unsigned long>()
{
    return H5T_NATIVE_ULONG;
}
template <>
inline hid_t HdfType<long long>()
{
    return H5T_NATIVE_LLONG;
}
template <>
inline hid_t HdfType<unsigned long long>()
{
    return H5T_NATIVE_ULLONG;
}

template <>
inline hid_t HdfType<char>()
{
    return H5T_NATIVE_CHAR;
}
template <>
inline hid_t HdfType<unsigned char>()
{
    return H5T_NATIVE_UCHAR;
}
template <>
inline hid_t HdfType<signed char>()
{
    return H5T_NATIVE_SCHAR;
}

template <>
inline hid_t HdfType<float>()
{
    return H5T_NATIVE_FLOAT;
}
template <>
inline hid_t HdfType<double>()
{
    return H5T_NATIVE_DOUBLE;
}
template <>
inline hid_t HdfType<long double>()
{
    return H5T_NATIVE_LDOUBLE;
}

template <typename T>
MPI_Datatype MpiType();

template <>
inline MPI_Datatype MpiType<short>()
{
    return MPI_SHORT;
}
template <>
inline MPI_Datatype MpiType<unsigned short>()
{
    return MPI_UNSIGNED_SHORT;
}
template <>
inline MPI_Datatype MpiType<int>()
{
    return MPI_INT;
}
template <>
inline MPI_Datatype MpiType<unsigned int>()
{
    return MPI_UINT32_T;
}
template <>
inline MPI_Datatype MpiType<long>()
{
    return MPI_LONG;
}
template <>
inline MPI_Datatype MpiType<unsigned long>()
{
    return MPI_UNSIGNED_LONG;
}
template <>
inline MPI_Datatype MpiType<long long>()
{
    return MPI_LONG_LONG_INT;
}
template <>
inline MPI_Datatype MpiType<unsigned long long>()
{
    return MPI_UNSIGNED_LONG_LONG;
}

template <>
inline MPI_Datatype MpiType<char>()
{
    return MPI_CHAR;
}
template <>
inline MPI_Datatype MpiType<unsigned char>()
{
    return MPI_UNSIGNED_CHAR;
}
template <>
inline MPI_Datatype MpiType<signed char>()
{
    return MPI_SIGNED_CHAR;
}

template <>
inline MPI_Datatype MpiType<float>()
{
    return MPI_FLOAT;
}
template <>
inline MPI_Datatype MpiType<double>()
{
    return MPI_DOUBLE;
}
template <>
inline MPI_Datatype MpiType<long double>()
{
    return MPI_LONG_DOUBLE;
}