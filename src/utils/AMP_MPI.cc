// This file impliments a wrapper class for MPI functions

// Include mpi.h if used
#ifdef USE_EXT_MPI
    #include "mpi.h"
#endif

// Include AMP headers
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"

// Include all other headers
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Utilities.h"

// Include timers (if MPI is not used)
#ifndef USE_EXT_MPI
    #if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        // We are using windows without mpi, we need windows timers
        #include <windows.h>
    #else
        // We are using linux without mpi, we need linux timers
        #include <sys/time.h>
    #endif
    // Global variable to track create new unique comms (dup and split)
    MPI_Comm uniqueGlobalComm=11;
#endif

// Platform-dependent routine to release the current thread for the scheduler
// This is useful if we are waiting for some event and want to allow other 
// threads/processes to run while we are waiting
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    // We are using windows, use Sleep(0)
    #include <windows.h>
    #define sched_yield() Sleep(0)
#else
    // We are using linux, use the sched_yield function
    #include <sched.h>
#endif

// Some special structs to work with MPI
struct IntIntStruct { int j; int i; };
struct LongIntStruct { long int j; int i; };
struct FloatIntStruct { float f; int i; };
struct DoubleIntStruct { double d; int i; };


namespace AMP{


// Initialized the static member variables
volatile unsigned int AMP_MPI::N_MPI_Comm_created=0;

// Static data for asyncronous communication without MPI
#ifndef USE_EXT_MPI
static const int mpi_max_tag = 0x003FFFFF;
struct Isendrecv_struct {
    const char *data;   // Pointer to data
    int status;         // Status: 1-sending, 2-recieving
};
std::map<MPI_Request,Isendrecv_struct>  global_isendrecv_list;
static MPI_Request getRequest( MPI_Comm comm, int tag )
{
    AMP_ASSERT(tag>=0&&tag<=mpi_max_tag);
    MPI_Request request = 0;
    if ( sizeof(MPI_Request)==4 && sizeof(MPI_Comm)==4 ) {
        AMP_ASSERT(sizeof(unsigned int)==4);
        unsigned int hash = comm*0x9E3779B9;                        // 2^32*0.5*(sqrt(5)-1)
        unsigned int key  = (hash&0xFFFFFFFF)>>22;                  // Get a key 0-1024
        request = (MPI_Request) tag + (key<<22);
    } else if ( sizeof(MPI_Request)==8 && sizeof(MPI_Comm)==4 ) {
        AMP_ASSERT(sizeof(unsigned long int)==8);
        request = (MPI_Request) ((unsigned long)tag + (((unsigned long)comm)<<32));
    } else if ( sizeof(MPI_Request)==4 && sizeof(MPI_Comm)==8 ) {
        AMP_ASSERT(sizeof(unsigned long int)==8);
        unsigned long int hash = comm*0x9E3779B97F4A7C15;           // 2^64*0.5*(sqrt(5)-1)
        unsigned long int key  = (hash&0xFFFFFFFF)>>22;             // Get a key 0-1024
        request = (MPI_Request) ((unsigned int)tag + (key<<32));
    } else if ( sizeof(MPI_Request)==8 && sizeof(MPI_Comm)==8 ) {
        AMP_ASSERT(sizeof(unsigned long int)==8);
        unsigned long int hash = comm*0x9E3779B97F4A7C15;           // 2^64*0.5*(sqrt(5)-1)
        unsigned long int key  = (hash&0xFFFFFFFF);                 // Get a key 0-2^32-1
        request = (MPI_Request) ((unsigned int)tag + (key<<32));
    } else {
        char text[50];
        sprintf(text,"Not Programmed (%i,%i)",(int)sizeof(MPI_Request),(int)sizeof(MPI_Comm));
        AMP_ERROR(std::string(text));
    }
    return request;
}
#endif


/************************************************************************
*  Empty constructor                                                    *
************************************************************************/
AMP_MPI::AMP_MPI() {
    // Initialize the data members to a defaul communicator of self
    #ifdef USE_EXT_MPI
        communicator = MPI_COMM_NULL;
        d_maxTag = 0x7FFFFFFF;
    #else
        communicator = AMP_COMM_NULL;
        d_maxTag = mpi_max_tag;
    #endif
    count = NULL;
    comm_rank = 0;
    comm_size = 1;
    d_isNull = true;
    call_abort_in_serial_instead_of_exit = true;
}


/************************************************************************
*  Empty deconstructor                                                  *
************************************************************************/
AMP_MPI::~AMP_MPI() {
    /*#ifdef USE_EXT_MPI
        // Temporary solution until memory bugs are sorted, remove!!
        if ( count == NULL ) {
        } else if ( *count == 1 ) {
            int flag;
            MPI_Finalized( &flag );
            if ( flag ) {
                printf("Warning: still destroying objects after MPI_Finialize\n");
                return;
            }
        }
    #endif*/
    // Check if the count is == 1
    if ( count == NULL ) {
        // We are not keeping a count, no need to free anything
    } else if ( *count > 1 ) {
        // We need to decrement count, but do not need to free anything
        (*count)--;
    } else if ( *count == 1 ) {
        // We are holding that last reference to the MPI_Comm object, we need to free it
        #ifdef USE_EXT_MPI
            *count = 0;
            delete count;
            count = NULL;
            int err = MPI_Comm_free(&communicator);
            if ( err != MPI_SUCCESS )
                AMP_ERROR("Problem free'ing MPI_Comm object");
            communicator = AMP_COMM_NULL;
        #else
            AMP_ERROR("Internal Error (why do we have a count in serial)");
        #endif
    } else {
        // The count is invalid
        AMP_ERROR("Invalid count");
    }
    comm_rank = 0;
    comm_size = 1;
    d_maxTag = 0;
    d_isNull = true;
    call_abort_in_serial_instead_of_exit = true;
}


/************************************************************************
*  Constructor from existing AMP_MPI object                             *
************************************************************************/
AMP_MPI::AMP_MPI( const AMP::AMP_MPI& comm ) {
    // Initialize the data members to the existing AMP_MPI object
    communicator = comm.communicator;
    comm_rank = comm.comm_rank;
    comm_size = comm.comm_size;
    d_isNull = comm.d_isNull;
    d_maxTag = comm.d_maxTag;
    call_abort_in_serial_instead_of_exit = comm.call_abort_in_serial_instead_of_exit;
    // Set and increment the count
    count = comm.count;
    if ( count != NULL )
        (*count)++;
}


/************************************************************************
*  Constructor from existing MPI communicator                           *
************************************************************************/
AMP_MPI::AMP_MPI( MPI_Comm comm ) {
    #ifdef USE_EXT_MPI
        // We are using MPI, use the MPI communicator to initialize the data
        if ( comm==AMP_COMM_WORLD ) {
            communicator = AMP::AMPManager::comm_world.communicator;
        } else if ( comm==AMP_COMM_SELF ) {
            communicator = MPI_COMM_SELF;
        } else if ( comm==AMP_COMM_NULL ) {
            communicator = MPI_COMM_NULL;
        } else {
            communicator = comm;
        }
        if ( communicator!=MPI_COMM_NULL) {
            // Attach the error handler
            MPI_Comm_set_errhandler( communicator, AMP::AMPManager::mpierr );
            // Get the communicator properties
            MPI_Comm_rank(communicator, &comm_rank);
            MPI_Comm_size(communicator, &comm_size);
            int flag, *val;
            int ierr = MPI_Comm_get_attr(communicator,MPI_TAG_UB,&val,&flag);
            AMP_ASSERT(ierr==MPI_SUCCESS);
            if ( flag==0 ) { 
                d_maxTag = 0x7FFFFFFF;     // The tag is not a valid attribute (set to 2^31-1)
            } else {
                d_maxTag = *val;
                if ( d_maxTag<0 ) { d_maxTag = 0x7FFFFFFF; }    // The maximum tag is > a signed int (set to 2^31-1)
                AMP_INSIST(d_maxTag>=0x7FFF,"maximum tag size is < MPI standard");
            }
        } else {
            comm_rank = 1;
            comm_size = 1;
            d_maxTag = 0x7FFFFFFF;
        }
        d_isNull = communicator==MPI_COMM_NULL;
    #else
        // We are not using MPI, intialize based on the communicator
        communicator = comm;
        comm_rank = 0;
        comm_size = 1;
        d_maxTag = mpi_max_tag;
        d_isNull = communicator==AMP_COMM_NULL;
    #endif
    call_abort_in_serial_instead_of_exit = true;
    // We are creating a AMP_MPI comm from an MPI_Comm, the user is responsible for freeing the MPI_Comm object
    count = NULL;
}


/************************************************************************
*  Intersect two communicators                                          *
************************************************************************/
#ifdef USE_EXT_MPI
AMP_MPI AMP_MPI::intersect( const AMP_MPI &comm1, const AMP_MPI &comm2 ) {
    MPI_Group group1, group2, group12;
    MPI_Comm_group ( comm1.communicator, &group1 );
    MPI_Comm_group ( comm2.communicator, &group2 );
    MPI_Group_intersection( group1, group2, &group12 );
    int compare1, compare2;
    MPI_Group_compare ( group1, group12, &compare1 );
    MPI_Group_compare ( group2, group12, &compare2 );
    AMP_MPI new_comm(AMP_COMM_NULL);
    if ( compare1!=MPI_UNEQUAL ) {
        new_comm = comm1;
    } else if ( compare2!=MPI_UNEQUAL ) {
        new_comm = comm2;
    } else {
        MPI_Comm  new_MPI_comm;
        MPI_Comm_create( comm1.communicator, group12, &new_MPI_comm );
        int size;
        MPI_Group_size( group12, &size );
        if ( size > 0 )
            new_comm = AMP_MPI( new_MPI_comm );
    }
    MPI_Group_free( &group1 );
    MPI_Group_free( &group2 );
    MPI_Group_free( &group12 );
    return new_comm;
}
#else
AMP_MPI AMP_MPI::intersect( const AMP_MPI &comm1, const AMP_MPI &comm2 ) {
    if ( comm1.isNull() || comm2.isNull() )
        return AMP_MPI(AMP_COMM_NULL);
    AMP_ASSERT(comm1.comm_size==1&&comm2.comm_size==1);
    return comm1;
}
#endif

/************************************************************************
*  Assignment operator                                                  *
************************************************************************/
AMP_MPI& AMP_MPI::operator=(const AMP::AMP_MPI& comm) {
    if (this == &comm) // protect against invalid self-assignment
        return *this;
    // Destroy a previous AMP_MPI object
    this->~AMP_MPI();
    // Initialize the data members to the existing AMP_MPI object
    this->communicator = comm.communicator;
    this->comm_rank = comm.comm_rank;
    this->comm_size = comm.comm_size;
    this->d_isNull = comm.d_isNull;
    this->d_maxTag = comm.d_maxTag;
    this->call_abort_in_serial_instead_of_exit = comm.call_abort_in_serial_instead_of_exit;
    // Set and increment the count
    this->count = comm.count;
    if ( this->count != NULL )
        (*(this->count))++;
    return *this;
}


/************************************************************************
*  Split an exisiting AMP_MPI object                                    *
************************************************************************/
AMP_MPI AMP_MPI::split( int color, int key ) const {
    MPI_Comm  new_MPI_comm;
    #ifdef USE_EXT_MPI
        // USE MPI to split the communicator
        if ( color==-1 )
            MPI_Comm_split(communicator,MPI_UNDEFINED,key,&new_MPI_comm);
        else
            MPI_Comm_split(communicator,color,key,&new_MPI_comm);
    #else
        if ( color==-1 ) {
            new_MPI_comm = AMP_COMM_NULL;
        } else {
            new_MPI_comm = uniqueGlobalComm;
            uniqueGlobalComm++;
        }
    #endif
    // Create the AMP_MPI object
    AMP_MPI new_comm(new_MPI_comm);
    #ifdef USE_EXT_MPI
        new_comm.d_isNull = new_comm.communicator==MPI_COMM_NULL;
    #else
        new_comm.d_isNull = new_comm.communicator==AMP_COMM_NULL;
    #endif
    new_comm.call_abort_in_serial_instead_of_exit = call_abort_in_serial_instead_of_exit;
    // Create the count
    #ifdef USE_EXT_MPI
        if ( new_comm.communicator != MPI_COMM_NULL ) {
            new_comm.count = new int;
            *(new_comm.count) = 1;
            ++N_MPI_Comm_created;
        }
    #endif
    return new_comm;
}


/************************************************************************
*  Duplicate an exisiting AMP_MPI object                                *
************************************************************************/
AMP_MPI AMP_MPI::dup( ) const {
    MPI_Comm  new_MPI_comm;
    #ifdef USE_EXT_MPI
        // USE MPI to duplicate the communicator
        MPI_Comm_dup(communicator,&new_MPI_comm);
    #else
        new_MPI_comm = uniqueGlobalComm;
        uniqueGlobalComm++;
    #endif
    // Create the AMP_MPI comm
    AMP_MPI new_comm(new_MPI_comm);
    new_comm.d_isNull = d_isNull;
    new_comm.call_abort_in_serial_instead_of_exit = call_abort_in_serial_instead_of_exit;
    // Create the count
    #ifdef USE_EXT_MPI
        if ( new_comm.communicator != MPI_COMM_NULL ) {
            new_comm.count = new int;
            *(new_comm.count) = 1;
            ++N_MPI_Comm_created;
        }
    #endif
    return new_comm;
}


/************************************************************************
*  Overload operator ==                                                 *
************************************************************************/
bool AMP_MPI::operator==(const AMP_MPI &comm) const {
    return communicator==comm.communicator;
}


/************************************************************************
*  Overload operator !=                                                 *
************************************************************************/
bool AMP_MPI::operator!=(const AMP_MPI &comm) const {
    return communicator!=comm.communicator;
}


/************************************************************************
*  Overload operator <                                                  *
************************************************************************/
bool AMP_MPI::operator<(const AMP_MPI &comm) const {
    AMP_ASSERT( !this->d_isNull && !comm.d_isNull );
    bool flag = true;
    // First check if either communicator is NULL
    #ifdef USE_EXT_MPI
        if ( communicator==MPI_COMM_NULL )
            return false;
        if ( comm.communicator==MPI_COMM_NULL )
            flag = false;
    #else
        if ( communicator==AMP_COMM_NULL )
            return false;
        if ( comm.communicator==AMP_COMM_NULL )
            flag = false;
    #endif
    // Check that the size of the other communicator is > the current communicator size
    if ( comm_size >= comm.comm_size )
        flag = false;
    // Check the union of the communicator groups
    // this is < comm iff this group is a subgroup of comm's group
    #ifdef USE_EXT_MPI
        MPI_Group group1, group2, group12;
        MPI_Comm_group ( communicator, &group1 );
        MPI_Comm_group ( comm.communicator, &group2 );
        MPI_Group_union ( group1, group2, &group12 );
        int compare;
        MPI_Group_compare ( group2, group12, &compare );
        if ( compare==MPI_UNEQUAL )
            flag = false;
        MPI_Group_free( &group1 );
        MPI_Group_free( &group2 );
        MPI_Group_free( &group12 );
    #endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce(flag);
}


/************************************************************************
*  Overload operator <=                                                 *
************************************************************************/
bool AMP_MPI::operator<=(const AMP_MPI &comm) const {
    AMP_ASSERT( !this->d_isNull && !comm.d_isNull );
    bool flag = true;
    // First check if either communicator is NULL
    #ifdef USE_EXT_MPI
        if ( communicator==MPI_COMM_NULL )
            return false;
        if ( comm.communicator==MPI_COMM_NULL )
            flag = false;
    #else
        if ( communicator==AMP_COMM_NULL )
            return false;
        if ( comm.communicator==AMP_COMM_NULL )
            flag = false;
    #endif
    // Check that the size of the other communicator is > the current communicator size
    // this is <= comm iff this group is a subgroup of comm's group
    if ( comm_size > comm.comm_size )
        flag = false;
    // Check the unnion of the communicator groups
    #ifdef USE_EXT_MPI
        MPI_Group group1, group2, group12;
        MPI_Comm_group ( communicator, &group1 );
        MPI_Comm_group ( comm.communicator, &group2 );
        MPI_Group_union ( group1, group2, &group12 );
        int compare;
        MPI_Group_compare ( group2, group12, &compare );
        if ( compare==MPI_UNEQUAL )
            flag = false;
        MPI_Group_free( &group1 );
        MPI_Group_free( &group2 );
        MPI_Group_free( &group12 );
    #endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce(flag);
}


/************************************************************************
*  Overload operator >                                                  *
************************************************************************/
bool AMP_MPI::operator>(const AMP_MPI &comm) const {
    bool flag = true;
    // First check if either communicator is NULL
    #ifdef USE_EXT_MPI
        if ( communicator==MPI_COMM_NULL )
            return false;
        if ( comm.communicator==MPI_COMM_NULL )
            flag = false;
    #else
        if ( communicator==AMP_COMM_NULL )
            return false;
        if ( comm.communicator==AMP_COMM_NULL )
            flag = false;
    #endif
    // Check that the size of the other communicator is > the current communicator size
    if ( comm_size <= comm.comm_size )
        flag = false;
    // Check the unnion of the communicator groups
    // this is > comm iff comm's group is a subgroup of this group
    #ifdef USE_EXT_MPI
        MPI_Group group1, group2, group12;
        MPI_Comm_group ( communicator, &group1 );
        MPI_Comm_group ( comm.communicator, &group2 );
        MPI_Group_union ( group1, group2, &group12 );
        int compare;
        MPI_Group_compare ( group1, group12, &compare );
        if ( compare==MPI_UNEQUAL )
            flag = false;
        MPI_Group_free( &group1 );
        MPI_Group_free( &group2 );
        MPI_Group_free( &group12 );
    #endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce(flag);
}


/************************************************************************
*  Overload operator >=                                                 *
************************************************************************/
bool AMP_MPI::operator>=(const AMP_MPI &comm) const {
    bool flag = true;
    // First check if either communicator is NULL
    #ifdef USE_EXT_MPI
        if ( communicator==MPI_COMM_NULL )
            return false;
        if ( comm.communicator==MPI_COMM_NULL )
            flag = false;
    #else
        if ( communicator==AMP_COMM_NULL )
            return false;
        if ( comm.communicator==AMP_COMM_NULL )
            flag = false;
    #endif
    // Check that the size of the other communicator is > the current communicator size
    if ( comm_size < comm.comm_size )
        flag = false;
    // Check the unnion of the communicator groups
    // this is >= comm iff comm's group is a subgroup of this group
    #ifdef USE_EXT_MPI
        MPI_Group group1, group2, group12;
        MPI_Comm_group ( communicator, &group1 );
        MPI_Comm_group ( comm.communicator, &group2 );
        MPI_Group_union ( group1, group2, &group12 );
        int compare;
        MPI_Group_compare ( group1, group12, &compare );
        if ( compare==MPI_UNEQUAL )
            flag = false;
        MPI_Group_free( &group1 );
        MPI_Group_free( &group2 );
        MPI_Group_free( &group12 );
    #endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce(flag);
}


/************************************************************************
*  Duplicate an exisiting AMP_MPI object                                *
************************************************************************/
int AMP_MPI::compare(const AMP_MPI &comm) const {
    if ( communicator==comm.communicator )
        return 1;
    #ifdef USE_EXT_MPI
        int result;
        MPI_Comm_compare( communicator, comm.communicator, &result );
        if ( result==MPI_IDENT )
            return 2;
        else if ( result==MPI_CONGRUENT )
            return 3;
        else if ( result==MPI_SIMILAR )
            return 4;
        else if ( result==MPI_UNEQUAL )
            return 0;
        AMP_ERROR("Unknown results from AMP_Comm_compare");
    #else
        if ( comm.communicator==AMP_COMM_NULL || communicator==AMP_COMM_NULL )
            return 0;
        else 
            return 3;
    #endif
    return 0;
}


/************************************************************************
*  Abort the program.                                                   *
************************************************************************/
void AMP_MPI::setCallAbortInSerialInsteadOfExit(bool flag)
{
    call_abort_in_serial_instead_of_exit = flag;
}
void AMP_MPI::abort() const
{
    #ifdef USE_EXT_MPI
        int initialized=0, finalized=0;
        MPI_Initialized(&initialized);
        MPI_Finalized(&finalized);
        MPI_Comm comm = communicator;
        if ( comm == MPI_COMM_NULL )
            comm = MPI_COMM_WORLD;
        if ( initialized==0 || finalized!=0 ) {
            // MPI is not availible
            exit(-1);
        } else if ( comm_size > 1 ) {
            MPI_Abort(comm, -1);
        } else if ( call_abort_in_serial_instead_of_exit ) {
            MPI_Abort(comm, -1);
        } else {
            exit(-1);
        }
    #else
        exit(-1);
    #endif
}


/************************************************************************
*  allReduce                                                            *
************************************************************************/
bool AMP_MPI::allReduce(const bool value) const {
    bool ret = value;
    if ( comm_size > 1 ) {
        #ifdef USE_EXT_MPI
            MPI_Allreduce( (void*) &value, (void*) &ret, 1, MPI_UNSIGNED_CHAR, MPI_MIN, communicator);
        #else
            AMP_ERROR("This shouldn't be possible");
        #endif
    }
    return ret;
}


/************************************************************************
*  anyReduce                                                            *
************************************************************************/
bool AMP_MPI::anyReduce(const bool value) const {
    bool ret = value;
    if ( comm_size > 1 ) {
        #ifdef USE_EXT_MPI
            MPI_Allreduce( (void*) &value, (void*) &ret, 1, MPI_UNSIGNED_CHAR, MPI_MAX, communicator);
        #else
            AMP_ERROR("This shouldn't be possible");
        #endif
    }
    return ret;
}


/************************************************************************
*  call_sumReduce                                                       *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
// unsigned char
template <>
void AMP_MPI::call_sumReduce<unsigned char>(const unsigned char *send, unsigned char *recv, const int n) const {
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED_CHAR, MPI_SUM, communicator);
}
template <>
void AMP_MPI::call_sumReduce<unsigned char>(unsigned char *x, const int n) const {
    unsigned char *send = x;
    unsigned char *recv = new unsigned char[n];
    MPI_Allreduce( send, recv, n, MPI_UNSIGNED_CHAR, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
}
// char
template <>
void AMP_MPI::call_sumReduce<char>(const char *send, char *recv, const int n) const {
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_SIGNED_CHAR, MPI_SUM, communicator);
}
template <>
void AMP_MPI::call_sumReduce<char>(char *x, const int n) const {
    char *send = x;
    char *recv = new char[n];
    MPI_Allreduce( send, recv, n, MPI_SIGNED_CHAR, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
}
// unsigned int
template <>
void AMP_MPI::call_sumReduce<unsigned int>(const unsigned int *send, unsigned int *recv, const int n) const {
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED, MPI_SUM, communicator);
}
template <>
void AMP_MPI::call_sumReduce<unsigned int>(unsigned int *x, const int n) const {
    unsigned int *send = x;
    unsigned int *recv = new unsigned int[n];
    MPI_Allreduce( send, recv, n, MPI_UNSIGNED, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
}
// int
template <>
void AMP_MPI::call_sumReduce<int>(const int *send, int *recv, const int n) const {
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_INT, MPI_SUM, communicator);
}
template <>
void AMP_MPI::call_sumReduce<int>(int *x, const int n) const {
    int *send = x;
    int *recv = new int[n];
    MPI_Allreduce( send, recv, n, MPI_INT, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
}
// long int
template <>
void AMP_MPI::call_sumReduce<long int>(const long int *send, long int *recv, const int n) const {
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_LONG, MPI_SUM, communicator);
}
template <>
void AMP_MPI::call_sumReduce<long int>(long int *x, const int n) const {
    long int *send = x;
    long int *recv = new long int[n];
    MPI_Allreduce( send, recv, n, MPI_LONG, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
}
// unsigned long int
template <>
void AMP_MPI::call_sumReduce<unsigned long>(const unsigned long *send, unsigned long *recv, const int n) const {
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED_LONG, MPI_SUM, communicator);
}
template <>
void AMP_MPI::call_sumReduce<unsigned long>(unsigned long *x, const int n) const {
    unsigned long int *send = x;
    unsigned long int *recv = new unsigned long int[n];
    MPI_Allreduce( send, recv, n, MPI_UNSIGNED_LONG, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
}
// float
template <>
void AMP_MPI::call_sumReduce<float>(const float *send, float *recv, const int n) const {
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_FLOAT, MPI_SUM, communicator);
}
template <>
void AMP_MPI::call_sumReduce<float>(float *x, const int n) const {
    float *send = x;
    float *recv = new float[n];
    MPI_Allreduce( send, recv, n, MPI_FLOAT, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
}
// double
template <>
void AMP_MPI::call_sumReduce<double>(const double *send, double *recv, const int n) const {
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_DOUBLE, MPI_SUM, communicator);
}
template <>
void AMP_MPI::call_sumReduce<double>(double *x, const int n) const {
    double *send = x;
    double *recv = new double[n];
    MPI_Allreduce( send, recv, n, MPI_DOUBLE, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
}
// std::complex<double>
template <>
void AMP_MPI::call_sumReduce< std::complex<double> >(const std::complex<double> *x, std::complex<double> *y, const int n) const {
    double *send = new double[2*n];
    double *recv = new double[2*n];
    for (int i=0; i<n; i++) {
        send[2*i+0] = real(x[i]);
        send[2*i+1] = imag(x[i]);
    }
    MPI_Allreduce( (void*) send, (void*) recv, 2*n, MPI_DOUBLE, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        y[i] = std::complex<double>(recv[2*i+0],recv[2*i+1]);
    delete [] send;
    delete [] recv;
}
template <>
void AMP_MPI::call_sumReduce< std::complex<double> >(std::complex<double> *x, const int n) const {
    double *send = new double[2*n];
    double *recv = new double[2*n];
    for (int i=0; i<n; i++) {
        send[2*i+0] = real(x[i]);
        send[2*i+1] = imag(x[i]);
    }
    MPI_Allreduce( send, recv, 2*n, MPI_DOUBLE, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = std::complex<double>(recv[2*i+0],recv[2*i+1]);
    delete [] send;
    delete [] recv;
}
#endif


/************************************************************************
*  call_minReduce                                                       *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
// unsigned char
template <>
void AMP_MPI::call_minReduce<unsigned char>(const unsigned char *send, unsigned char *recv, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) 
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED_CHAR, MPI_MIN, communicator);
    else
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
}
template <>
void AMP_MPI::call_minReduce<unsigned char>(unsigned char *x, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        unsigned char *send = x;
        unsigned char *recv = new unsigned char[n];
        MPI_Allreduce( send, recv, n, MPI_UNSIGNED_CHAR, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
    }
}
// char
template <>
void AMP_MPI::call_minReduce<char>(const char *send, char *recv, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) 
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_SIGNED_CHAR, MPI_MIN, communicator);
    else
         AMP_ERROR("Returning the rank of min with char is not supported yet");
}
template <>
void AMP_MPI::call_minReduce<char>(char *x, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        char *send = x;
        char *recv = new char[n];
        MPI_Allreduce( send, recv, n, MPI_SIGNED_CHAR, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         AMP_ERROR("Returning the rank of min with char is not supported yet");
    }
}
// unsigned int
template <>
void AMP_MPI::call_minReduce<unsigned int>(const unsigned int *send, unsigned int *recv, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) 
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED, MPI_MIN, communicator);
    else
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
}
template <>
void AMP_MPI::call_minReduce<unsigned int>(unsigned int *x, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        unsigned int *send = x;
        unsigned int *recv = new unsigned int[n];
        MPI_Allreduce( send, recv, n, MPI_UNSIGNED, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         AMP_ERROR("Returning the rank of min with unsigned int is not supported yet");
    }
}
// int
template <>
void AMP_MPI::call_minReduce<int>(const int *x, int *y, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_INT, MPI_MIN, communicator);
    } else {
         IntIntStruct *recv = new IntIntStruct[n];
         IntIntStruct *send = new IntIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_2INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
template <>
void AMP_MPI::call_minReduce<int>(int *x, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        int *send = x;
        int *recv = new int[n];
        MPI_Allreduce( send, recv, n, MPI_INT, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         IntIntStruct *recv = new IntIntStruct[n];
         IntIntStruct *send = new IntIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_2INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
// unsigned long int
template <>
void AMP_MPI::call_minReduce<unsigned long int>(const unsigned long int *send, unsigned long int *recv, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) 
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED_LONG, MPI_MIN, communicator);
    else
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
}
template <>
void AMP_MPI::call_minReduce<unsigned long int>(unsigned long int *x, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        unsigned long int *send = x;
        unsigned long int *recv = new unsigned long int[n];
        MPI_Allreduce( send, recv, n, MPI_UNSIGNED_LONG, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         AMP_ERROR("Returning the rank of min with unsigned long int is not supported yet");
    }
}

// long int
template <>
void AMP_MPI::call_minReduce<long int>(const long int *x, long int *y, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_LONG, MPI_MIN, communicator);
    } else {
         LongIntStruct *recv = new LongIntStruct[n];
         LongIntStruct *send = new LongIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_LONG_INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
template <>
void AMP_MPI::call_minReduce<long int>(long int *x, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        long int *send = x;
        long int *recv = new long int[n];
        MPI_Allreduce( send, recv, n, MPI_LONG, MPI_MIN, communicator);
        for (long int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         LongIntStruct *recv = new LongIntStruct[n];
         LongIntStruct *send = new LongIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_LONG_INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
// float
template <>
void AMP_MPI::call_minReduce<float>(const float *x, float *y, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_INT, MPI_MIN, communicator);
    } else {
         FloatIntStruct *recv = new FloatIntStruct[n];
         FloatIntStruct *send = new FloatIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].f = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_FLOAT_INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].f;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
template <>
void AMP_MPI::call_minReduce<float>(float *x, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        float *send = x;
        float *recv = new float[n];
        MPI_Allreduce( send, recv, n, MPI_FLOAT, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         FloatIntStruct *recv = new FloatIntStruct[n];
         FloatIntStruct *send = new FloatIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].f = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_FLOAT_INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].f;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
// double
template <>
void AMP_MPI::call_minReduce<double>(const double *x, double *y, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_DOUBLE, MPI_MIN, communicator);
    } else {
         DoubleIntStruct *recv = new DoubleIntStruct[n];
         DoubleIntStruct *send = new DoubleIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].d = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_DOUBLE_INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].d;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
template <>
void AMP_MPI::call_minReduce<double>(double *x, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        double *send = x;
        double *recv = new double[n];
        MPI_Allreduce( send, recv, n, MPI_DOUBLE, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         DoubleIntStruct *recv = new DoubleIntStruct[n];
         DoubleIntStruct *send = new DoubleIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].d = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_DOUBLE_INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].d;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
#endif


/************************************************************************
*  call_maxReduce                                                    *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
// unsigned char
template <>
void AMP_MPI::call_maxReduce<unsigned char>(const unsigned char *send, unsigned char *recv, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) 
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED_CHAR, MPI_MAX, communicator);
    else
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
}template <>
void AMP_MPI::call_maxReduce<unsigned char>(unsigned char *x, const int n, int *comm_rank_of_max) const {
    if ( comm_rank_of_max==NULL ) {
        unsigned char *send = x;
        unsigned char *recv = new unsigned char[n];
        MPI_Allreduce( send, recv, n, MPI_UNSIGNED_CHAR, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
    }
}
// char
template <>
void AMP_MPI::call_maxReduce<char>(const char *send, char *recv, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) 
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_SIGNED_CHAR, MPI_MAX, communicator);
    else
         AMP_ERROR("Returning the rank of min with char is not supported yet");
}
template <>
void AMP_MPI::call_maxReduce<char>(char *x, const int n, int *comm_rank_of_max) const {
    if ( comm_rank_of_max==NULL ) {
        char *send = x;
        char *recv = new char[n];
        MPI_Allreduce( send, recv, n, MPI_SIGNED_CHAR, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         AMP_ERROR("Returning the rank of min with char is not supported yet");
    }
}
// unsigned int
template <>
void AMP_MPI::call_maxReduce<unsigned int>(const unsigned int *send, unsigned int *recv, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) 
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED, MPI_MAX, communicator);
    else
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
}
template <>
void AMP_MPI::call_maxReduce<unsigned int>(unsigned int *x, const int n, int *comm_rank_of_max) const {
    if ( comm_rank_of_max==NULL ) {
        unsigned int *send = x;
        unsigned int *recv = new unsigned int[n];
        MPI_Allreduce( send, recv, n, MPI_UNSIGNED, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         AMP_ERROR("Returning the rank of min with unsigned int is not supported yet");
    }
}
// int
template <>
void AMP_MPI::call_maxReduce<int>(const int *x, int *y, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_INT, MPI_MAX, communicator);
    } else {
         IntIntStruct *recv = new IntIntStruct[n];
         IntIntStruct *send = new IntIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_2INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
template <>
void AMP_MPI::call_maxReduce<int>(int *x, const int n, int *comm_rank_of_max) const {
    if ( comm_rank_of_max==NULL ) {
        int *send = x;
        int *recv = new int[n];
        MPI_Allreduce( send, recv, n, MPI_INT, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         IntIntStruct *recv = new IntIntStruct[n];
         IntIntStruct *send = new IntIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_2INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].j;
            comm_rank_of_max[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
// long int
template <>
void AMP_MPI::call_maxReduce<long int>(const long int *x, long int *y, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_LONG, MPI_MAX, communicator);
    } else {
         LongIntStruct *recv = new LongIntStruct[n];
         LongIntStruct *send = new LongIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_LONG_INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
template <>
void AMP_MPI::call_maxReduce<long int>(long int *x, const int n, int *comm_rank_of_max) const {
    if ( comm_rank_of_max==NULL ) {
        long int *send = x;
        long int *recv = new long int[n];
        MPI_Allreduce( send, recv, n, MPI_LONG, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         LongIntStruct *recv = new LongIntStruct[n];
         LongIntStruct *send = new LongIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_LONG_INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].j;
            comm_rank_of_max[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
// unsigned long int
template <>
void AMP_MPI::call_maxReduce<unsigned long int>(const unsigned long int *send, unsigned long int *recv, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) 
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED_LONG, MPI_MAX, communicator);
    else
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
}
template <>
void AMP_MPI::call_maxReduce<unsigned long int>(unsigned long int *x, const int n, int *comm_rank_of_max) const {
    if ( comm_rank_of_max==NULL ) {
        unsigned long int *send = x;
        unsigned long int *recv = new unsigned long int[n];
        MPI_Allreduce( send, recv, n, MPI_UNSIGNED_LONG, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         AMP_ERROR("Returning the rank of min with unsigned long int is not supported yet");
    }
}
// float
template <>
void AMP_MPI::call_maxReduce<float>(const float *x, float *y, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_FLOAT, MPI_MAX, communicator);
    } else {
         FloatIntStruct *recv = new FloatIntStruct[n];
         FloatIntStruct *send = new FloatIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].f = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_FLOAT_INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].f;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
template <>
void AMP_MPI::call_maxReduce<float>(float *x, const int n, int *comm_rank_of_max) const {
    if ( comm_rank_of_max==NULL ) {
        float *send = x;
        float *recv = new float[n];
        MPI_Allreduce( send, recv, n, MPI_FLOAT, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         FloatIntStruct *recv = new FloatIntStruct[n];
         FloatIntStruct *send = new FloatIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].f = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_FLOAT_INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].f;
            comm_rank_of_max[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
// double
template <>
void AMP_MPI::call_maxReduce<double>(const double *x, double *y, const int n, int *comm_rank_of_min) const {
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_DOUBLE, MPI_MAX, communicator);
    } else {
         DoubleIntStruct *recv = new DoubleIntStruct[n];
         DoubleIntStruct *send = new DoubleIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].d = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_DOUBLE_INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].d;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
template <>
void AMP_MPI::call_maxReduce<double>(double *x, const int n, int *comm_rank_of_max) const {
    if ( comm_rank_of_max==NULL ) {
        double *send = x;
        double *recv = new double[n];
        MPI_Allreduce( send, recv, n, MPI_DOUBLE, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         DoubleIntStruct *recv = new DoubleIntStruct[n];
         DoubleIntStruct *send = new DoubleIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].d = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_DOUBLE_INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].d;
            comm_rank_of_max[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
}
#endif


/************************************************************************
*  bcast                                                                *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
// char
template <>
void AMP_MPI::call_bcast<char>(char *x, const int n, const int root) const {
    MPI_Bcast( x, n, MPI_CHAR, root, communicator);
}
// int
template <>
void AMP_MPI::call_bcast<int>(int *x, const int n, const int root) const {
    MPI_Bcast( x, n, MPI_INT, root, communicator);
}
// float
template <>
void AMP_MPI::call_bcast<float>(float *x, const int n, const int root) const {
    MPI_Bcast( x, n, MPI_FLOAT, root, communicator);
}
// double
template <>
void AMP_MPI::call_bcast<double>(double *x, const int n, const int root) const {
    MPI_Bcast( x, n, MPI_DOUBLE, root, communicator);
}
#else
// We need a concrete instantiation of bcast<char>(x,n,root);
template <>
void AMP_MPI::call_bcast<char>(char *x, const int n, const int root) const {
    AMP_ERROR("Internal error in AMP_MPI (bcast) ");
}
#endif


/************************************************************************
*  Perform a global barrier across all processors.                      *
************************************************************************/
void AMP_MPI::barrier()
{
    #ifdef USE_EXT_MPI
        MPI_Barrier(communicator);
    #endif
}


/************************************************************************
*  Send data array to another processor.                                *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
// char
template <>
void AMP_MPI::send<char>(const char *buf, const int length, 
    const int recv_proc_number, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    // Send the data 
    MPI_Send((void*)buf, length, MPI_CHAR, recv_proc_number, tag, communicator);
}
// int
template <>
void AMP_MPI::send<int>(const int *buf, const int length, 
    const int recv_proc_number, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    // Send the data 
    MPI_Send((void*)buf, length, MPI_INT, recv_proc_number, tag, communicator);
}
// float
template <>
void AMP_MPI::send<float>(const float *buf, const int length, 
    const int recv_proc_number, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    // Send the data 
    MPI_Send((void*)buf, length, MPI_FLOAT, recv_proc_number, tag, communicator);
}
// double
template <>
void AMP_MPI::send<double>(const double *buf, const int length, 
    const int recv_proc_number, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    // Send the data 
    MPI_Send((void*)buf, length, MPI_DOUBLE, recv_proc_number, tag, communicator);
}
#else
// We need a concrete instantiation of send for USE_EXT_MPI=false
template <>
void AMP_MPI::send<char>(const char *buf, const int length, 
    const int recv_proc_number, int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request id = getRequest( communicator, tag );
    std::map<MPI_Request,Isendrecv_struct>::iterator it = global_isendrecv_list.find(id);
    AMP_INSIST(it==global_isendrecv_list.end(),"send must be paired with a previous call to irecv in serial");
    AMP_ASSERT(it->second.status==2);
    memcpy((char*)it->second.data,buf,length);
    global_isendrecv_list.erase( it );
}
#endif


/************************************************************************
*  Non-blocking send data array to another processor.                   *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
// char
template <>
MPI_Request AMP_MPI::Isend<char>(const char *buf, const int length, const int recv_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    MPI_Isend((void*)buf, length, MPI_CHAR, recv_proc, tag, communicator, &request);
    return request;
}
// int
template <>
MPI_Request AMP_MPI::Isend<int>(const int *buf, const int length, const int recv_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    MPI_Isend((void*)buf, length, MPI_INT, recv_proc, tag, communicator, &request);
    return request;
}
// float
template <>
MPI_Request AMP_MPI::Isend<float>(const float *buf, const int length, const int recv_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    MPI_Isend((void*)buf, length, MPI_FLOAT, recv_proc, tag, communicator, &request);
    return request;
}
// double
template <>
MPI_Request AMP_MPI::Isend<double>(const double *buf, const int length, const int recv_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    MPI_Isend((void*)buf, length, MPI_DOUBLE, recv_proc, tag, communicator, &request);
    return request;
}
#else
// We need a concrete instantiation of send for USE_EXT_MPI=false
template <>
MPI_Request AMP_MPI::Isend<char>(const char *buf, const int length, const int recv_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request id = getRequest( communicator, tag );
    std::map<MPI_Request,Isendrecv_struct>::iterator it = global_isendrecv_list.find(id);
    if ( it==global_isendrecv_list.end() ) {
        // We are calling isend first
        Isendrecv_struct data;
        data.data = buf;
        data.status = 1;
        global_isendrecv_list.insert( std::pair<MPI_Request,Isendrecv_struct>(id,data) );
    } else {
        // We called irecv first
        AMP_ASSERT(it->second.status==2);
        memcpy((char*)it->second.data,buf,length);
        global_isendrecv_list.erase( it );
    }
    return id;
}
#endif


/************************************************************************
*  Send byte array to another processor.                                *
************************************************************************/
void AMP_MPI::sendBytes(const void *buf, const int number_bytes, 
    const int recv_proc_number, int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    send<char>((const char*)buf,number_bytes,recv_proc_number,tag);
}


/************************************************************************
*  Non-blocking send byte array to another processor.                   *
************************************************************************/
MPI_Request AMP_MPI::IsendBytes(const void *buf, const int number_bytes, const int recv_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    return Isend<char>((const char*)buf,number_bytes,recv_proc,tag);
}


/************************************************************************
*  Recieve data array to another processor.                             *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
// char
template <>
void AMP_MPI::recv<char>(char *buf, int &length, 
    const int send_proc_number, const bool get_length, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    // Get the recieve length if necessary
    if (get_length) {
        int bytes = this->probe( send_proc_number, tag );
        int recv_length = bytes/sizeof(char);
        AMP_INSIST(length>=recv_length,"Recived length is larger than allocated array");
        length = recv_length;
    }
    // Send the data 
    MPI_Status status;
    MPI_Recv((void*)buf, length, MPI_CHAR, send_proc_number, tag, communicator, &status);
}
// int
template <>
void AMP_MPI::recv<int>(int *buf, int &length, 
    const int send_proc_number, const bool get_length, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    // Get the recieve length if necessary
    if (get_length) {
        int bytes = this->probe( send_proc_number, tag );
        int recv_length = bytes/sizeof(int);
        AMP_INSIST(length>=recv_length,"Recived length is larger than allocated array");
        length = recv_length;
    }
    // Send the data 
    MPI_Status status;
    MPI_Recv((void*)buf, length, MPI_INT, send_proc_number, tag, communicator, &status);
}
// float
template <>
void AMP_MPI::recv<float>(float *buf, int &length, 
    const int send_proc_number, const bool get_length, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    // Get the recieve length if necessary
    if (get_length) {
        int bytes = this->probe( send_proc_number, tag );
        int recv_length = bytes/sizeof(float);
        AMP_INSIST(length>=recv_length,"Recived length is larger than allocated array");
        length = recv_length;
    }
    // Send the data 
    MPI_Status status;
    MPI_Recv((void*)buf, length, MPI_FLOAT, send_proc_number, tag, communicator, &status);
}
// double
template <>
void AMP_MPI::recv<double>(double *buf, int &length, 
    const int send_proc_number, const bool get_length, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    // Get the recieve length if necessary
    if (get_length) {
        int bytes = this->probe( send_proc_number, tag );
        int recv_length = bytes/sizeof(double);
        AMP_INSIST(length>=recv_length,"Recived length is larger than allocated array");
        length = recv_length;
    }
    // Send the data 
    MPI_Status status;
    MPI_Recv((void*)buf, length, MPI_DOUBLE, send_proc_number, tag, communicator, &status);
}
#else
// We need a concrete instantiation of send for USE_EXT_MPI=false
template <>
void AMP_MPI::recv<char>(char *buf, int &length, 
    const int send_proc_number, const bool get_length, int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request id = getRequest( communicator, tag );
    std::map<MPI_Request,Isendrecv_struct>::iterator it = global_isendrecv_list.find(id);
    AMP_INSIST(it!=global_isendrecv_list.end(),"recv must be paired with a previous call to isend in serial");
    AMP_ASSERT(it->second.status==1);
    memcpy(buf,it->second.data,length);
    global_isendrecv_list.erase( it );
}
#endif


/************************************************************************
*  Non-blocking recieve data array to another processor.                *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
// char
template <>
MPI_Request AMP_MPI::Irecv<char>(char *buf, const int length, const int send_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    MPI_Irecv((void*)buf, length, MPI_CHAR, send_proc, tag, communicator, &request);
    return request;
}
// int
template <>
MPI_Request AMP_MPI::Irecv<int>(int *buf, const int length, const int send_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    MPI_Irecv((void*)buf, length, MPI_INT, send_proc, tag, communicator, &request);
    return request;
}
// float
template <>
MPI_Request AMP_MPI::Irecv<float>(float *buf, const int length, const int send_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    MPI_Irecv((void*)buf, length, MPI_FLOAT, send_proc, tag, communicator, &request);
    return request;
}
// double
template <>
MPI_Request AMP_MPI::Irecv<double>(double *buf, const int length, const int send_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    MPI_Irecv((void*)buf, length, MPI_DOUBLE, send_proc, tag, communicator, &request);
    return request;
}
#else
// We need a concrete instantiation of send for USE_EXT_MPI=false
template <>
MPI_Request AMP_MPI::Irecv<char>(char *buf, const int length, const int send_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request id = getRequest( communicator, tag );
    std::map<MPI_Request,Isendrecv_struct>::iterator it = global_isendrecv_list.find(id);
    if ( it==global_isendrecv_list.end() ) {
        // We are calling Irecv first
        Isendrecv_struct data;
        data.data = buf;
        data.status = 2;
        global_isendrecv_list.insert( std::pair<MPI_Request,Isendrecv_struct>(id,data) );
    } else {
        // We called Isend first
        AMP_ASSERT(it->second.status==1);
        memcpy(buf,it->second.data,length);
        global_isendrecv_list.erase( it );
    }
    return id;
}
#endif


/************************************************************************
*  Recieve byte array to another processor.                             *
************************************************************************/
void AMP_MPI::recvBytes(void *buf, int &number_bytes, const int send_proc, int tag) const
{
    recv<char>((char*)buf,number_bytes,send_proc,false,tag);
}


/************************************************************************
*  Recieve byte array to another processor.                             *
************************************************************************/
MPI_Request AMP_MPI::IrecvBytes(void *buf, const int number_bytes, const int send_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    return Irecv<char>((char*)buf,number_bytes,send_proc,tag);
}


/************************************************************************
*  allGather                                                            *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
// unsigned char
template <>
void AMP_MPI::call_allGather<unsigned char>(const unsigned char x_in, unsigned char *x_out) const {
    MPI_Allgather( (void*) &x_in, 1, MPI_UNSIGNED_CHAR, (void*) x_out, 1, MPI_UNSIGNED_CHAR, communicator );
}
template <>
void AMP_MPI::call_allGather<unsigned char>(const unsigned char *x_in, int size_in, unsigned char *x_out, int *size_out, int *disp_out) const {
    MPI_Allgatherv( (void*) x_in, size_in, MPI_CHAR, (void*) x_out, size_out, disp_out, MPI_CHAR, communicator );
}
// char
template <>
void AMP_MPI::call_allGather<char>(const char x_in, char *x_out) const {
    MPI_Allgather( (void*) &x_in, 1, MPI_CHAR, (void*) x_out, 1, MPI_CHAR, communicator );
}
template <>
void AMP_MPI::call_allGather<char>(const char *x_in, int size_in, char *x_out, int *size_out, int *disp_out) const {
    MPI_Allgatherv( (void*) x_in, size_in, MPI_CHAR, (void*) x_out, size_out, disp_out, MPI_CHAR, communicator );
}
// unsigned int
template <>
void AMP_MPI::call_allGather<unsigned int>(const unsigned int x_in, unsigned int *x_out) const {
    MPI_Allgather( (void*) &x_in, 1, MPI_UNSIGNED, (void*) x_out, 1, MPI_UNSIGNED, communicator );
}
template <>
void AMP_MPI::call_allGather<unsigned int>(const unsigned int *x_in, int size_in, unsigned int *x_out, int *size_out, int *disp_out) const {
    MPI_Allgatherv( (void*) x_in, size_in, MPI_UNSIGNED, (void*) x_out, size_out, disp_out, MPI_UNSIGNED, communicator );
}
// int
template <>
void AMP_MPI::call_allGather<int>(const int x_in, int *x_out) const {
    MPI_Allgather( (void*) &x_in, 1, MPI_INT, (void*) x_out, 1, MPI_INT, communicator );
}
template <>
void AMP_MPI::call_allGather<int>(const int *x_in, int size_in, int *x_out, int *size_out, int *disp_out) const {
    MPI_Allgatherv( (void*) x_in, size_in, MPI_INT, (void*) x_out, size_out, disp_out, MPI_INT, communicator );
}
// unsigned long int
template <>
void AMP_MPI::call_allGather<unsigned long int>(const unsigned long int x_in, unsigned long int *x_out) const {
    MPI_Allgather( (void*) &x_in, 1, MPI_UNSIGNED_LONG, (void*) x_out, 1, MPI_UNSIGNED_LONG, communicator );
}
template <>
void AMP_MPI::call_allGather<unsigned long int>(const unsigned long int *x_in, int size_in, unsigned long int *x_out, int *size_out, int *disp_out) const {
    MPI_Allgatherv( (void*) x_in, size_in, MPI_UNSIGNED_LONG, (void*) x_out, size_out, disp_out, MPI_UNSIGNED_LONG, communicator );
}
// long int
template <>
void AMP_MPI::call_allGather<long int>(const long int x_in, long int *x_out) const {
    MPI_Allgather( (void*) &x_in, 1, MPI_LONG, (void*) x_out, 1, MPI_LONG, communicator );
}
template <>
void AMP_MPI::call_allGather<long int>(const long int *x_in, int size_in, long int *x_out, int *size_out, int *disp_out) const {
    MPI_Allgatherv( (void*) x_in, size_in, MPI_LONG, (void*) x_out, size_out, disp_out, MPI_LONG, communicator );
}
// float
template <>
void AMP_MPI::call_allGather<float>(const float x_in, float *x_out) const {
    MPI_Allgather( (void*) &x_in, 1, MPI_FLOAT, (void*) x_out, 1, MPI_FLOAT, communicator );
}
template <>
void AMP_MPI::call_allGather<float>(const float *x_in, int size_in, float *x_out, int *size_out, int *disp_out) const {
    MPI_Allgatherv( (void*) x_in, size_in, MPI_FLOAT, (void*) x_out, size_out, disp_out, MPI_FLOAT, communicator );
}
// double
template <>
void AMP_MPI::call_allGather<double>(const double x_in, double *x_out) const {
    MPI_Allgather( (void*) &x_in, 1, MPI_DOUBLE, (void*) x_out, 1, MPI_DOUBLE, communicator );
}
template <>
void AMP_MPI::call_allGather<double>(const double *x_in, int size_in, double *x_out, int *size_out, int *disp_out) const {
    MPI_Allgatherv( (void*) x_in, size_in, MPI_DOUBLE, (void*) x_out, size_out, disp_out, MPI_DOUBLE, communicator );
}
#else
// We need a concrete instantiation of call_allGather<char>(x_in,size_in,x_out,size_out)
template <>
void AMP_MPI::call_allGather<char>(const char *x_in, int size_in, char *x_out, int *size_out, int *disp_out) const {
    AMP_ERROR("Internal error in AMP_MPI (allGather) ");
}
#endif


/************************************************************************
*  allToAll                                                             *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
template <> void AMP_MPI::allToAll<unsigned char>(const int n, const unsigned char *send, unsigned char *recv ) const {
    MPI_Alltoall( (void*) send, n, MPI_UNSIGNED_CHAR, (void*) recv, n, MPI_UNSIGNED_CHAR, communicator);
}
template <> void AMP_MPI::allToAll<char>(const int n, const char *send, char *recv ) const {
    MPI_Alltoall( (void*) send, n, MPI_CHAR, (void*) recv, n, MPI_CHAR, communicator);
}
template <> void AMP_MPI::allToAll<unsigned int>(const int n, const unsigned int *send, unsigned int *recv ) const {
    MPI_Alltoall( (void*) send, n, MPI_UNSIGNED, (void*) recv, n, MPI_UNSIGNED, communicator);
}
template <> void AMP_MPI::allToAll<int>(const int n, const int *send, int *recv ) const {
    MPI_Alltoall( (void*) send, n, MPI_INT, (void*) recv, n, MPI_INT, communicator);
}
template <> void AMP_MPI::allToAll<unsigned long int>(const int n, const unsigned long int *send, unsigned long int *recv ) const {
    MPI_Alltoall( (void*) send, n, MPI_UNSIGNED_LONG, (void*) recv, n, MPI_UNSIGNED_LONG, communicator);
}
template <> void AMP_MPI::allToAll<long int>(const int n, const long int *send, long int *recv ) const {
    MPI_Alltoall( (void*) send, n, MPI_LONG, (void*) recv, n, MPI_LONG, communicator);
}
template <> void AMP_MPI::allToAll<float>(const int n, const float *send, float *recv ) const {
    MPI_Alltoall( (void*) send, n, MPI_FLOAT, (void*) recv, n, MPI_FLOAT, communicator);
}
template <> void AMP_MPI::allToAll<double>(const int n, const double *send, double *recv ) const {
    MPI_Alltoall( (void*) send, n, MPI_DOUBLE, (void*) recv, n, MPI_DOUBLE, communicator);
}
#endif


/************************************************************************
*  call_allToAll                                                        *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
// unsigned char
template <>
void AMP_MPI::call_allToAll<unsigned char>(const unsigned char *send_data, const int send_cnt[], 
        const int send_disp[], unsigned char *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_UNSIGNED_CHAR, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_UNSIGNED_CHAR, communicator ); 
}
// char
template <>
void AMP_MPI::call_allToAll<char>(const char *send_data, const int send_cnt[], 
        const int send_disp[], char *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_CHAR, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_CHAR, communicator ); 
}
// unsigned int
template <>
void AMP_MPI::call_allToAll<unsigned int>(const unsigned int *send_data, const int send_cnt[], 
        const int send_disp[], unsigned int *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_UNSIGNED, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_UNSIGNED, communicator ); 
}
// int
template <>
void AMP_MPI::call_allToAll<int>(const int *send_data, const int send_cnt[], 
        const int send_disp[], int *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_INT, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_INT, communicator ); 
}
// unsigned long int
template <>
void AMP_MPI::call_allToAll<unsigned long int>(const unsigned long int *send_data, const int send_cnt[], 
        const int send_disp[], unsigned long int *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_UNSIGNED_LONG, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_UNSIGNED_LONG, communicator ); 
}
// long int
template <>
void AMP_MPI::call_allToAll<long int>(const long int *send_data, const int send_cnt[], 
        const int send_disp[], long int *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_LONG, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_LONG, communicator ); 
}
// float
template <>
void AMP_MPI::call_allToAll<float>(const float *send_data, const int send_cnt[], 
        const int send_disp[], float *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_FLOAT, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_FLOAT, communicator ); 
}
// double
template <>
void AMP_MPI::call_allToAll<double>(const double *send_data, const int send_cnt[], 
        const int send_disp[], double *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_DOUBLE, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_DOUBLE, communicator ); 
}
#else
// Default instatiation of unsigned char
template <>
void AMP_MPI::call_allToAll<char>(const char *send_data, const int send_cnt[], 
        const int send_disp[], char *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    AMP_ERROR("Should not reach this point");
}
#endif


/************************************************************************
*  call_sumScan                                                         *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
// unsigned char
template <>
void AMP_MPI::call_sumScan<unsigned char>(const unsigned char *send, unsigned char *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED_CHAR, MPI_SUM, communicator);
}
// char
template <>
void AMP_MPI::call_sumScan<char>(const char *send, char *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_SIGNED_CHAR, MPI_SUM, communicator);
}
// unsigned int
template <>
void AMP_MPI::call_sumScan<unsigned int>(const unsigned int *send, unsigned int *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED, MPI_SUM, communicator);
}
// int
template <>
void AMP_MPI::call_sumScan<int>(const int *send, int *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_INT, MPI_SUM, communicator);
}
// long int
template <>
void AMP_MPI::call_sumScan<long int>(const long int *send, long int *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_LONG, MPI_SUM, communicator);
}
// unsigned long int
template <>
void AMP_MPI::call_sumScan<unsigned long>(const unsigned long *send, unsigned long *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED_LONG, MPI_SUM, communicator);
}
// float
template <>
void AMP_MPI::call_sumScan<float>(const float *send, float *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_FLOAT, MPI_SUM, communicator);
}
// double
template <>
void AMP_MPI::call_sumScan<double>(const double *send, double *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_DOUBLE, MPI_SUM, communicator);
}
// std::complex<double>
template <>
void AMP_MPI::call_sumScan< std::complex<double> >(const std::complex<double> *x, std::complex<double> *y, int n) const {
    double *send = new double[2*n];
    double *recv = new double[2*n];
    for (int i=0; i<n; i++) {
        send[2*i+0] = real(x[i]);
        send[2*i+1] = imag(x[i]);
    }
    MPI_Scan( (void*) send, (void*) recv, 2*n, MPI_DOUBLE, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        y[i] = std::complex<double>(recv[2*i+0],recv[2*i+1]);
    delete [] send;
    delete [] recv;
}
#endif


/************************************************************************
*  call_minScan                                                         *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
// unsigned char
template <>
void AMP_MPI::call_minScan<unsigned char>(const unsigned char *send, unsigned char *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED_CHAR, MPI_MIN, communicator);
}
// char
template <>
void AMP_MPI::call_minScan<char>(const char *send, char *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_SIGNED_CHAR, MPI_MIN, communicator);
}
// unsigned int
template <>
void AMP_MPI::call_minScan<unsigned int>(const unsigned int *send, unsigned int *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED, MPI_MIN, communicator);
}
// int
template <>
void AMP_MPI::call_minScan<int>(const int *send, int *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_INT, MPI_MIN, communicator);
}
// unsigned long int
template <>
void AMP_MPI::call_minScan<unsigned long int>(const unsigned long int *send, unsigned long int *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED_LONG, MPI_MIN, communicator);
}
// long int
template <>
void AMP_MPI::call_minScan<long int>(const long int *send, long int *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_LONG, MPI_MIN, communicator);
}
// float
template <>
void AMP_MPI::call_minScan<float>(const float *send, float *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_FLOAT, MPI_MIN, communicator);
}
// double
template <>
void AMP_MPI::call_minScan<double>(const double *send, double *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_DOUBLE, MPI_MIN, communicator);
}
#endif


/************************************************************************
*  call_maxScan                                                         *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_EXT_MPI
// unsigned char
template <>
void AMP_MPI::call_maxScan<unsigned char>(const unsigned char *send, unsigned char *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED_CHAR, MPI_MAX, communicator);
}
// char
template <>
void AMP_MPI::call_maxScan<char>(const char *send, char *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_SIGNED_CHAR, MPI_MAX, communicator);
}
// unsigned int
template <>
void AMP_MPI::call_maxScan<unsigned int>(const unsigned int *send, unsigned int *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED, MPI_MAX, communicator);
}
// int
template <>
void AMP_MPI::call_maxScan<int>(const int *send, int *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_INT, MPI_MAX, communicator);
}
// long int
template <>
void AMP_MPI::call_maxScan<long int>(const long int *send, long int *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_LONG, MPI_MAX, communicator);
}
// unsigned long int
template <>
void AMP_MPI::call_maxScan<unsigned long int>(const unsigned long int *send, unsigned long int *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED_LONG, MPI_MAX, communicator);
}
// float
template <>
void AMP_MPI::call_maxScan<float>(const float *send, float *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_INT, MPI_MAX, communicator);
}
// double
template <>
void AMP_MPI::call_maxScan<double>(const double *send, double *recv, int n) const {
    MPI_Scan( (void*) send, (void*) recv, n, MPI_DOUBLE, MPI_MAX, communicator);
}
#endif


/************************************************************************
*  Wait functions                                                       *
************************************************************************/
#ifdef USE_EXT_MPI
void AMP_MPI::wait( MPI_Request request) {
    MPI_Status  status;
    int flag = 0;
    int err = MPI_Test( &request, &flag, &status );
    AMP_ASSERT(err==MPI_SUCCESS);   // Check that the first call is valid
    while ( !flag ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Test( &request, &flag, &status );
    }
}
int AMP_MPI::waitAny( int count, MPI_Request *request) {
    if ( count==0 ) 
        return -1;
    int index = -1;
    int flag = 0;
    MPI_Status *status = new MPI_Status[count];
    int err = MPI_Testany( count, request, &index, &flag, status );
    AMP_ASSERT(err==MPI_SUCCESS);   // Check that the first call is valid
    while ( !flag ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Testany( count, request, &index, &flag, status );
    }
    AMP_ASSERT(index>=0);   // Check that the index is valid
    delete [] status;
    return index;
}
void AMP_MPI::waitAll( int count, MPI_Request *request) {
    if ( count==0 ) 
        return;
    int flag = 0;
    MPI_Status *status = new MPI_Status[count];
    int err = MPI_Testall( count, request, &flag, status );
    AMP_ASSERT(err==MPI_SUCCESS);   // Check that the first call is valid
    while ( !flag ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Testall( count, request, &flag, status );
    }
    delete [] status;
}
#else
void AMP_MPI::wait( MPI_Request request) {
    while ( 1 ) {
        // Check if the request is in our list
        if ( global_isendrecv_list.find(request)==global_isendrecv_list.end() )
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
}
int AMP_MPI::waitAny( int count, MPI_Request *request) {
    int index = 0;
    while ( 1 ) {
        // Check if the request is in our list
        bool found_any = false;
        for (int i=0; i<count; i++) {
            if ( global_isendrecv_list.find(request[i])==global_isendrecv_list.end() ) {
                found_any = true;
                index = i;
            }
        }
        if ( found_any )
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
    return index;
}
void AMP_MPI::waitAll( int count, MPI_Request *request) {
    while ( 1 ) {
        // Check if the request is in our list
        bool found_all = true;
        for (int i=0; i<count; i++) {
            if ( global_isendrecv_list.find(request[i])!=global_isendrecv_list.end() )
                found_all = false;
        }
        if ( found_all )
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
}
#endif


/************************************************************************
*  Probe functions                                                      *
************************************************************************/
#ifdef USE_EXT_MPI
int AMP_MPI::Iprobe( int source, int tag) const {
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Status status;
    int flag = 0;
    MPI_Iprobe(source,tag,communicator,&flag,&status);
    if ( flag==0 )
        return -1;
    int count;
    MPI_Get_count(&status,MPI_BYTE,&count);
    AMP_ASSERT(count>=0);
    return count;
}
int AMP_MPI::probe( int source, int tag) const {
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Status status;
    MPI_Probe(source,tag,communicator,&status);
    int count;
    MPI_Get_count(&status,MPI_BYTE,&count);
    AMP_ASSERT(count>=0);
    return count;
}
#else
int AMP_MPI::Iprobe( int source, int tag) const {
    AMP_ERROR("Not implimented for serial codes (Iprobe)");
    return 0;
}
int AMP_MPI::probe( int source, int tag) const {
    AMP_ERROR("Not implimented for serial codes (probe)");
    return 0;
}
#endif



/************************************************************************
*  Timer functions                                                      *
************************************************************************/
#ifdef USE_EXT_MPI
    double AMP_MPI::time() { 
        return MPI_Wtime();
    }
    double AMP_MPI::tick() {
        return MPI_Wtick();
    }
#else
    #if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        double AMP_MPI::time() { 
            LARGE_INTEGER end, f;
            QueryPerformanceFrequency(&f);
            QueryPerformanceCounter(&end);       
            double time = ((double)end.QuadPart)/((double)f.QuadPart);
            return time;
        }
        double AMP_MPI::tick() {
            LARGE_INTEGER f;
            QueryPerformanceFrequency(&f);
            double resolution = ((double)1.0)/((double)f.QuadPart);
            return resolution;
        }
    #else
        double AMP_MPI::time() { 
            timeval current_time;
            gettimeofday(&current_time,NULL);
            double time = ((double)current_time.tv_sec)+1e-6*((double)current_time.tv_usec);
            return time;
        }
        double AMP_MPI::tick() {
            timeval start, end;
            gettimeofday(&start,NULL);
            gettimeofday(&end,NULL);
            while ( end.tv_sec==start.tv_sec &&  end.tv_usec==start.tv_usec )
                gettimeofday(&end,NULL);
            double resolution = ((double)(end.tv_sec-start.tv_sec))+1e-6*((double)(end.tv_usec-start.tv_usec));
            return resolution;
        }
    #endif
#endif


} // namespace AMP

