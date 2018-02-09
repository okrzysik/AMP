#include "AMP/utils/StackTrace.h"

#include <algorithm>
#include <csignal>
#include <cstring>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>


// Detect the OS
// clang-format off
#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 ) || defined( _MSC_VER )
    #define USE_WINDOWS
    #define NOMINMAX
#elif defined( __APPLE__ )
    #define USE_MAC
    #define USE_NM
#elif defined( __linux ) || defined( __unix ) || defined( __posix )
    #define USE_LINUX
    #define USE_NM
#else
    #error Unknown OS
#endif
// clang-format on


// Include system dependent headers
// clang-format off
// Detect the OS and include system dependent headers
#ifdef USE_WINDOWS
    #include <windows.h>
    #include <dbghelp.h>
    #include <DbgHelp.h>
    #include <TlHelp32.h>
    #include <Psapi.h>
    #include <process.h>
    #include <stdio.h>
    #include <tchar.h>
    #pragma comment( lib, "version.lib" ) // for "VerQueryValue"
#else
    #include <dlfcn.h>
    #include <execinfo.h>
    #include <sched.h>
    #include <sys/time.h>
    #include <ctime>
    #include <unistd.h>
    #include <sys/syscall.h>
#endif
#ifdef USE_MAC
    #include <mach-o/dyld.h>
    #include <mach/mach.h>
    #include <sys/sysctl.h>
    #include <sys/types.h>
#endif
// clang-format on


#ifdef __GNUC__
#define USE_ABI
#include <cxxabi.h>
#endif


// Check for and include MPI
// clang-format off
#if defined(USE_MPI)
    #include "mpi.h"
#elif defined(__has_include)
    #if __has_include("mpi.h")
        #include "mpi.h"
    #endif
    #define USE_MPI
#elif defined(USE_EXT_MPI)
    #include "mpi.h"
    #define USE_MPI
#endif
// clang-format on


#ifndef NULL_USE
#define NULL_USE( variable )                 \
    do {                                     \
        if ( 0 ) {                           \
            char *temp = (char *) &variable; \
            temp++;                          \
        }                                    \
    } while ( 0 )
#endif


namespace AMP {


// Set the callstack signal
#ifdef SIGRTMIN
#define CALLSTACK_SIG SIGRTMIN + 4
#else
#define CALLSTACK_SIG SIGUSR1
#define SIGRTMIN SIGUSR1
#define SIGRTMAX SIGUSR1
#endif


// Utility to call system command and return output
#if defined( USE_LINUX ) || defined( USE_MAC )
static inline std::string exec( const std::string &cmd )
{
    signal( SIGCHLD, SIG_DFL ); // Clear child exited
    FILE *pipe = popen( cmd.c_str(), "r" );
    if ( pipe == nullptr )
        return std::string();
    std::string result = "";
    result.reserve( 1024 );
    while ( !feof( pipe ) ) {
        char buffer[257];
        buffer[256] = 0;
        if ( fgets( buffer, 128, pipe ) != nullptr )
            result += buffer;
    }
    pclose( pipe );
    return result;
}
#endif


// Utility to break a string by a newline
static inline std::vector<std::string> breakString( const std::string &str )
{
    std::vector<std::string> strvec;
    size_t i1 = 0;
    size_t i2 = std::min( str.find( '\n', i1 ), str.length() );
    while ( i1 < str.length() ) {
        strvec.push_back( str.substr( i1, i2 - i1 ) );
        i1 = i2 + 1;
        i2 = std::min( str.find( '\n', i1 ), str.length() );
    }
    return strvec;
}


// Utility to strip the path from a filename
static inline std::string stripPath( const std::string &filename )
{
    if ( filename.empty() ) {
        return std::string();
    }
    int i = 0;
    for ( i = (int) filename.size() - 1; i >= 0 && filename[i] != 47 && filename[i] != 92; i-- ) {
    }
    i = std::max( 0, i + 1 );
    return filename.substr( i );
}


// Inline function to subtract two addresses returning the absolute difference
static inline void *subtractAddress( void *a, void *b )
{
    return reinterpret_cast<void *>(
        std::abs( reinterpret_cast<long long int>( a ) - reinterpret_cast<long long int>( b ) ) );
}


#ifdef USE_WINDOWS
static BOOL __stdcall readProcMem( HANDLE hProcess,
                                   DWORD64 qwBaseAddress,
                                   PVOID lpBuffer,
                                   DWORD nSize,
                                   LPDWORD lpNumberOfBytesRead )
{
    SIZE_T st;
    BOOL bRet = ReadProcessMemory( hProcess, (LPVOID) qwBaseAddress, lpBuffer, nSize, &st );
    *lpNumberOfBytesRead = (DWORD) st;
    return bRet;
}
static inline std::string getCurrentDirectory()
{
    char temp[1024] = { 0 };
    GetCurrentDirectoryA( sizeof( temp ), temp );
    return temp;
}
namespace StackTrace {
BOOL GetModuleListTH32( HANDLE hProcess, DWORD pid );
BOOL GetModuleListPSAPI( HANDLE hProcess );
DWORD LoadModule( HANDLE hProcess, LPCSTR img, LPCSTR mod, DWORD64 baseAddr, DWORD size );
void LoadModules();
}; // namespace StackTrace
#endif


/****************************************************************************
 *  stack_info                                                               *
 ****************************************************************************/
std::string StackTrace::stack_info::print() const
{
    char tmp[32];
    sprintf( tmp, "0x%016llx:  ", reinterpret_cast<unsigned long long int>( address ) );
    std::string stack( tmp );
    sprintf( tmp, "%i", line );
    std::string line_str( tmp );
    stack += stripPath( object );
    stack.resize( std::max<size_t>( stack.size(), 38 ), ' ' );
    stack += "  " + function;
    if ( !filename.empty() && line > 0 ) {
        stack.resize( std::max<size_t>( stack.size(), 72 ), ' ' );
        stack += "  " + stripPath( filename ) + ":" + line_str;
    } else if ( !filename.empty() ) {
        stack.resize( std::max<size_t>( stack.size(), 72 ), ' ' );
        stack += "  " + stripPath( filename );
    } else if ( line > 0 ) {
        stack += " : " + line_str;
    }
    return stack;
}
/*static int maxDepth( const StackTrace::multi_stack_info& stack )
{
    int depth = 0;
    for ( auto child : stack.children )
        depth = std::max<int>( depth, maxDepth( child ) );
    return depth+1;
}*/
std::vector<std::string> StackTrace::multi_stack_info::print( const std::string &prefix ) const
{
    std::vector<std::string> text;
    // auto depth = maxDepth( *this );
    // std::string line = prefix + "[" + std::to_string( N ) + "] ";
    // for (auto i=1; i<depth; i++)
    //    line += "--";
    // line += stack.print();
    std::string line = prefix + "[" + std::to_string( N ) + "] " + stack.print();
    text.push_back( line );
    std::string prefix2 = prefix + "  ";
    for ( size_t i = 0; i < children.size(); i++ ) {
        const auto &child = children[i];
        auto tmp          = child.print();
        for ( size_t j = 0; j < tmp.size(); j++ ) {
            std::string line = prefix2 + tmp[j];
            if ( children.size() > 1 && j > 0 && i < children.size() - 1 )
                line[prefix2.size()] = '|';
            text.push_back( line );
        }
    }
    return text;
}


/****************************************************************************
 *  Function to find an entry                                                *
 ****************************************************************************/
template<class TYPE>
inline size_t findfirst( const std::vector<TYPE> &X, TYPE Y )
{
    if ( X.empty() )
        return 0;
    size_t lower = 0;
    size_t upper = X.size() - 1;
    if ( X[lower] >= Y )
        return lower;
    if ( X[upper] < Y )
        return upper;
    while ( ( upper - lower ) != 1 ) {
        size_t value = ( upper + lower ) / 2;
        if ( X[value] >= Y )
            upper = value;
        else
            lower = value;
    }
    return upper;
}


/****************************************************************************
 * Function to get symbols for the executable from nm (if availible)         *
 * Note: this function maintains an internal cached copy to prevent          *
 *    exccessive calls to nm.  This function also uses a lock to ensure      *
 *    thread safety.                                                         *
 ****************************************************************************/
std::mutex getSymbols_mutex;
struct global_symbols_struct {
    std::vector<void *> address;
    std::vector<char> type;
    std::vector<std::string> obj;
    int error;
} global_symbols;
std::string StackTrace::getExecutable()
{
    std::string exe;
    try {
#ifdef USE_LINUX
        auto *buf = new char[0x10000];
        int len   = ::readlink( "/proc/self/exe", buf, 0x10000 );
        if ( len != -1 ) {
            buf[len] = '\0';
            exe      = std::string( buf );
        }
        delete[] buf;
#elif defined( USE_MAC )
        uint32_t size = 0x10000;
        char *buf     = new char[size];
        memset( buf, 0, size );
        if ( _NSGetExecutablePath( buf, &size ) == 0 )
            exe = std::string( buf );
        delete[] buf;
#elif defined( USE_WINDOWS )
        DWORD size = 0x10000;
        char *buf  = new char[size];
        memset( buf, 0, size );
        GetModuleFileName( nullptr, buf, size );
        exe = std::string( buf );
        delete[] buf;
#endif
    } catch ( ... ) {
    }
    return exe;
}
std::string global_exe_name = StackTrace::getExecutable();
static const global_symbols_struct &getSymbols2()
{
    static bool loaded = false;
    static global_symbols_struct data;
    // Load the symbol tables if they have not been loaded
    if ( !loaded ) {
        getSymbols_mutex.lock();
        if ( !loaded ) {
            loaded = true;
#ifdef USE_NM
            try {
                char cmd[1024];
#ifdef USE_LINUX
                sprintf( cmd, "nm -n --demangle %s", global_exe_name.c_str() );
#elif defined( USE_MAC )
                sprintf( cmd, "nm -n %s | c++filt", global_exe_name.c_str() );
#else
#error Unknown OS using nm
#endif
                auto output = breakString( exec( cmd ) );
                for ( const auto &line : output ) {
                    if ( line.empty() )
                        continue;
                    if ( line[0] == ' ' )
                        continue;
                    auto *a = const_cast<char *>( line.c_str() );
                    char *b = strchr( a, ' ' );
                    if ( b == nullptr )
                        continue;
                    b[0] = 0;
                    b++;
                    char *c = strchr( b, ' ' );
                    if ( c == nullptr )
                        continue;
                    c[0] = 0;
                    c++;
                    char *d = strchr( c, '\n' );
                    if ( d )
                        d[0] = 0;
                    size_t add = strtoul( a, nullptr, 16 );
                    data.address.push_back( reinterpret_cast<void *>( add ) );
                    data.type.push_back( b[0] );
                    data.obj.emplace_back( c );
                }
            } catch ( ... ) {
                data.error = -3;
            }
            data.error = 0;
#else
            data.error = -1;
#endif
        }
        getSymbols_mutex.unlock();
    }
    return data;
}
int StackTrace::getSymbols( std::vector<void *> &address,
                            std::vector<char> &type,
                            std::vector<std::string> &obj )
{
    const global_symbols_struct &data = getSymbols2();
    address                           = data.address;
    type                              = data.type;
    obj                               = data.obj;
    return data.error;
}


/****************************************************************************
 *  Function to get call stack info                                          *
 ****************************************************************************/
#ifdef USE_MAC
static void *loadAddress( const std::string &object )
{
    static std::map<std::string, void *> obj_map;
    if ( obj_map.empty() ) {
        uint32_t numImages = _dyld_image_count();
        for ( uint32_t i = 0; i < numImages; i++ ) {
            const struct mach_header *header = _dyld_get_image_header( i );
            const char *name                 = _dyld_get_image_name( i );
            const char *p                    = strrchr( name, '/' );
            struct mach_header *address      = const_cast<struct mach_header *>( header );
            obj_map.insert( std::pair<std::string, void *>( p + 1, address ) );
            // printf("   module=%s, address=%p\n", p + 1, header);
        }
    }
    auto it       = obj_map.find( object );
    void *address = 0;
    if ( it != obj_map.end() ) {
        address = it->second;
    } else {
        it = obj_map.find( stripPath( object ) );
        if ( it != obj_map.end() )
            address = it->second;
    }
    // printf("%s: 0x%016llx\n",object.c_str(),address);
    return address;
}
static std::tuple<std::string, std::string, std::string, int> split_atos( const std::string &buf )
{
    if ( buf.empty() )
        return std::tuple<std::string, std::string, std::string, int>();
    // Get the function
    size_t index = buf.find( " (in " );
    if ( index == std::string::npos )
        return std::make_tuple(
            buf.substr( 0, buf.length() - 1 ), std::string(), std::string(), 0 );
    std::string fun = buf.substr( 0, index );
    std::string tmp = buf.substr( index + 5 );
    // Get the object
    index           = tmp.find( ')' );
    std::string obj = tmp.substr( 0, index );
    tmp             = tmp.substr( index + 1 );
    // Get the filename and line number
    size_t p1 = tmp.find( '(' );
    size_t p2 = tmp.find( ')' );
    tmp       = tmp.substr( p1 + 1, p2 - p1 - 1 );
    index     = tmp.find( ':' );
    std::string file;
    int line = 0;
    if ( index != std::string::npos ) {
        file = tmp.substr( 0, index );
        line = std::stoi( tmp.substr( index + 1 ) );
    } else if ( p1 != std::string::npos ) {
        file = tmp;
    }
    return std::make_tuple( fun, obj, file, line );
}
#endif
#ifdef USE_LINUX
using uint_p = uint64_t;
#elif defined( USE_MAC )
typedef unsigned long uint_p;
#endif
#if defined( USE_LINUX ) || defined( USE_MAC )
static inline std::string generateCmd( const std::string &s1,
                                       const std::string &s2,
                                       const std::string &s3,
                                       std::vector<void *> addresses,
                                       const std::string &s4 )
{
    std::string cmd = s1 + s2 + s3;
    for ( auto &addresse : addresses ) {
        char tmp[32];
        sprintf( tmp, "%lx ", reinterpret_cast<uint_p>( addresse ) );
        cmd += tmp;
    }
    cmd += s4;
    return cmd;
}
#endif
// clang-format off
static void getFileAndLineObject( std::vector<StackTrace::stack_info*> &info )
{
    if ( info.empty() )
        return;
    // This gets the file and line numbers for multiple stack lines in the same object
    #if defined( USE_LINUX )
        // Create the call command
        std::vector<void*> address_list(info.size(),nullptr);
        for (size_t i=0; i<info.size(); i++) {
            address_list[i] = info[i]->address;
            if ( info[i]->object.find( ".so" ) != std::string::npos )
                address_list[i] = info[i]->address2; 
        }
        std::string cmd = generateCmd( "addr2line -C -e ", info[0]->object,
            " -f -i ", address_list, " 2> /dev/null" );
        // Get the function/line/file
        auto cmd_output = exec( cmd );
        auto output = breakString( cmd_output );
        if ( output.size() != 2*info.size() )
            return;
        // Add the results to info
        for (size_t i=0; i<info.size(); i++) {
            // get function name
            if ( info[i]->function.empty() ) {
                info[i]->function = output[2*i+0];
                info[i]->function.resize( std::max<size_t>( info[i]->function.size(), 1 ) - 1 );
            }
            // get file and line
            const char *buf = output[2*i+1].c_str();
            if ( buf[0] != '?' && buf[0] != 0 ) {
                size_t j = 0;
                for ( j = 0; j < 4095 && buf[j] != ':'; j++ ) {
                }
                info[i]->filename = std::string( buf, j );
                info[i]->line     = atoi( &buf[j + 1] );
            }
        }
    #elif defined( USE_MAC )
        // Create the call command
        void* load_address = loadAddress( info[0]->object );
        if ( load_address == nullptr )
            return;
        std::vector<void*> address_list(info.size(),nullptr);
        for (size_t i=0; i<info.size(); i++)
            address_list[i] = info[i]->address;
        // Call atos to get the object info
        char tmp[64];
        sprintf( tmp, " -l %lx ", (uint_p) load_address );
        std::string cmd = generateCmd( "atos -o ", info[0]->object,
            tmp, address_list, " 2> /dev/null" );
        // Get the function/line/file
        auto cmd_output = exec( cmd );
        auto output = breakString( cmd_output );
        if ( output.size() != info.size() )
            return;
        // Parse the output for function, file and line info
        for ( size_t i=0; i<info.size(); i++) {
            auto data = split_atos( output[i] );
            if ( info[i]->function.empty() )
                info[i]->function = std::get<0>(data);
            if ( info[i]->object.empty() )
                info[i]->object = std::get<1>(data);
            if ( info[i]->filename.empty() )
                info[i]->filename = std::get<2>(data);
            if ( info[i]->line==0 )
                info[i]->line = std::get<3>(data);
        }
    #endif
}
static void getFileAndLine( std::vector<StackTrace::stack_info> &info )
{
    // Build a list of stack elements for each object
    std::map<std::string,std::vector<StackTrace::stack_info*>> obj_map;
    for (auto & i : info) {
        auto& list = obj_map[i.object];
        list.emplace_back( &i );
    }
    // For each object, get the file/line numbers for all entries
    for ( auto& entry : obj_map ) 
        getFileAndLineObject( entry.second );
}
// Try to use the global symbols to decode info about the stack
static void getDataFromGlobalSymbols( StackTrace::stack_info &info )
{
    const global_symbols_struct &data = getSymbols2();
    if ( data.error == 0 ) {
        size_t index = findfirst( global_symbols.address, info.address );
        if ( index > 0 )
            info.object = global_symbols.obj[index - 1];
        else
            info.object = global_exe_name;
    }
}
static void signal_handler( int sig )
{
    printf("Signal caught acquiring stack (%i)\n",sig);
    StackTrace::setErrorHandlers( [](std::string,StackTrace::terminateType) { exit( -1 ); } );
}
StackTrace::stack_info StackTrace::getStackInfo( void *address )
{
    return getStackInfo( std::vector<void*>(1,address) )[0];
}
std::vector<StackTrace::stack_info> StackTrace::getStackInfo( const std::vector<void*>& address )
{
    // Temporarily handle signals to prevent recursion on the stack
    auto prev_handler = signal( SIGINT, signal_handler );
    // Get the detailed stack info
    std::vector<StackTrace::stack_info> info(address.size());
    try {
        #ifdef USE_WINDOWS
            IMAGEHLP_SYMBOL64 pSym[1024];
            memset( pSym, 0, sizeof( pSym ) );
            pSym->SizeOfStruct  = sizeof( IMAGEHLP_SYMBOL64 );
            pSym->MaxNameLength = 1024;

            IMAGEHLP_MODULE64 Module;
            memset( &Module, 0, sizeof( Module ) );
            Module.SizeOfStruct = sizeof( Module );

            HANDLE pid = GetCurrentProcess();

            for (size_t i=0; i<address.size(); i++) {
                info[i].address = address[i];
                DWORD64 address2 = reinterpret_cast<DWORD64>( address[i] );
                DWORD64 offsetFromSmybol;
                if ( SymGetSymFromAddr( pid, address2, &offsetFromSmybol, pSym ) != FALSE ) {
                    char name[8192]={0};
                    DWORD rtn = UnDecorateSymbolName( pSym->Name, name, sizeof(name)-1, UNDNAME_COMPLETE );
                    if ( rtn == 0 )
                        info[i].function = std::string(pSym->Name);
                    else
                        info[i].function = std::string(name);
                } else {
                    printf( "ERROR: SymGetSymFromAddr (%d,%p)\n", GetLastError(), address2 );
                }

                // Get line number
                IMAGEHLP_LINE64 Line;
                memset( &Line, 0, sizeof( Line ) );
                Line.SizeOfStruct = sizeof( Line );
                DWORD offsetFromLine;
                if ( SymGetLineFromAddr64( pid, address2, &offsetFromLine, &Line ) != FALSE ) {
                    info[i].line     = Line.LineNumber;
                    info[i].filename = std::string( Line.FileName );
                } else {
                    info[i].line     = 0;
                    info[i].filename = std::string();
                }

                // Get the object
                if ( SymGetModuleInfo64( pid, address2, &Module ) != FALSE ) {
                    //info[i].object = std::string( Module.ModuleName );
                    info[i].object = std::string( Module.LoadedImageName );
                    //info[i].baseOfImage = Module.BaseOfImage;
                }
            }
        #else
            for (size_t i=0; i<address.size(); i++) {
                info[i].address = address[i];
                #if defined(_GNU_SOURCE) || defined(USE_MAC)
                    Dl_info dlinfo;
                    if ( !dladdr( info[i].address, &dlinfo ) ) {
                        getDataFromGlobalSymbols( info[i] );
                        continue;
                    }
                    info[i].address2 = subtractAddress( info[i].address, dlinfo.dli_fbase );
                    info[i].object   = std::string( dlinfo.dli_fname );
                    #if defined( USE_ABI )
                        int status;
                        char *demangled = abi::__cxa_demangle( dlinfo.dli_sname, nullptr, nullptr, &status );
                        if ( status == 0 && demangled != nullptr ) {
                            info[i].function = std::string( demangled );
                        } else if ( dlinfo.dli_sname != nullptr ) {
                            info[i].function = std::string( dlinfo.dli_sname );
                        }
                        free( demangled );
                    #else
                        if ( dlinfo.dli_sname != NULL )
                            info[i].function = std::string( dlinfo.dli_sname );
                    #endif
                #else
                    getDataFromGlobalSymbols( info[i] );
                #endif
            }
            // Get the filename / line numbers for each item on the stack
            getFileAndLine( info );
        #endif
    } catch ( ... ) {
    }
    signal( SIGINT, prev_handler ) ;
    return info;
}


/****************************************************************************
*  Function to get the backtrace                                            *
****************************************************************************/
#if defined( USE_LINUX ) || defined( USE_MAC )
static std::vector<void*> thread_backtrace;
static bool thread_backtrace_finished;
static std::mutex thread_backtrace_mutex;
static void _callstack_signal_handler( int, siginfo_t*, void* )
{
    thread_backtrace = StackTrace::backtrace( );
    thread_backtrace_finished = true;
}
#endif
std::vector<void*> StackTrace::backtrace( std::thread::native_handle_type tid )
{
    std::vector<void*> trace;
    #if defined( USE_LINUX ) || defined( USE_MAC )
        // Get the trace
        if ( tid == pthread_self() ) {
            trace.resize(1000,nullptr);
            int trace_size = ::backtrace( trace.data(), trace.size() );
            trace.resize (trace_size );
        } else {
            // Note: this will get the backtrace, but terminates the thread in the process!!!
            thread_backtrace_mutex.lock();
            struct sigaction sa;
            sigfillset(&sa.sa_mask);
            sa.sa_flags = SA_SIGINFO;
            sa.sa_sigaction = _callstack_signal_handler;
            sigaction(CALLSTACK_SIG, &sa, nullptr);
            thread_backtrace_finished = false;
            pthread_kill( tid, CALLSTACK_SIG );
            auto t1 = std::chrono::high_resolution_clock::now();
            auto t2 = std::chrono::high_resolution_clock::now();
            while ( !thread_backtrace_finished && std::chrono::duration<double>(t2-t1).count()<0.1 ) {
                std::this_thread::yield();
                t2 = std::chrono::high_resolution_clock::now();
            }
            std::swap( trace, thread_backtrace );
            thread_backtrace_finished = false;
            thread_backtrace_mutex.unlock();
        }
    #elif defined( USE_WINDOWS )
        #if defined(DBGHELP)

            // Load the modules for the stack trace
            LoadModules();

            // Initialize stackframe for first call
            ::CONTEXT context;
            memset( &context, 0, sizeof( context ) );
            context.ContextFlags = CONTEXT_FULL;
            RtlCaptureContext( &context );
            STACKFRAME64 frame; // in/out stackframe
            memset( &frame, 0, sizeof( frame ) );
            #ifdef _M_IX86
                DWORD imageType = IMAGE_FILE_MACHINE_I386;
                frame.AddrPC.Offset    = context.Eip;
                frame.AddrPC.Mode      = AddrModeFlat;
                frame.AddrFrame.Offset = context.Ebp;
                frame.AddrFrame.Mode   = AddrModeFlat;
                frame.AddrStack.Offset = context.Esp;
                frame.AddrStack.Mode   = AddrModeFlat;
            #elif _M_X64
                DWORD imageType = IMAGE_FILE_MACHINE_AMD64;
                frame.AddrPC.Offset    = context.Rip;
                frame.AddrPC.Mode      = AddrModeFlat;
                frame.AddrFrame.Offset = context.Rsp;
                frame.AddrFrame.Mode   = AddrModeFlat;
                frame.AddrStack.Offset = context.Rsp;
                frame.AddrStack.Mode   = AddrModeFlat;
            #elif _M_IA64
                DWORD imageType = IMAGE_FILE_MACHINE_IA64;
                frame.AddrPC.Offset     = context.StIIP;
                frame.AddrPC.Mode       = AddrModeFlat;
                frame.AddrFrame.Offset  = context.IntSp;
                frame.AddrFrame.Mode    = AddrModeFlat;
                frame.AddrBStore.Offset = context.RsBSP;
                frame.AddrBStore.Mode   = AddrModeFlat;
                frame.AddrStack.Offset  = context.IntSp;
                frame.AddrStack.Mode    = AddrModeFlat;
            #else
                #error "Platform not supported!"
            #endif

            trace.reserve( 1000 );
            auto pid = GetCurrentProcess();
            for ( int frameNum = 0; frameNum<1024; ++frameNum ) {
                BOOL rtn = StackWalk64( imageType, pid, tid, &frame, &context, readProcMem,
                                        SymFunctionTableAccess, SymGetModuleBase64, NULL );
                if ( !rtn ) {
                    printf( "ERROR: StackWalk64 (%p)\n", frame.AddrPC.Offset );
                    break;
                }

                if ( frame.AddrPC.Offset != 0 )
                    trace.push_back( reinterpret_cast<void*>( frame.AddrPC.Offset ) );

                if ( frame.AddrReturn.Offset == 0 )
                    break;
            }
            SetLastError( ERROR_SUCCESS );
        #endif
    #else
        #warning Stack trace is not supported on this compiler/OS
    #endif
    return trace;
}
std::vector<void*> StackTrace::backtrace()
{
    std::vector<void*> trace = backtrace( thisThread() );
    return trace;
}


/****************************************************************************
*  Function to get the list of all active threads                           *
****************************************************************************/
#if defined( USE_LINUX )
static std::thread::native_handle_type thread_handle;
static void _activeThreads_signal_handler( int )
{
    auto handle = StackTrace::thisThread( );
    thread_handle = handle;
    thread_backtrace_finished = true;
}
static inline int get_tid( int pid, const std::string& line )
{
    char buf2[128];
    int i1 = 0;
    while ( line[i1]==' ' && line[i1]!=0 ) { i1++; }
    int i2 = i1;
    while ( line[i2]!=' ' && line[i2]!=0 ) { i2++; }
    memcpy(buf2,&line[i1],i2-i1);
    buf2[i2-i1+1] = 0;
    int pid2 = atoi(buf2);
    if ( pid2 != pid )
        return -1;
    i1 = i2;
    while ( line[i1]==' ' && line[i1]!=0 ) { i1++; }
    i2 = i1;
    while ( line[i2]!=' ' && line[i2]!=0 ) { i2++; }
    memcpy(buf2,&line[i1],i2-i1);
    buf2[i2-i1+1] = 0;
    int tid = atoi(buf2);
    return tid;
}
#endif
std::thread::native_handle_type StackTrace::thisThread( )
{
    #if defined( USE_LINUX ) || defined( USE_MAC )
        return pthread_self();
    #elif defined( USE_WINDOWS )
        return GetCurrentThread();
    #else
        #warning Stack trace is not supported on this compiler/OS
        return std::thread::native_handle_type();
    #endif
}
std::vector<std::thread::native_handle_type> StackTrace::activeThreads( )
{
    std::set<std::thread::native_handle_type> threads;
    #if defined( USE_LINUX )
        std::set<int> tid;
        int pid = getpid();
        char cmd[128];
        sprintf( cmd, "ps -T -p %i", pid );
        signal( SIGCHLD, SIG_DFL );     // Clear child exited
        auto output = breakString( exec( cmd ) );
        for ( const auto& line : output ) {
            int tid2 = get_tid( pid, line );
            if ( tid2 != -1 )
                tid.insert( tid2 );
        }
        tid.erase( syscall(SYS_gettid) );
        signal( CALLSTACK_SIG, _activeThreads_signal_handler );
        for ( auto tid2 : tid ) {
            thread_backtrace_mutex.lock();
            thread_backtrace_finished = false;
            thread_handle = thisThread();
            syscall( SYS_tgkill, pid, tid2, CALLSTACK_SIG );
            auto t1 = std::chrono::high_resolution_clock::now();
            auto t2 = std::chrono::high_resolution_clock::now();
            while ( !thread_backtrace_finished && std::chrono::duration<double>(t2-t1).count()<0.1 ) {
                std::this_thread::yield();
                t2 = std::chrono::high_resolution_clock::now();
            }
            threads.insert( thread_handle );
            thread_backtrace_mutex.unlock();
        }
    #elif defined( USE_MAC )
        printf("activeThreads not finished\n");
    #elif defined( USE_WINDOWS )
        HANDLE hThreadSnap = CreateToolhelp32Snapshot( TH32CS_SNAPTHREAD, 0 ); 
        if( hThreadSnap != INVALID_HANDLE_VALUE ) {
            // Fill in the size of the structure before using it
            THREADENTRY32 te32
            te32.dwSize = sizeof(THREADENTRY32 );
            // Retrieve information about the first thread, and exit if unsuccessful
            if( !Thread32First( hThreadSnap, &te32 ) ) {
                printError( TEXT("Thread32First") );    // Show cause of failure
                CloseHandle( hThreadSnap );             // Must clean up the snapshot object!
                return( FALSE );
            }
            // Now walk the thread list of the system
            do { 
                if ( te32.th32OwnerProcessID == dwOwnerPID )
                    threads.insert( te32.th32ThreadID );
            } while( Thread32Next(hThreadSnap, &te32 ) );
            CloseHandle( hThreadSnap );                 // Must clean up the snapshot object!
        }
    #else
        #warning activeThreads is not yet supported on this compiler/OS
    #endif
    std::vector<std::thread::native_handle_type> threads2;
    threads2.push_back( thisThread() );
    for ( const auto& tmp : threads )
        threads2.push_back( tmp );
    return threads2;
}
// clang-format on


/****************************************************************************
 *  Function to get the current call stack                                   *
 ****************************************************************************/
std::vector<StackTrace::stack_info> StackTrace::getCallStack()
{
    auto trace = StackTrace::backtrace();
    auto info  = getStackInfo( trace );
    return info;
}
std::vector<StackTrace::stack_info> StackTrace::getCallStack( std::thread::native_handle_type id )
{
    auto trace = StackTrace::backtrace( id );
    auto info  = getStackInfo( trace );
    return info;
}
static void fillMultiStackInfoHelper( std::vector<StackTrace::multi_stack_info> &stack,
                                      std::map<void *, StackTrace::stack_info> &stack_data )
{
    for ( auto &i : stack ) {
        i.stack = stack_data[i.stack.address];
        fillMultiStackInfoHelper( i.children, stack_data );
    }
}
static void getAddresses( const std::vector<StackTrace::multi_stack_info> &stack,
                          std::set<void *> &addresses )
{
    for ( const auto &i : stack ) {
        addresses.insert( i.stack.address );
        getAddresses( i.children, addresses );
    }
}
static void fillMultiStackInfo( std::vector<StackTrace::multi_stack_info> &stack )
{
    // Get the list of addresses we encountered
    std::set<void *> addresses;
    getAddresses( stack, addresses );
    // Get the stack info
    std::vector<void *> addresses2( addresses.begin(), addresses.end() );
    auto stack_data = StackTrace::getStackInfo( addresses2 );
    std::map<void *, StackTrace::stack_info> map_data;
    for ( size_t i = 0; i < addresses2.size(); i++ )
        map_data.insert( std::make_pair( addresses2[i], stack_data[i] ) );
    // Fill the data
    fillMultiStackInfoHelper( stack, map_data );
}
static std::vector<StackTrace::multi_stack_info>
generateMultiStack( const std::vector<std::vector<void *>> &thread_backtrace )
{
    std::vector<StackTrace::multi_stack_info> stack;
    for ( const auto &trace : thread_backtrace ) {
        if ( trace.empty() )
            continue;
        std::vector<StackTrace::multi_stack_info> *parent = &stack;
        for ( int i = trace.size() - 1; i >= 0; i-- ) {
            void *ptr                           = trace[i];
            StackTrace::multi_stack_info *child = nullptr;
            for ( auto &tmp : *parent ) {
                if ( tmp.stack.address == ptr ) {
                    child = &tmp;
                    break;
                }
            }
            if ( child == nullptr ) {
                parent->resize( parent->size() + 1 );
                child                = &( parent->back() );
                child->N             = 0;
                child->stack.address = ptr;
            }
            ( child->N )++;
            parent = &( child->children );
        }
    }
    return stack;
}
std::vector<StackTrace::multi_stack_info> StackTrace::getAllCallStacks()
{
    // Get the list of threads
    auto threads = activeThreads();
    // Get the backtrace of each thread
    std::vector<std::vector<void *>> thread_backtrace;
    for ( auto thread : threads )
        thread_backtrace.push_back( backtrace( thread ) );
    // Create the multi-stack strucutre
    auto stack = generateMultiStack( thread_backtrace );
    // Fill the data for the stacks
    fillMultiStackInfo( stack );
    return stack;
}


/****************************************************************************
 *  Function to get system search paths                                      *
 ****************************************************************************/
std::string StackTrace::getSymPaths()
{
    std::string paths;
#ifdef USE_WINDOWS
    // Create the path list (seperated by ';' )
    paths = std::string( ".;" );
    paths.reserve( 1000 );
    // Add the current directory
    paths += getCurrentDirectory() + ";";
    // Now add the path for the main-module:
    char temp[1024];
    memset( temp, 0, sizeof( temp ) );
    if ( GetModuleFileNameA( nullptr, temp, sizeof( temp ) - 1 ) > 0 ) {
        for ( char *p = ( temp + strlen( temp ) - 1 ); p >= temp; --p ) {
            // locate the rightmost path separator
            if ( ( *p == '\\' ) || ( *p == '/' ) || ( *p == ':' ) ) {
                *p = 0;
                break;
            }
        }
        if ( strlen( temp ) > 0 ) {
            paths += temp;
            paths += ";";
        }
    }
    memset( temp, 0, sizeof( temp ) );
    if ( GetEnvironmentVariableA( "_NT_SYMBOL_PATH", temp, sizeof( temp ) - 1 ) > 0 ) {
        paths += temp;
        paths += ";";
    }
    memset( temp, 0, sizeof( temp ) );
    if ( GetEnvironmentVariableA( "_NT_ALTERNATE_SYMBOL_PATH", temp, sizeof( temp ) - 1 ) > 0 ) {
        paths += temp;
        paths += ";";
    }
    memset( temp, 0, sizeof( temp ) );
    if ( GetEnvironmentVariableA( "SYSTEMROOT", temp, sizeof( temp ) - 1 ) > 0 ) {
        paths += temp;
        paths += ";";
        // also add the "system32"-directory:
        paths += temp;
        paths += "\\system32;";
    }
    memset( temp, 0, sizeof( temp ) );
    if ( GetEnvironmentVariableA( "SYSTEMDRIVE", temp, sizeof( temp ) - 1 ) > 0 ) {
        paths += "SRV*;" + std::string( temp ) +
                 "\\websymbols*http://msdl.microsoft.com/download/symbols;";
    } else {
        paths += "SRV*c:\\websymbols*http://msdl.microsoft.com/download/symbols;";
    }
#endif
    return paths;
}


/****************************************************************************
 *  Load modules for windows                                                 *
 ****************************************************************************/
#ifdef USE_WINDOWS
BOOL StackTrace::GetModuleListTH32( HANDLE hProcess, DWORD pid )
{
    // CreateToolhelp32Snapshot()
    typedef HANDLE( __stdcall * tCT32S )( DWORD dwFlags, DWORD th32ProcessID );
    // Module32First()
    typedef BOOL( __stdcall * tM32F )( HANDLE hSnapshot, LPMODULEENTRY32 lpme );
    // Module32Next()
    typedef BOOL( __stdcall * tM32N )( HANDLE hSnapshot, LPMODULEENTRY32 lpme );

    // try both dlls...
    const TCHAR *dllname[] = { _T("kernel32.dll"), _T("tlhelp32.dll") };
    HINSTANCE hToolhelp    = nullptr;
    tCT32S pCT32S          = nullptr;
    tM32F pM32F            = nullptr;
    tM32N pM32N            = nullptr;

    HANDLE hSnap;
    MODULEENTRY32 me;
    me.dwSize = sizeof( me );

    for ( size_t i = 0; i < ( sizeof( dllname ) / sizeof( dllname[0] ) ); i++ ) {
        hToolhelp = LoadLibrary( dllname[i] );
        if ( hToolhelp == nullptr )
            continue;
        pCT32S = (tCT32S) GetProcAddress( hToolhelp, "CreateToolhelp32Snapshot" );
        pM32F  = (tM32F) GetProcAddress( hToolhelp, "Module32First" );
        pM32N  = (tM32N) GetProcAddress( hToolhelp, "Module32Next" );
        if ( ( pCT32S != nullptr ) && ( pM32F != nullptr ) && ( pM32N != nullptr ) )
            break; // found the functions!
        FreeLibrary( hToolhelp );
        hToolhelp = nullptr;
    }

    if ( hToolhelp == nullptr )
        return FALSE;

    hSnap = pCT32S( TH32CS_SNAPMODULE, pid );
    if ( hSnap == (HANDLE) -1 ) {
        FreeLibrary( hToolhelp );
        return FALSE;
    }

    bool keepGoing = !!pM32F( hSnap, &me );
    int cnt        = 0;
    while ( keepGoing ) {
        LoadModule( hProcess, me.szExePath, me.szModule, (DWORD64) me.modBaseAddr, me.modBaseSize );
        cnt++;
        keepGoing = !!pM32N( hSnap, &me );
    }
    CloseHandle( hSnap );
    FreeLibrary( hToolhelp );
    if ( cnt <= 0 )
        return FALSE;
    return TRUE;
}
DWORD StackTrace::LoadModule(
    HANDLE hProcess, LPCSTR img, LPCSTR mod, DWORD64 baseAddr, DWORD size )
{
    CHAR *szImg  = _strdup( img );
    CHAR *szMod  = _strdup( mod );
    DWORD result = ERROR_SUCCESS;
    if ( ( szImg == nullptr ) || ( szMod == nullptr ) ) {
        result = ERROR_NOT_ENOUGH_MEMORY;
    } else {
        if ( SymLoadModule( hProcess, 0, szImg, szMod, baseAddr, size ) == 0 )
            result = GetLastError();
    }
    ULONGLONG fileVersion = 0;
    if ( szImg != nullptr ) {
        // try to retrive the file-version:
        VS_FIXEDFILEINFO *fInfo = nullptr;
        DWORD dwHandle;
        DWORD dwSize = GetFileVersionInfoSizeA( szImg, &dwHandle );
        if ( dwSize > 0 ) {
            LPVOID vData = malloc( dwSize );
            if ( vData != nullptr ) {
                if ( GetFileVersionInfoA( szImg, dwHandle, dwSize, vData ) != 0 ) {
                    UINT len;
                    TCHAR szSubBlock[] = _T("\\");
                    if ( VerQueryValue( vData, szSubBlock, (LPVOID *) &fInfo, &len ) == 0 ) {
                        fInfo = nullptr;
                    } else {
                        fileVersion = ( (ULONGLONG) fInfo->dwFileVersionLS ) +
                                      ( (ULONGLONG) fInfo->dwFileVersionMS << 32 );
                    }
                }
                free( vData );
            }
        }

        // Retrive some additional-infos about the module
        IMAGEHLP_MODULE64 Module;
        Module.SizeOfStruct = sizeof( IMAGEHLP_MODULE64 );
        SymGetModuleInfo64( hProcess, baseAddr, &Module );
        LPCSTR pdbName = Module.LoadedImageName;
        if ( Module.LoadedPdbName[0] != 0 )
            pdbName = Module.LoadedPdbName;
    }
    if ( szImg != nullptr )
        free( szImg );
    if ( szMod != nullptr )
        free( szMod );
    return result;
}
BOOL StackTrace::GetModuleListPSAPI( HANDLE hProcess )
{
    DWORD cbNeeded;
    HMODULE hMods[1024];
    char tt[8192];
    char tt2[8192];
    if ( !EnumProcessModules( hProcess, hMods, sizeof( hMods ), &cbNeeded ) ) {
        return false;
    }
    if ( cbNeeded > sizeof( hMods ) ) {
        printf( "Insufficient memory allocated in GetModuleListPSAPI\n" );
        return false;
    }
    int cnt = 0;
    for ( DWORD i = 0; i < cbNeeded / sizeof( hMods[0] ); i++ ) {
        // base address, size
        MODULEINFO mi;
        GetModuleInformation( hProcess, hMods[i], &mi, sizeof( mi ) );
        // image file name
        tt[0] = 0;
        GetModuleFileNameExA( hProcess, hMods[i], tt, sizeof( tt ) );
        // module name
        tt2[0] = 0;
        GetModuleBaseNameA( hProcess, hMods[i], tt2, sizeof( tt2 ) );
        DWORD dwRes = LoadModule( hProcess, tt, tt2, (DWORD64) mi.lpBaseOfDll, mi.SizeOfImage );
        if ( dwRes != ERROR_SUCCESS )
            printf( "ERROR: LoadModule (%d)\n", dwRes );
        cnt++;
    }

    return cnt != 0;
}
void StackTrace::LoadModules()
{
    static bool modules_loaded = false;
    if ( !modules_loaded ) {
        modules_loaded = true;

        // Get the search paths for symbols
        std::string paths = StackTrace::getSymPaths();

        // Initialize the symbols
        if ( SymInitialize( GetCurrentProcess(), paths.c_str(), FALSE ) == FALSE )
            printf( "ERROR: SymInitialize (%d)\n", GetLastError() );

        DWORD symOptions = SymGetOptions();
        symOptions |= SYMOPT_LOAD_LINES | SYMOPT_FAIL_CRITICAL_ERRORS;
        symOptions     = SymSetOptions( symOptions );
        char buf[1024] = { 0 };
        if ( SymGetSearchPath( GetCurrentProcess(), buf, sizeof( buf ) ) == FALSE )
            printf( "ERROR: SymGetSearchPath (%d)\n", GetLastError() );

        // First try to load modules from toolhelp32
        BOOL loaded = StackTrace::GetModuleListTH32( GetCurrentProcess(), GetCurrentProcessId() );

        // Try to load from Psapi
        if ( !loaded )
            loaded = StackTrace::GetModuleListPSAPI( GetCurrentProcess() );
    }
}
#endif


/****************************************************************************
 *  Get the signal name                                                      *
 ****************************************************************************/
std::string StackTrace::signalName( int sig ) { return std::string( strsignal( sig ) ); }
std::vector<int> StackTrace::allSignalsToCatch()
{
    std::set<int> signals;
    for ( int i = 1; i < 32; i++ )
        signals.insert( i );
    for ( int i = SIGRTMIN; i <= SIGRTMAX; i++ )
        signals.insert( i );
    signals.erase( SIGKILL );
    signals.erase( SIGSTOP );
    return std::vector<int>( signals.begin(), signals.end() );
}
std::vector<int> StackTrace::defaultSignalsToCatch()
{
    auto tmp = allSignalsToCatch();
    std::set<int> signals( tmp.begin(), tmp.end() );
    signals.erase( SIGWINCH ); // Don't catch window changed by default
    signals.erase( SIGCONT );  // Don't catch continue by default
    return std::vector<int>( signals.begin(), signals.end() );
}


/****************************************************************************
 *  Set the signal handlers                                                  *
 ****************************************************************************/
static std::function<void( std::string, StackTrace::terminateType )> abort_fun;
static std::string rethrow()
{
    std::string last_message;
#ifdef USE_LINUX
    try {
        static int tried_throw = 0;
        if ( tried_throw == 0 ) {
            tried_throw = 1;
            throw;
        }
        // No active exception
    } catch ( const std::exception &err ) {
        // Caught a std::runtime_error
        last_message = err.what();
    } catch ( ... ) {
        // Caught an unknown exception
        last_message = "unknown exception occurred.";
    }
#endif
    return last_message;
}
static void term_func_abort( int sig )
{
    std::string msg( "Caught signal: " );
    msg += StackTrace::signalName( sig );
    abort_fun( msg, StackTrace::terminateType::signal );
}
static std::set<int> signals_set = std::set<int>();
static void term_func()
{
    std::string last_message = rethrow();
    StackTrace::clearSignals();
    abort_fun( "Unhandled exception:\n" + last_message, StackTrace::terminateType::signal );
}
void StackTrace::clearSignal( int sig )
{
    if ( signals_set.find( sig ) != signals_set.end() ) {
        signal( sig, SIG_DFL );
        signals_set.erase( sig );
    }
}
void StackTrace::clearSignals()
{
    for ( auto sig : signals_set )
        signal( sig, SIG_DFL );
    signals_set.clear();
}
void StackTrace::setSignals( const std::vector<int> &signals, void ( *handler )( int ) )
{
    for ( auto sig : signals ) {
        signal( sig, handler );
        signals_set.insert( sig );
    }
}
void StackTrace::setErrorHandlers(
    std::function<void( std::string, StackTrace::terminateType )> abort )
{
    abort_fun = abort;
    std::set_terminate( term_func );
    setSignals( defaultSignalsToCatch(), &term_func_abort );
    std::set_unexpected( term_func );
}


} // namespace AMP
