#include "Utilities.h"

#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/Logger.h"
#include "utils/PIO.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <signal.h>

// Detect the OS and include system dependent headers
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64) || defined(_MSC_VER)
    // Note: windows has not been testeds
    #define USE_WINDOWS
    #include <windows.h>
    #include <stdio.h>   
    #include <tchar.h>
    #include <psapi.h>
    #include <DbgHelp.h>
    #define mkdir(path, mode) _mkdir(path)
    //#pragma comment(lib, psapi.lib) //added
    //#pragma comment(linker, /DEFAULTLIB:psapi.lib)
#elif defined(__APPLE__)
    #define USE_MAC
    #include <signal.h>
    #include <execinfo.h>
    #include <cxxabi.h>
    #include <dlfcn.h>
    #include <mach/mach.h>
    #include <stdint.h>
    #include <sys/types.h>
    #include <sys/sysctl.h>
#elif defined(__linux) || defined(__unix) || defined(__posix)
    #define USE_LINUX
    #define USE_NM
    #include <sys/time.h>
    #include <execinfo.h>
    #include <dlfcn.h>
    #include <malloc.h>
#else
    #error Unknown OS
#endif


#ifdef __GNUC__
    #define USE_ABI
    #include <cxxabi.h>
#endif


namespace AMP{

// Routine to check if a file exists
bool Utilities::fileExists( const std::string& filename )
{
   std::ifstream ifile(filename.c_str());
   return ifile!=NULL;
}

// Routine to rename a file.
void Utilities::renameFile(const std::string& old_filename, const std::string& new_filename)
{
   AMP_ASSERT(!old_filename.empty());
   AMP_ASSERT(!new_filename.empty());
   rename(old_filename.c_str(), new_filename.c_str());
}

// Routine to delete a file.
void Utilities::deleteFile(const std::string& filename)
{
   AMP_ASSERT(!filename.empty());
   if ( fileExists( filename ) ) {
      int error = remove(filename.c_str());
      AMP_INSIST(error==0,"Error deleting file");
   }
}

/*
 * Routine to recursively construct directories based on a relative path name.
 */
void Utilities::recursiveMkdir(
   const std::string& path, 
   mode_t mode,
   bool only_node_zero_creates)
{
   AMP_MPI comm = AMP_MPI(AMP_COMM_WORLD);
   if ( (!only_node_zero_creates) || (comm.getRank() == 0)) {
      int length = (int) path.length();
      char *path_buf= new char[length+1];
      sprintf(path_buf,"%s",path.c_str());
      struct stat status;
      int pos = length - 1;
   
      /* find part of path that has not yet been created */
      while ( (stat(path_buf,&status) != 0) && (pos >= 0) ) {
   
         /* slide backwards in string until next slash found */
         bool slash_found = false;
         while ( (!slash_found) && (pos >= 0) ) {
           if ( path_buf[pos]=='/' || path_buf[pos]==92 ) {
              slash_found = true;
              if (pos >= 0) path_buf[pos] = '\0';
           } else pos--;
         }
      } 

      /* 
       * if there is a part of the path that already exists make sure
       * it is really a directory
       */
      if (pos >= 0) {
         if ( !S_ISDIR(status.st_mode) ) {
            AMP_ERROR("Error in Utilities::recursiveMkdir...\n"
               << "    Cannot create directories in path = " << path
               << "\n    because some intermediate item in path exists and"
               << "is NOT a directory" << std::endl);
         }
      }
   
      /* make all directories that do not already exist */
   
      /* 
       * if (pos < 0), then there is no part of the path that
       * already exists.  Need to make the first part of the 
       * path before sliding along path_buf.
       */
      if (pos < 0) {
	 if(mkdir(path_buf,mode) != 0) {
	    AMP_ERROR("Error in Utilities::recursiveMkdir...\n"
		       << "    Cannot create directory  = " 
		       << path_buf << std::endl);
	 }
	 pos = 0;
      }
   
      /* make rest of directories */
      do {
   
         /* slide forward in string until next '\0' found */
         bool null_found = false;
         while ( (!null_found) && (pos < length) ) {
           if (path_buf[pos] == '\0') {
              null_found = true;
              path_buf[pos] = '/';
           }
           pos++;
         }
   
         /* make directory if not at end of path */
	 if (pos < length) {
	    if(mkdir(path_buf,mode) != 0) {
	       AMP_ERROR("Error in Utilities::recursiveMkdir...\n"
			  << "    Cannot create directory  = " 
			  << path_buf << std::endl);
	    }
	 }
      } while (pos < length);

      delete [] path_buf;
   }

   /* 
    * Make sure all processors wait until node zero creates 
    * the directory structure.
    */
   if (only_node_zero_creates) {
     comm.barrier();
   }
}

/*
 * Routine to convert an integer to a string. 
 */
std::string Utilities::intToString(int num, int min_width)
{
   int tmp_width = ( min_width > 0 ? min_width : 1 );
   std::ostringstream os;
   if ( num < 0 ) {
      os << '-' << std::setw(tmp_width-1) << std::setfill('0') << -num;
   } else {
      os << std::setw(tmp_width) << std::setfill('0') << num;
   }
   os << std::flush;

   return(os.str()); //returns the string form of the stringstream object
}

std::string Utilities::nodeToString(int num) 
{
   return intToString(num, 5);
}

std::string Utilities::processorToString(int num) 
{
   return intToString(num, 5);
}


std::string Utilities::patchToString(int num) {
   return intToString(num, 4);
}


std::string Utilities::levelToString(int num) {
   return intToString(num, 4);
}


std::string Utilities::blockToString(int num) {
   return intToString(num, 4);
}


/*
 * Routine that calls abort and prints calling location to error stream.
 */
void Utilities::abort(const std::string &message, 
	              const std::string &filename, 
	              const int line) 
{
    if ( AMP::AMPManager::use_MPI_Abort==true) {
        // Print the call stack and memory usage
        long long unsigned int N_bytes = getMemoryUsage();
        printf("Bytes used = %llu\n",N_bytes);
        std::vector<std::string> stack = getCallStack();
        printf("Stack Trace:\n");
        for (size_t i=0; i<stack.size(); i++)
            printf("   %s\n",stack[i].c_str());
        printf("\n");
        // Log the abort message
        Logger::getInstance() -> logAbort(message, filename, line);
        // Use MPI_abort (will terminate all processes)
        AMP_MPI comm = AMP_MPI(AMP_COMM_WORLD);
        comm.abort();
    } else {
        // Throw and standard exception (allows the use of try, catch)
        // std::stringstream  stream;
        // stream << message << std::endl << "  " << filename << ":  " << line;
        // std::cout << stream.str() << std::endl;
        throw std::logic_error(message);
    }
}


// Function to create a 32-bit hash key from a character array
unsigned int Utilities::hash_char(const char* name)
{
    AMP_INSIST(sizeof(unsigned int)==4,"Need unsigned 32-bit int");
    unsigned int hash = 5381;
    unsigned char c;
    while((c = *name++)) {
        // hash = hash * 33 ^ c
        hash = ((hash << 5) + hash) ^ c;
    }
    return hash;
}


/****************************************************************************
*  Function to get the memory usage                                         *
*  Note: this function should be thread-safe                                *
****************************************************************************/
#if defined(USE_MAC) || defined(USE_LINUX)
    // Get the page size on mac or linux
    static size_t page_size = static_cast<size_t>(sysconf(_SC_PAGESIZE));
#endif
static size_t N_bytes_initialization = Utilities::getMemoryUsage();
size_t Utilities::getSystemMemory()
{
    size_t N_bytes = 0;
    #if defined(USE_LINUX)
        static long pages = sysconf(_SC_PHYS_PAGES);
        N_bytes = pages * page_size;
    #elif defined(USE_MAC)
        int mib[2] = { CTL_HW, HW_MEMSIZE };
        u_int namelen = sizeof(mib) / sizeof(mib[0]);
        uint64_t size;
        size_t len = sizeof(size);
        if (sysctl(mib, namelen, &size, &len, NULL, 0) == 0)
            N_bytes = size;
    #elif defined(USE_WINDOWS)
        MEMORYSTATUSEX status;
        status.dwLength = sizeof(status);
        GlobalMemoryStatusEx(&status);
        N_bytes = status.ullTotalPhys;
    #endif
    return N_bytes;
}
size_t Utilities::getMemoryUsage()
{
    size_t N_bytes = 0;
    #if defined(USE_LINUX)
        struct mallinfo meminfo = mallinfo();
        size_t size_hblkhd = static_cast<size_t>( meminfo.hblkhd );
        size_t size_uordblks = static_cast<size_t>( meminfo.uordblks );
        N_bytes = static_cast<size_t>( size_hblkhd + size_uordblks );
    #elif defined(USE_MAC)
        struct task_basic_info t_info;
        mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
        if (KERN_SUCCESS != task_info(mach_task_self(),
                              TASK_BASIC_INFO, (task_info_t)&t_info, 
                              &t_info_count)) {
            return 0;
        }
        N_bytes = t_info.virtual_size;
    #elif defined(USE_WINDOWS)
        PROCESS_MEMORY_COUNTERS memCounter;
        GetProcessMemoryInfo( GetCurrentProcess(), &memCounter, sizeof(memCounter) );
        N_bytes = memCounter.WorkingSetSize;
    #endif
    return N_bytes;
}


/****************************************************************************
*  Function to get the current call stack                                   *
****************************************************************************/
struct global_symbols_struct {
    bool loaded;
    std::vector<void*> address;
    std::vector<std::string> obj;
    global_symbols_struct(): loaded(false) {}
} global_symbols;
std::string print_address( void *address )
{
    char tmp[20];
    unsigned long long int address2 = reinterpret_cast<size_t>(address);
    sprintf(tmp,"0x%016llx",address2);
    return std::string(tmp);
}
std::string remove_path(const std::string& var)
{
    const char* tmp = std::max(strrchr(var.c_str(),47),strrchr(var.c_str(),92));
    if ( tmp==0 )
        tmp = var.c_str();
    else
        tmp++;
    return std::string(tmp);
}
std::vector<std::string>  Utilities::getCallStack()
{
    std::vector<std::string>  stack_list;
    #if defined(USE_LINUX) || defined(USE_MAC)
        void *trace[100];
        memset(trace,0,100*sizeof(void*));
        Dl_info dlinfo;
        int status;
        const char *symname;
        char *demangled=NULL;
        int trace_size = backtrace(trace,100);
        char **names = backtrace_symbols(trace,100);
        for (int i=0; i<trace_size; ++i) {  
            if(!dladdr(trace[i], &dlinfo))
                continue;
            symname = dlinfo.dli_sname;
            #if defined(USE_ABI)
                demangled = abi::__cxa_demangle(symname, NULL, 0, &status);
                if(status == 0 && demangled)
                    symname = demangled;
            #endif
            std::string object = std::string(dlinfo.dli_fname);
            std::string function = "";
            if ( symname!=NULL )
                function = std::string(symname);
            // Create the stack item
            object = remove_path(object);
            std::string stack_item = print_address(trace[i]) + ":   " + object + ":   ";
            while ( stack_item.size() < 40 )
                stack_item.push_back(' ');
            stack_item += function;
            stack_list.push_back(stack_item);
            if ( demangled!=NULL ) {
                free(demangled);
                demangled=NULL;
            }
        }
        free(names);
        if ( stack_list.empty() ) {
            // Try to load the global symbols and use them
            if ( !global_symbols.loaded ) {
                global_symbols.loaded = true;
                std::vector<char> type;
                Utilities::get_symbols( global_symbols.address, type, global_symbols.obj );
            }
            if ( !global_symbols.address.empty() ) {
                for (int i=0; i<trace_size; ++i) {  
                    // Search through the global symbols to find the address
                    size_t index = findfirst(global_symbols.address,trace[i])-1;
                    std::string stack_item = print_address(trace[i]) + ":   " + global_symbols.obj[index];
                    stack_list.push_back(stack_item);
                }
            }
        }
        if ( stack_list.empty() ) {
            // We don't have access to the symbols, print the addresses
            // We can then use nm or some other method to convert the addresses to the function
            for (int i=0; i<trace_size; ++i) {  
                std::string stack_item = print_address(trace[i]);
                stack_list.push_back(stack_item);
            }
        }
    #elif defined(USE_WINDOWS)
        ::CONTEXT lContext;
        ::ZeroMemory( &lContext, sizeof( ::CONTEXT ) );
        ::RtlCaptureContext( &lContext );
        ::STACKFRAME64 lFrameStack;
        ::ZeroMemory( &lFrameStack, sizeof( ::STACKFRAME64 ) );
        lFrameStack.AddrPC.Offset = lContext.Rip;
        lFrameStack.AddrFrame.Offset = lContext.Rbp;
        lFrameStack.AddrStack.Offset = lContext.Rsp;
        lFrameStack.AddrPC.Mode = lFrameStack.AddrFrame.Mode = lFrameStack.AddrStack.Mode = AddrModeFlat;
        #ifdef _M_IX86
            DWORD MachineType = IMAGE_FILE_MACHINE_I386;
        #endif
        #ifdef _M_X64
            DWORD MachineType = IMAGE_FILE_MACHINE_AMD64;
        #endif
        #ifdef _M_IA64
            DWORD MachineType = IMAGE_FILE_MACHINE_IA64;
        #endif
        while ( 1 ) {
            int rtn = ::StackWalk64( MachineType, ::GetCurrentProcess(), ::GetCurrentThread(), 
                &lFrameStack, MachineType == IMAGE_FILE_MACHINE_I386 ? 0 : &lContext,
                NULL, &::SymFunctionTableAccess64, &::SymGetModuleBase64, NULL );
            if( !rtn )
                break;
            if( lFrameStack.AddrPC.Offset == 0 )
                break;
            ::MEMORY_BASIC_INFORMATION lInfoMemory;
            ::VirtualQuery( ( ::PVOID )lFrameStack.AddrPC.Offset, &lInfoMemory, sizeof( lInfoMemory ) );
            if ( lInfoMemory.Type==MEM_PRIVATE )
                continue;
            ::DWORD64 lBaseAllocation = reinterpret_cast< ::DWORD64 >( lInfoMemory.AllocationBase );
            ::TCHAR lNameModule[ 1024 ];
            ::HMODULE hBaseAllocation = reinterpret_cast< ::HMODULE >( lBaseAllocation );
            ::GetModuleFileName( hBaseAllocation, lNameModule, 1024 );
            PIMAGE_DOS_HEADER lHeaderDOS = reinterpret_cast<PIMAGE_DOS_HEADER>( lBaseAllocation );
            if ( lHeaderDOS==NULL )
                continue;
            PIMAGE_NT_HEADERS lHeaderNT = reinterpret_cast<PIMAGE_NT_HEADERS>( lBaseAllocation + lHeaderDOS->e_lfanew );
            PIMAGE_SECTION_HEADER lHeaderSection = IMAGE_FIRST_SECTION( lHeaderNT );
            ::DWORD64 lRVA = lFrameStack.AddrPC.Offset - lBaseAllocation;
            ::DWORD64 lNumberSection = ::DWORD64();
            ::DWORD64 lOffsetSection = ::DWORD64();
            for( int lCnt = ::DWORD64(); lCnt < lHeaderNT->FileHeader.NumberOfSections; lCnt++, lHeaderSection++ ) {
                ::DWORD64 lSectionBase = lHeaderSection->VirtualAddress;
                ::DWORD64 lSectionEnd = lSectionBase + max( lHeaderSection->SizeOfRawData, lHeaderSection->Misc.VirtualSize );
                if( ( lRVA >= lSectionBase ) && ( lRVA <= lSectionEnd ) ) {
                    lNumberSection = lCnt + 1;
                    lOffsetSection = lRVA - lSectionBase;
                    break;
                }
            }
            std::stringstream stream;
            stream << lNameModule << " : 000" << lNumberSection << " : " << reinterpret_cast<void*>(lOffsetSection);
            stack_list.push_back(stream.str());
        }
    #else
        #warning Stack trace is not supported on this compiler/OS
    #endif
    return stack_list;
}


/****************************************************************************
*  Function to get symbols for the executable from nm (if availible)        *
****************************************************************************/
int Utilities::get_symbols( std::vector<void*>& address, std::vector<char>& type, 
    std::vector<std::string>& obj )
{
    address.clear();
    type.clear();
    obj.clear();
    #ifdef USE_NM
        try { 
            char buf[1024];
            int len = ::readlink("/proc/self/exe",buf,sizeof(buf)-1);
            if ( len==-1 )
                return -2;
            buf[len] = '\0';
            char cmd[1024];
            sprintf(cmd,"nm --demangle --numeric-sort %s",buf);
            FILE *in = popen(cmd,"r");
            if ( in==NULL )
                return -2;
            while ( fgets(buf,sizeof(buf)-1,in)!=NULL ) {
                if ( buf[0] == ' ' )
                    continue;
                char *a = buf;
                char *b = strchr(a,' ');  b[0] = 0;  b++;
                char *c = strchr(b,' ');  c[0] = 0;  c++;
                char *d = strchr(c,'\n');  if ( d ) { d[0]=0; }
                size_t add = strtoul(a,NULL,16);
                address.push_back( reinterpret_cast<void*>(add) );
                type.push_back( b[0] );
                obj.push_back( std::string(c) );
        	}
            pclose(in);
        } catch (...) {
            return -3;
        }
        return 0;
    #else
        return -1;
    #endif
}


// Print AMP Banner
void Utilities::printBanner()
{
    std::ostringstream banner;
    banner << std::endl;
    banner << "            _____                    _____                    _____"           << std::endl;
    banner << "           /\\    \\                  /\\    \\                  /\\    \\ "         << std::endl;
    banner << "          /::\\    \\                /::\\____\\                /::\\    \\"         << std::endl;        
    banner << "         /::::\\    \\              /::::|   |               /::::\\    \\"        << std::endl;       
    banner << "        /::::::\\    \\            /:::::|   |              /::::::\\    \\"       << std::endl;      
    banner << "       /:::/\\:::\\    \\          /::::::|   |             /:::/\\:::\\    \\"      << std::endl;     
    banner << "      /:::/__\\:::\\    \\        /:::/|::|   |            /:::/__\\:::\\    \\"     << std::endl;    
    banner << "     /::::\\   \\:::\\    \\      /:::/ |::|   |           /::::\\   \\:::\\    \\"    << std::endl;   
    banner << "    /::::::\\   \\:::\\    \\    /:::/  |::|___|______    /::::::\\   \\:::\\    \\"   << std::endl;  
    banner << "   /:::/\\:::\\   \\:::\\    \\  /:::/   |::::::::\\    \\  /:::/\\:::\\   \\:::\\____\\"  << std::endl; 
    banner << "  /:::/  \\:::\\   \\:::\\____\\/:::/    |:::::::::\\____\\/:::/  \\:::\\   \\:::|    |" << std::endl;
    banner << "  \\::/    \\:::\\  /:::/    /\\::/    / ~~~~~/:::/    /\\::/    \\:::\\  /:::|____|" << std::endl;
    banner << "   \\/____/ \\:::\\/:::/    /  \\/____/      /:::/    /  \\/_____/\\:::\\/:::/    /"  << std::endl; 
    banner << "            \\::::::/    /               /:::/    /            \\::::::/    /"   << std::endl;  
    banner << "             \\::::/    /               /:::/    /              \\::::/    /"    << std::endl;   
    banner << "             /:::/    /               /:::/    /                \\::/____/"     << std::endl;    
    banner << "            /:::/    /               /:::/    /"                               << std::endl;                            
    banner << "           /:::/    /               /:::/    /"                                << std::endl;                               
    banner << "          /:::/    /               /:::/    /"                                 << std::endl;                                
    banner << "          \\::/    /                \\::/    /"                                  << std::endl;                                 
    banner << "           \\/____/                  \\/____/"                                   << std::endl;                                  
    banner << std::endl << std::endl;

    AMP::pout << banner.str();

}

// Factor a number into it's prime factors
std::vector<int> Utilities::factor(size_t number)
{
    if ( number<=3 ) 
        return std::vector<int>(1,(int)number);
    size_t i, n, n_max;
    bool factor_found;
    // Compute the maximum number of factors
    int N_primes_max = 1;
    n = number;
    while (n >>= 1) ++N_primes_max;
    // Initialize n, factors 
    n = number;
    std::vector<int> factors;
    factors.reserve(N_primes_max);
    while ( 1 ) {
        // Check if n is a trivial prime number
        if ( n==2 || n==3 || n==5 ) {
            factors.push_back( (int) n );
            break;
        } 
        // Check if n is divisible by 2
        if ( n%2 == 0 ) {
            factors.push_back( 2 );
            n/=2;
            continue;
        } 
        // Check each odd number until a factor is reached
        n_max = (size_t) floor(sqrt((double) n));
        factor_found = false;
        for (i=3; i<=n_max; i+=2) {
            if ( n%i == 0 ) {
                factors.push_back( i );
                n/=i;
                factor_found = true;
                break;
            } 
        }
        if ( factor_found )
            continue;
        // No factors were found, the number must be prime
        factors.push_back( (int) n );
        break;
    }
    // Sort the factors
    AMP::Utilities::quicksort(factors);
    return factors;
}


// Function to perform linear interpolation
double Utilities::linear( const std::vector<double>& x, 
    const std::vector<double>& f, double xi )
{
    size_t Nx = x.size();
    AMP_ASSERT( Nx>1 );
    AMP_ASSERT( f.size()==Nx );
    size_t i = AMP::Utilities::findfirst( x, xi );
    if ( i==0 ) { i=1; }
    if ( i==x.size() ) { i=x.size()-1; }
    double dx = (xi-x[i-1])/(x[i]-x[i-1]);
    return dx*f[i] + (1.0-dx)*f[i-1];
}


// Function to perform bi-linear interpolation
double Utilities::bilinear( const std::vector<double>& x, const std::vector<double>& y, 
    const std::vector<double>& f, double xi, double yi )
{
    size_t Nx = x.size();
    size_t Ny = y.size();
    AMP_ASSERT( Nx>1 && Ny>1 );
    AMP_ASSERT( f.size() == Nx*Ny );
    size_t i = AMP::Utilities::findfirst( x, xi );
    size_t j = AMP::Utilities::findfirst( y, yi );
    if ( i==0 ) { i=1; }
    if ( j==0 ) { j=1; }
    if ( i==x.size() ) { i=x.size()-1; }
    if ( j==y.size() ) { j=y.size()-1; }
    double dx = (xi-x[i-1])/(x[i]-x[i-1]);
    double dy = (yi-y[j-1])/(y[j]-y[j-1]);
    double f1 = f[ i-1 + (j-1)*Nx ];
    double f2 = f[ i   + (j-1)*Nx ];
    double f3 = f[ i-1 + j*Nx     ];
    double f4 = f[ i   + j*Nx     ];
    double dx2 = 1.0-dx;
    double dy2 = 1.0-dy;
    return (dx*f2 + dx2*f1)*dy2 + (dx*f4 + dx2*f3)*dy;
}


// Function to perform tri-linear interpolation
double Utilities::trilinear( const std::vector<double>& x, const std::vector<double>& y, 
    const std::vector<double>& z, const std::vector<double>& f, double xi, double yi, double zi )
{
    size_t Nx = x.size();
    size_t Ny = y.size();
    size_t Nz = z.size();
    AMP_ASSERT( Nx>1 && Ny>1 && Nz>1 );
    AMP_ASSERT( f.size() == Nx*Ny*Nz );
    size_t i = AMP::Utilities::findfirst( x, xi );
    size_t j = AMP::Utilities::findfirst( y, yi );
    size_t k = AMP::Utilities::findfirst( z, zi );
    if ( i==0 ) { i=1; }
    if ( j==0 ) { j=1; }
    if ( k==0 ) { k=1; }
    if ( i==x.size() ) { i=x.size()-1; }
    if ( j==y.size() ) { j=y.size()-1; }
    if ( k==z.size() ) { k=z.size()-1; }
    double dx = (xi-x[i-1])/(x[i]-x[i-1]);
    double dy = (yi-y[j-1])/(y[j]-y[j-1]);
    double dz = (zi-z[k-1])/(z[k]-z[k-1]);
    double f1 = f[ i-1 + (j-1)*Nx + (k-1)*Nx*Ny ];
    double f2 = f[ i   + (j-1)*Nx + (k-1)*Nx*Ny ];
    double f3 = f[ i-1 + j*Nx     + (k-1)*Nx*Ny ];
    double f4 = f[ i   + j*Nx     + (k-1)*Nx*Ny ];
    double f5 = f[ i-1 + (j-1)*Nx + k*Nx*Ny     ];
    double f6 = f[ i   + (j-1)*Nx + k*Nx*Ny     ];
    double f7 = f[ i-1 + j*Nx     + k*Nx*Ny     ];
    double f8 = f[ i   + j*Nx     + k*Nx*Ny     ];
    double dx2 = 1.0-dx;
    double dy2 = 1.0-dy;
    double dz2 = 1.0-dz;
    double h0  = (dx*f2 + dx2*f1)*dy2 + (dx*f4 + dx2*f3)*dy;
    double h1  = (dx*f6 + dx2*f5)*dy2 + (dx*f8 + dx2*f7)*dy;
    return h0*dz2 + h1*dz;
}


}

