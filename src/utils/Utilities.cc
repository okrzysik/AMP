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
    #include<mach/mach.h>
#elif defined(__linux) || defined(__unix) || defined(__posix)
    #define USE_LINUX
    #include <signal.h>
    #include <execinfo.h>
    #include <cxxabi.h>
    #include <dlfcn.h>
#else
    #error Unknown OS
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
            printf("   %s",stack[i].c_str());
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


// Function to get the memory usage
// Note: this function should be thread-safe
#if defined(USE_LINUX)
    // Get the page size on linux
    size_t page_size = static_cast<size_t>(sysconf(_SC_PAGESIZE));
    // Open the /proc/self/status so we can cache the access
    static void fidDeleter(FILE* fid) { fclose(fid); };
    static boost::shared_ptr<FILE> proc_fid;
    static size_t N_bytes_initialization = Utilities::getMemoryUsage();
#endif
size_t Utilities::getMemoryUsage()
{
    size_t N_bytes = 0;
    #if defined(USE_LINUX)
        /* Use /proc/self/statm
         *     Provides information about memory usage, measured in pages.
         *     The columns are:
         *         size       (1) total program size (same as VmSize in /proc/[pid]/status)
         *         resident   (2) resident set size  (same as VmRSS in /proc/[pid]/status)
         *         share      (3) shared pages (i.e., backed by a file)
         *         text       (4) text (code)
         *         lib        (5) library (unused in Linux 2.6)
         *         data       (6) data + stack
         *         dt         (7) dirty pages (unused in Linux 2.6)
         */
        if ( proc_fid==NULL ) 
            proc_fid = boost::shared_ptr<FILE>(fopen("/proc/self/statm","rb"),fidDeleter);
        FILE* fid = proc_fid.get();
        if (fid==NULL) { return 0; }
        char data_text[128];
        flockfile(fid);     // lock file for thread safety
        int size = fread(data_text,1,127,fid);  // read the data
        rewind(fid);        // rewind the file so it is availible for the next read
        funlockfile(fid);   // unlock file
        char* tmp = data_text;
        long int data[10];
        int i = 0;
        while ( tmp-data_text < size-1 ) {
            data[i] = strtol(tmp,&tmp,10);
            i++;
        }
        if ( i!=7 ) {
            printf("Format for /proc/self/statm changed\n");
            return 0;
        }
        N_bytes = page_size*static_cast<size_t>(data[5]);
    #elif defined(USE_MAC)
        struct task_basic_info t_info;
        mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
        if (KERN_SUCCESS != task_info(mach_task_self(),
                              TASK_BASIC_INFO, (task_info_t)&t_info, 
                              &t_info_count)) {
            return 0;
        }
        N_bytes = t_info.resident_size;
    #elif defined(USE_WINDOWS)
        PROCESS_MEMORY_COUNTERS memCounter;
        GetProcessMemoryInfo( GetCurrentProcess(), &memCounter, sizeof(memCounter) );
        N_bytes = memCounter.WorkingSetSize;
    #endif
    return N_bytes;
}


//! Function to print the current call stack
std::vector<std::string>  Utilities::getCallStack()
{
    std::vector<std::string>  stack_list;
    #if defined(USE_LINUX) || defined(USE_MAC)
        void *trace[100];
        Dl_info dlinfo;
        int status;
        const char *symname;
        char *demangled=NULL;
        int trace_size = backtrace(trace,100);
        for (int i=0; i<trace_size; ++i) {  
            if(!dladdr(trace[i], &dlinfo))
                continue;
            symname = dlinfo.dli_sname;
            demangled = abi::__cxa_demangle(symname, NULL, 0, &status);
            if(status == 0 && demangled)
                symname = demangled;
            std::string object = std::string(dlinfo.dli_fname);
            std::string function = "";
            if ( symname!=NULL )
                function = std::string(symname);
            if ( i!=0 ) {  // Skip the current function
                std::string stack_item = object + ":   " + function + "\n";
                //stack_item = "object: " + object + "\n";
                //stack_item += "function: " + function + "\n";
                stack_list.push_back(stack_item);
            }
            if ( demangled==NULL ) {
                free(demangled);
                demangled=NULL;
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
            bool rtn = ::StackWalk64( MachineType, ::GetCurrentProcess(), ::GetCurrentThread(), &lFrameStack, MachineType == IMAGE_FILE_MACHINE_I386 ? 0 : &lContext,
                NULL, &::SymFunctionTableAccess64, &::SymGetModuleBase64, NULL );
            if( !rtn )
                break;
            if( lFrameStack.AddrPC.Offset == 0 )
                break;
            ::MEMORY_BASIC_INFORMATION lInfoMemory;
            ::VirtualQuery( ( ::PVOID )lFrameStack.AddrPC.Offset, &lInfoMemory, sizeof( lInfoMemory ) );
            ::DWORD64 lBaseAllocation = reinterpret_cast< ::DWORD64 >( lInfoMemory.AllocationBase );
            ::TCHAR lNameModule[ 1024 ];
            ::GetModuleFileName( reinterpret_cast< ::HMODULE >( lBaseAllocation ), lNameModule, 1024 );
            PIMAGE_DOS_HEADER lHeaderDOS = reinterpret_cast<PIMAGE_DOS_HEADER>( lBaseAllocation );
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

            /*DWORD lineDisplacement = 0;
            IMAGEHLP_LINE64  line;
            bool bLine = SymGetLineFromAddr64( ::GetCurrentProcess(), lFrameStack.AddrPC.Offset, &lineDisplacement, &line );
            PDWORD64 symDisplacement = 0;
            enum { emMaxNameLength = 512 };
            union {
                SYMBOL_INFO symb;
                BYTE symbolBuffer[ sizeof(SYMBOL_INFO) + emMaxNameLength ];
            } u;
            PSYMBOL_INFO pSymbol = & u.symb;
            char buf[100];
            DWORD64 pFrame = lFrameStack.AddrFrame.Offset;
            if ( SymFromAddr( ::GetCurrentProcess(), lFrameStack.AddrPC.Offset,
                                symDisplacement, pSymbol) )
            {
                if( bLine ) {
                    sprintf( buf, "   %s() line %d\n",
                        lFrameStack.AddrPC.Offset, pFrame,
                        pSymbol->Name, line.LineNumber );
                } else {
                    sprintf( buf, "  %s() + %X\n",
                        lFrameStack.AddrPC.Offset, pFrame,
                        pSymbol->Name, symDisplacement );
                }

            }
            else    // No symbol found.  Print out the logical address instead.
            {
                DWORD err = GetLastError();
                //FIXED_ARRAY( szModule , TCHAR, MAX_PATH );
                //char szModule = '\0';
                DWORD section = 0, offset = 0;

                //GetLogicalAddress(  (PVOID)lFrameStack.AddrPC.Offset,
                //                    szModule, sizeof(szModule), section, offset );

                char szModule=0;
                sprintf( buf, "  %04X:%08X %s (err = %d)\n",
                    lFrameStack.AddrPC.Offset, pFrame,
                    section, offset, szModule, err );
            }*/

/*            // Save line
            size_t l = strlen(buf);
            if( i_line >= m_Levels || i_buf + l >= m_Bytes ) {
                // We have saved all of the stack we can save
                break;
            }
            buf[ l - 1 ] = '\0';    // Remove trailing '\n'
            char * s = & m_Buffer[ i_buf ];
            m_Lines[ i_line++ ] = s;
            strncpy( s, buf, l );
            i_buf += l;*/
        }
    #endif
    return stack_list;
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

