//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/base/Utilities.C $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2039 $
// Modified:    $LastChangedDate: 2008-03-11 13:23:52 -0700 (Tue, 11 Mar 2008) $
// Description: Utility functions for error reporting, file manipulation, etc.
//

#include "Utilities.h"

#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/Logger.h"
#include "utils/PIO.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <stdlib.h>
#include <stdexcept>

// Detect the OS and include system dependent headers
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    // Note: windows has not been testeds
    #define USE_WINDOWS
    #include "windows.h"
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

/*
 * Routine to rename a file.
 */
void Utilities::renameFile(const std::string& old_filename,
                           const std::string& new_filename)
{
   AMP_ASSERT(!old_filename.empty());
   AMP_ASSERT(!new_filename.empty());
   rename(old_filename.c_str(), new_filename.c_str());
}

/*
 * Routine to recursively construct directories based on a relative path name.
 */
void Utilities::recursiveMkdir(
   const std::string& path, 
   mode_t mode,
   bool only_node_zero_creates)
{

#ifdef _MSC_VER
   const char seperator = '/';
#define mkdir(path, mode) mkdir(path)
#else
   const char seperator = '/';
#endif

   AMP_MPI comm = AMP_MPI(AMP_COMM_WORLD);
   if ( (!only_node_zero_creates) || (comm.getRank() == 0)) {
      int length = path.length();
      char *path_buf= new char[length+1];
      sprintf(path_buf,"%s",path.c_str());
      struct stat status;
      int pos = length - 1;
   
      /* find part of path that has not yet been created */
      while ( (stat(path_buf,&status) != 0) && (pos >= 0) ) {
   
         /* slide backwards in string until next slash found */
         bool slash_found = false;
         while ( (!slash_found) && (pos >= 0) ) {
           if (path_buf[pos] == seperator) {
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
              path_buf[pos] = seperator;
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
        std::stringstream  stream;
        stream << message << std::endl << "  " << filename << ":  " << line;
        std::cout << stream.str() << std::endl;
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
size_t Utilities::getMemoryUsage()
{
    size_t N_bytes = 0;
    #if defined(USE_LINUX)
        std::string mem;
        std::ifstream proc("/proc/self/status");
        std::string s;
        while(getline(proc, s), !proc.fail()) {
            if(s.substr(0, 6) == "VmSize") {
                mem = s.substr(7);
                break;
            }
        }
        for (size_t i=0; i<mem.size(); i++) {
            if ( mem[i]==9 || mem[i]==10 || mem[i]==12 || mem[i]==13 )
                mem[i] = 32;
        }
        size_t mult = 1;
        if ( mem.find("kB")!=std::string::npos ) {
            mult = 0x400;
            mem.erase(mem.find("kB"),2);
        } else if ( mem.find("MB")!=std::string::npos ) {
            mult = 0x100000;
            mem.erase(mem.find("MB"),2);
        } else if ( mem.find("GB")!=std::string::npos ) {
            mult = 0x40000000;
            mem.erase(mem.find("GB"),2);
        }
        for (size_t i=0; i<mem.size(); i++) {
            if ( ( mem[i]<48 || mem[i]>57 ) && mem[i]!=32 ) {
                printf("Unable to get size from string: %s\n",s.c_str());
                return 0;
            }
        }
        long int N = atol(mem.c_str());
        N_bytes = ((size_t)N)*mult;
    #elif defined(USE_MAC)
        struct task_basic_info t_info;
        mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
        if (KERN_SUCCESS != task_info(mach_task_self(),
                              TASK_BASIC_INFO, (task_info_t)&t_info, 
                              &t_info_count)) {
            return 0;
        }
        N_bytes = t_info.resident_size + t_info.virtual_size;
    #endif
    return N_bytes;
}


//! Function to print the current call stack
std::vector<std::string>  Utilities::getCallStack()
{
    std::vector<std::string>  stack;
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
                stack.push_back(stack_item);
            }
            if ( demangled==NULL ) {
                free(demangled);
                demangled=NULL;
            }
        } 
    #endif
    return stack;
}

}

