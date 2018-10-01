#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
exe_name=$0
exe_dir=`dirname "$0"`
echo "------------------------------------------"
if [ "x$1" = "x" ]; then
  echo Usage:
  echo    $0 hierBAPS_CMD args
else
  echo hierBAPS directory:  $exe_dir
  echo Setting up environment variables
  MCRROOT="/home/lcheng/MATLAB/MATLAB_Compiler_Runtime/v84"; #set your MCR installation directory here.
  
  echo ---
  LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
  export LD_LIBRARY_PATH;
  echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
  
  echo "Finnished."
  echo "------------------------------------------\n"  

  echo "${exe_dir}"/"$@"
  "${exe_dir}"/"$@"
  
fi
exit

