  ########################################################################################
  # Usage function... Make sure we are called safely before proceeeding...               # 
  ########################################################################################
  # Check the arguments to make sure we got the start date... 
  if [ "$#" -ne 1 ]
  then
    echo "Warning, no input start time was specified. Please see the invocation usage below for an example."
    echo "Example Usage: $0 20170326"
    echo "Aborting. Please try again after correcting usage."
    exit -1
  fi

 
  echo "*******************************************************************************"
  echo "Starting DeepCurrent execution..."
  echo "*******************************************************************************"

  # Show how the command was launched for debugging purposes...
  echo "Run was invoked with: \"$*\""

  # This path addition should now handled by /etc/profile.d/netCDF.sh... 
  #export $LD_LIBRARY_PATH=/usr/local/netcdf/lib:$LD_LIBRARY_PATH
  

  # Set the number of processes to use for this invocation (per-node)
  nproc=30

  # Specify the input file configuring the model
  INPFIL="EFDC.INP"
 
  # Set the start date to the argument to this script... 
  #start_date="20170328"   # start day of form YYYYMMDD (assume UTC)
  start_date=$1

  # Simulation duration in days, we do a 24 hour simulation by default...
  exectime="1"       # simulation time in days

  echo "Editing the DeepCurrent input files with correct start date..."
  #FIXME: This should be rewritten with a template file then an eye-catcher for sed changes 

  # edit EFDC input files with start date and run_time
  awk -v exectime="$exectime" 'NR==175{$1=exectime}1;' "$INPFIL" > foobar    && mv foobar "$INPFIL"; 
  awk -v start_date="$start_date" 'NR==191{$2=start_date}1;' "$INPFIL" > foobar   && mv foobar "$INPFIL";
  echo "File configuration complete."

  # Specify path here to the DeepThunder outputs since it is not static
  # Generally it depends on the start date
  # DT configured to write files to output directory of the form 20161115-00z
  #
  pathdt="/jp2/tjpwm/runoff_test/Runoff/20160501-00z-Fearghal/ncOutput_MERGED"
#  sed -i "2s/.*/$pathdt/"  'DTCOUPLE.INP'
  echo "pathdt="$pathdt
  awk -v  pathdt="$pathdt" 'NR==2{$0=pathdt}1;' "DTCOUPLE.INP" > foobar && mv foobar "DTCOUPLE.INP"
  
  # Ceate directory for outputfiles if doesn't already exist
  # TODO Need to find way to pass this to EFDC without having to read
 
  # For some weird reason jpnx-io1 sees root directory as /jp2
  # while jp-mapping-bz sees it as /jp2-jpnx
  # Hence, create the directory dirwrite for DC to write to and 
  # dirout to create a symbolic link
 

  
  # Actually invoke MPI to run the job on the distributed nodes...  
  echo "Begin DC simulation "
  # Run under valgrind...
  mpirun -np $nproc  ./EFDC
 echo "*******************************************************************************"
  echo " DeepCurrent execution Ended ..."
  echo "*******************************************************************************"
