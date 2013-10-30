#!/bin/csh
# -------------------------------------- 
# Script generated  at Mon Nov 14 02:19:57 EST 2011 by the STAR scheduler and submitted with
# cd /star/u/fengzhao/ucla/prog/run11auau/v0_27GeV/batchjobs; condor_submit /star/u/fengzhao/ucla/prog/run11auau/v0_27GeV/batchjobs/schedA702E860D9CCE027669EF26E4CA97157_151_250.condor
# --------------------------------------

onintr  eviction



/bin/echo 'We are starting on node: '`/bin/hostname`


# Preparing environment variables

setenv FILEBASENAME "no_name"
setenv FILELIST "/star/u/fengzhao/ucla/prog/run11auau/v0_27GeV/batchjobs/schedA702E860D9CCE027669EF26E4CA97157_239.list"
setenv FILELIST_ALL "/star/u/fengzhao/ucla/prog/run11auau/v0_27GeV/batchjobs/schedA702E860D9CCE027669EF26E4CA97157.list"
setenv INPUTFILECOUNT "5"
setenv JOBID "A702E860D9CCE027669EF26E4CA97157_239"
setenv JOBINDEX "239"
setenv LOGGING "STD"
setenv REQUESTID "A702E860D9CCE027669EF26E4CA97157"
setenv SCRATCH "/tmp/$USER/$JOBID"
setenv SUBMITATTEMPT "1"
setenv SUBMITTINGDIRECTORY "/star/u/fengzhao/ucla/prog/run11auau/v0_27GeV/batchjobs"
setenv SUBMITTINGNODE "rcas6019.rcf.bnl.gov"
setenv SUBMIT_TIME "2011-11-14 07:19:57"
setenv SUMS_AUTHENTICATED_USER "fengzhao@rhic.bnl.gov"
setenv SUMS_USER "fengzhao"
setenv SUMS_nProcesses "2851"
setenv SUMS_name "null"
setenv USEXROOTD "1"
setenv INPUTFILE0 "root://xrdstar.rcf.bnl.gov:1095//home/starreco/reco/AuAu27_production_2011/FullField/P11id/2011/179/12179038/st_physics_12179038_raw_5030001.MuDst.root"
setenv INPUTFILE1 "root://xrdstar.rcf.bnl.gov:1095//home/starreco/reco/AuAu27_production_2011/FullField/P11id/2011/179/12179059/st_physics_12179059_raw_1040001.MuDst.root"
setenv INPUTFILE2 "root://xrdstar.rcf.bnl.gov:1095//home/starreco/reco/AuAu27_production_2011/FullField/P11id/2011/179/12179059/st_physics_12179059_raw_5020001.MuDst.root"
setenv INPUTFILE3 "root://xrdstar.rcf.bnl.gov:1095//home/starreco/reco/AuAu27_production_2011/FullField/P11id/2011/179/12179094/st_physics_12179094_raw_4010001.MuDst.root"
setenv INPUTFILE4 "root://xrdstar.rcf.bnl.gov:1095//home/starreco/reco/AuAu27_production_2011/FullField/P11id/2011/178/12178020/st_physics_adc_12178020_raw_3520002.MuDst.root"

# This big blobk is used by CONDOR where $HOME in not set

if ( ! $?USER ) then
    echo "USER is not defined"
    set USER=`id | sed "s/).*//" | sed "s/.*(//"`
endif
if ( ! $?HOME ) then
    echo "HOME is not defined"

    if ( -x /usr/bin/getent ) then
	# we have getent, should not be on aix, bsd, Tru64 however
	# will work for Linux
	echo "Using getent method"
	setenv HOME `/usr/bin/getent passwd $USER | /bin/sed 's|.*\:.*\:.*\:.*\:\([^\:]*\):.*|\1|'`
    else 
	set PTEST=`which perl`
	if ( "$PTEST" != "" ) then
	    echo "Using perl method"
	    # we have perl defined, we can get info from there
	    /bin/cat <<EOF >test$$.pl
my(\$user) = getpwuid(\$<);
@items = getpwnam(\$user);
print \$items[7];
EOF
	    setenv HOME `$PTEST test$$.pl` && /bin/rm  -f test$$.pl
	else
	    set CTEST=`which cc`
	    if ( "$CTEST" != "" ) then
		echo "Using C code method"
		# use C code for doing this
		/bin/cat <<EOF >test$$.c
#include <stdio.h>
#include <stdlib.h>
#include <pwd.h>

int main()
{
  struct passwd *info;
  /* assume \$USER is defined */
  const char *user=getenv("USER");
  info = getpwnam(user);
  (void) printf("%s",info->pw_dir);

}

EOF
		$CTEST -o test$$ test$$.c
		/bin/chmod +x test$$
		setenv HOME `./test$$` && /bin/rm -f test$$ test$$.c
	    else
		echo "We have no ways to define HOME and it is not defined"
		exit
	    endif
	endif
    endif
endif

echo "HOME is now $HOME"


# Default value for path if not defined.
if ( ! $?PATH ) then
   setenv PATH /usr/local/bin:/bin:/usr/bin
endif


/usr/bin/test -e $HOME/.cshrc && source $HOME/.cshrc

# Creating the scratch directory, return failure status
set SUMS_try=0
set SUMS_PAD=""

MKDTRY:
setenv SCRATCH "/tmp/$USER$SUMS_PAD/$JOBID"
/bin/mkdir -p $SCRATCH   >& /dev/null
set STS=$status
if (! -d $SCRATCH) then
       #test if porper UID
       set SUMS_IsUID=`(/usr/bin/test -O /tmp/$USER$SUMS_PAD && echo 1) || echo 0`
       if ( $SUMS_try == 0 &&  -e "/tmp/$USER$SUMS_PAD" && ! $SUMS_IsUID) then
            # First try, directory exists but is not owned by $USER
            # Create a different path and try again
            echo "Scheduler:: $SCRATCH not owned by $USER, trying alternative"
            @ try++
            set SUMS_seed=`/bin/date +%j%H%M%S`
            set SUMS_PAD=`expr $SUMS_seed  \* 25973 \% 100`
            goto MKDTRY
        else
		 echo "Scheduler:: Failed to create $SCRATCH on $HOST"
	         exit $STS
        endif

endif

#Note: The default directory in which jobs start has been fix to $SCRATCH
cd $SCRATCH

echo "--------------------------------"
echo 'STAR Unified Meta Scheduler 1.10.4 we are starting in $SCRATCH :' `/bin/pwd`
echo "--------------------------------"
/bin/ls -l

###################################################
# User command BEGIN ----------------------------->

	/star/u/fengzhao/ucla/prog/run11auau/v0_27GeV/runv0.csh $FILELIST $INPUTFILECOUNT 0 output/ $JOBID
   
# <------------------------------User command BEGIN
###################################################

# Copy output files (if any where specified)

# Delete the scratch directory
/bin/rm -fr $SCRATCH




# Function to pass batch system signals to the process, if need
eviction:
set PIDS=`ps  -o "%p" --no-headers --ppid $$`
foreach  process ($PIDS)
  set killit=`ps -o "%p" --no-headers -p $process`
  if ( "x$killit" !~ "x" ) then
     echo "evicting `ps  --no-headers -p$process`"
     kill -15 $process
  endif
end
