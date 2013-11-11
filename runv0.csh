#!/usr/local/bin/tcsh
setenv STAR_LEVEL SL11d
if ( -e $GROUP_DIR/group_env.csh ) then
        source $GROUP_DIR/group_env.csh
endif

cd /star/u/lwen1990/ucla/v0_200GeV_Liwen/
root4star -b -q Macro.C\(\"$1\",$2,\"$3\",\"$4\"\)
