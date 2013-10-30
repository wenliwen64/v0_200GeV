#!/usr/local/bin/tcsh
setenv STAR_LEVEL SL11d
if ( -e $GROUP_DIR/group_env.csh ) then
        source $GROUP_DIR/group_env.csh
endif

cd /star/u/fengzhao/ucla/prog/run11auau/v0_27GeV
root4star -b -q reconstruct_v0.C\(\"$1\",$2,$3,\"$4\",\"$5\"\)
