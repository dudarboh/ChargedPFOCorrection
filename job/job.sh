#!/bin/bash
source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh && \
export MARLIN_DLL=$MARLIN_DLL:/afs/desy.de/user/d/dudarboh/checkVertex/ChargedPFOCorrection/lib/libVertexAnalysis.so && \
rm -rf * && cp ${1} . && \
filename=`ls *.slcio` && \
Marlin /afs/desy.de/user/d/dudarboh/checkVertex/ChargedPFOCorrection/xml/steer.xml --global.LCIOInputFiles="${filename}" && \
mv -f before_refit.root ../output/before_${2}.root && \
mv -f after_refit.root ../output/after_${2}.root && \
rm -f *.slcio
