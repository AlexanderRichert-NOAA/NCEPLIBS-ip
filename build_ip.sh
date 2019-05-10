#!/bin/bash

 (( $# == 0 )) && {
   echo "*** Usage: $0 wcoss|dell|cray|theia|intel_general|gnu_general [debug|build] [[local]install[only]]"
   exit 1
 }

 sys=${1,,}
 [[ $sys == wcoss || $sys == dell || $sys == cray ||\
    $sys == theia || $sys == intel_general || $sys == gnu_general ]] || {
   echo "*** Usage: $0 wcoss|dell|cray|theia|intel_general|gnu_general [debug|build] [[local]install[only]]"
   exit 1
 }
 debg=false
 inst=false
 skip=false
 local=false
 (( $# > 1 )) && {
   [[ ${2,,} == build ]] && debg=false
   [[ ${2,,} == debug ]] && debg=true
   [[ ${2,,} == install ]] && inst=true
   [[ ${2,,} == localinstall ]] && { local=true; inst=true; }
   [[ ${2,,} == installonly ]] && { inst=true; skip=true; }
   [[ ${2,,} == localinstallonly ]] && { local=true; inst=true; skip=true; }
 }
 (( $# > 2 )) && {
   [[ ${3,,} == build ]] && debg=false
   [[ ${3,,} == debug ]] && debg=true
   [[ ${3,,} == install ]] && inst=true
   [[ ${3,,} == localinstall ]] && { local=true; inst=true; }
   [[ ${3,,} == installonly ]] && { inst=true; skip=true; }
   [[ ${3,,} == localinstallonly ]] && { local=true; inst=true; skip=true; }
 }
 if [[ ${sys} == "intel_general" ]]; then
   sys6=${sys:6}
   source ./Conf/Ip_${sys:0:5}_${sys6^}.sh
 elif [[ ${sys} == "gnu_general" ]]; then
   sys4=${sys:4}
   source ./Conf/Ip_${sys:0:3}_${sys4^}.sh
 else
   source ./Conf/Ip_intel_${sys^}.sh
 fi
 [[ -z $IP_VER || -z $IP_LIB4 ]] && {
   echo "??? IP: module/environment not set."
   exit 1
 }

 source ./Conf/Collect_info.sh
 source ./Conf/Gen_cfunction.sh

set -x
 ipLib4=$(basename $IP_LIB4)
 ipLib8=$(basename $IP_LIB8)
 ipLibd=$(basename $IP_LIBd)
 ipInc4=$(basename $IP_INC4)
 ipInc8=$(basename $IP_INC8)
 ipIncd=$(basename $IP_INCd)

#################
 cd src
#################

 $skip || {
#-------------------------------------------------------------------
# Start building libraries
#
 echo
 echo "   ... build (i4/r4) ip library ..."
 echo
   make clean LIB=$ipLib4 MOD=$ipInc4
   mkdir -p $ipInc4
   FFLAGS4="$I4R4 $FFLAGS ${MODPATH}$ipInc4"
   collect_info ip 4 OneLine4 LibInfo4
   ipInfo4=ip_info_and_log4.txt
   $debg && make debug CPPDEFS="-DLSIZE=4" FFLAGS="$FFLAGS4" LIB=$ipLib4 \
                                                             &> $ipInfo4 \
         || make build CPPDEFS="-DLSIZE=4" FFLAGS="$FFLAGS4" LIB=$ipLib4 \
                                                             &> $ipInfo4
   make message MSGSRC="$(gen_cfunction $ipInfo4 OneLine4 LibInfo4)" LIB=$ipLib4

 echo
 echo "   ... build (i8/r8) ip library ..."
 echo
   make clean LIB=$ipLib8 MOD=$ipInc8
   mkdir -p $ipInc8
   FFLAGS8="$I8R8 $FFLAGS ${MODPATH}$ipInc8"
   collect_info ip 8 OneLine8 LibInfo8
   ipInfo8=ip_info_and_log8.txt
   $debg && make debug CPPDEFS="-DLSIZE=8" FFLAGS="$FFLAGS8" LIB=$ipLib8 \
                                                             &> $ipInfo8 \
         || make build CPPDEFS="-DLSIZE=8" FFLAGS="$FFLAGS8" LIB=$ipLib8 \
                                                             &> $ipInfo8
   make message MSGSRC="$(gen_cfunction $ipInfo8 OneLine8 LibInfo8)" LIB=$ipLib8

 echo
 echo "   ... build (i4/r8) ip library ..."
 echo
   make clean LIB=$ipLibd MOD=$ipIncd
   mkdir -p $ipIncd
   FFLAGSd="$I4R8 $FFLAGS ${MODPATH}$ipIncd"
   collect_info ip d OneLined LibInfod
   ipInfod=ip_info_and_logd.txt
   $debg && make debug CPPDEFS="-DLSIZE=D" FFLAGS="$FFLAGSd" LIB=$ipLibd \
                                                             &> $ipInfod \
         || make build CPPDEFS="-DLSIZE=D" FFLAGS="$FFLAGSd" LIB=$ipLibd \
                                                             &> $ipInfod
   make message MSGSRC="$(gen_cfunction $ipInfod OneLined LibInfod)" LIB=$ipLibd
 }

 $inst && {
#
#     Install libraries and source files
#
   $local && {
              LIB_DIR4=..
              LIB_DIR8=..
              LIB_DIRd=..
             } || {
                   LIB_DIR4=$(dirname ${IP_LIB4})
                   LIB_DIR8=$(dirname ${IP_LIB8})
                   LIB_DIRd=$(dirname ${IP_LIBd})
                  }
   [ -d $LIB_DIR4 ] || mkdir -p $LIB_DIR4
   [ -d $LIB_DIR8 ] || mkdir -p $LIB_DIR8
   [ -d $LIB_DIRd ] || mkdir -p $LIB_DIRd
   INCP_DIR4=$(dirname $IP_INC4)
   [ -d $IP_INC4 ] && rm -rf $IP_INC4 || mkdir -p $INCP_DIR4
   INCP_DIR8=$(dirname $IP_INC8)
   [ -d $IP_INC8 ] && rm -rf $IP_INC8 || mkdir -p $INCP_DIR8
   INCP_DIRd=$(dirname $IP_INCd)
   [ -d $IP_INCd ] && rm -rf $IP_INCd || mkdir -p $INCP_DIRd
   SRC_DIR=$IP_SRC
   $local && SRC_DIR=
   [ -d $SRC_DIR ] || mkdir -p $SRC_DIR
   make clean LIB=
   make install LIB=$ipLib4 MOD=$ipInc4 \
                LIB_DIR=$LIB_DIR4 INC_DIR=$INCP_DIR4 SRC_DIR=
   make install LIB=$ipLib8 MOD=$ipInc8 \
                LIB_DIR=$LIB_DIR8 INC_DIR=$INCP_DIR8 SRC_DIR=
   make install LIB=$ipLibd MOD=$ipIncd \
                LIB_DIR=$LIB_DIRd INC_DIR=$INCP_DIRd SRC_DIR=$SRC_DIR
 }
