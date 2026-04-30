#!/bin/bash


echo "##############################################"
echo "####                                     ####"
echo "####            DM EVENT PULL            ####"
echo "####                                     ####"
echo "#############################################"

echo "     "
echo " ... Checking PSQL log for EPIC ID: "$1
echo " "
echo "  id   | logid |     instance      |   day    "
echo "-------+-------+-------------------+----------"


#psql -d postgres -U eewuser -c "select  id,logid,instance,day  from dmbkprod1.event where id IN ($1) and ver=0 and system='epic'"  > testsept.txt

OUTP=`psql -d postgres -U eewuser -c "select  id,logid,instance,day  from dmbkprod1.event where id IN ($1) and ver=0 and system='epic'" | sed  -n '3 p'`

echo $OUTP 
echo "-------+-------+-------------------+----------"

#######################################################
input_id=`echo $OUTP | awk -F'|' '{print $1}'`
logid=`echo $OUTP | awk -F'|' '{print $2}'`
instance=`echo $OUTP | awk -F'|' '{print $3}'`
day=`echo $OUTP | awk -F'|' '{print $4}'`
#####################################################



python pull_logs.py $instance $day $input_id $logid



