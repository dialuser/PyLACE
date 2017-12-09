#!/bin/bash

###$SQL settings####

USER=xxx
TBL=xxx
DB=xxx
SQL=psql -U $USER -h localhost $DB -c 

get_queued_running_jobs()
{
    QUERY="select id,tacc_job_id from $TBL where status=\"queued\" or status=\"running\""
    $SQL -u $USER -N -s -e "$QUERY" $DB
}

#########
## pull job information
##########
getJob()
{
GETJOB="select id, command from $TBL where status=\"submitted\""

RESULT=`$SQL -B -N -s -e "$GETJOB" $DB | tr '\t' ','`

ID=`echo $RESULT |  cut -d',' -f1`
CMD=`echo $RESULT | cut -d',' -f2`

#echo Running Query: $GETJOB 

if [ -z "$ID" ]
then
   echo ""
else 
   echo $ID, $CMD
fi

INSERT_JOB_LOG_FN="update $TBL set job_log_fn=\"JOB_LOG_$ID\" where id=\"$ID\""
$SQL -u $USER -N -s -e "$INSERT_JOB_LOG_FN" $DB

}

#########
## update function
##take two parameters
#########

update()
{

if [ -z "$1" ]
then
  echo "Please sepecify job ID, status"
  exit -1
fi
  
ID=$1
STATUS=$2

QUERY="update $TBL set status=\"$STATUS\", ${STATUS}_time=CURRENT_TIMESTAMP where id=\"$ID\""
#echo $QUERY
$SQL -u $USER -N -s -e "$QUERY" $DB

}

#########
## update job info function
##take two parameters
#########

update_job_info()
{

if [ -z "$1" ]
then
  echo "Please sepecify job ID, slurm job id, slurm output filename, job log filename"
  exit -1
fi
  
ID=$1
TACC_JOB_ID=$2
TACC_OUTPUT_FN=$3
JOB_LOG_FN=$4

QUERY="update $TBL set tacc_job_id=$TACC_JOB_ID,tacc_output_fn=\"$TACC_OUTPUT_FN\" where id=$ID"
#echo $QUERY
$SQL -u $USER -N -s -e "$QUERY" $DB

}

update_queued_running_jobs()
{
    if [ -z "$1" ]
    then
        echo "Please specify job id, status, running time"
        exit -1
    fi
    ID=$1
    STATUS=$2
    RUNNING_TIME=$3
    QUERY="update $TBL set status=lower(\"$STATUS\"), tacc_running_time=\"$RUNNING_TIME\" where id=$ID"
#echo $QUERY
    $SQL -u $USER -N -s -e "$QUERY" $DB
}

usage()
{
    echo $"USAGE: $1 {getjob|update id status|update_job_info id tacc_id tacc_output_fn} args"
}


#######
#main function
######


if [ -z "$1" ]
then
#echo $"USAGE: $0 {getjob|update id status} args"
  usage $0
  exit -1
fi  

case $1 in 
   getjob)
      getJob $2
      ;;
   update)
      update "$2" "$3"
      ;;
   update_job_info)
      update_job_info "$2" "$3" "$4" "$5"
      ;; 
   get_queued_running_jobs)
      get_queued_running_jobs
      ;;
   update_queued_running_jobs)
      update_queued_running_jobs "$2" "$3" "$4"
      ;;

   *)
 #     echo $"USAGE: $0 {getjob|update id status} args"
      usage $0
      exit 1
esac 

