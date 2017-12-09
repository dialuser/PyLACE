IMS_JOB_ID=$1
if [ -z $IMS_JOB_ID ]
then
    echo "Provide job id for tacc_jobs"
    $SSH_CONN "$JOB_SCRIPT update $IMS_JOB_ID stopped (no job id)"
    exit -1
fi
