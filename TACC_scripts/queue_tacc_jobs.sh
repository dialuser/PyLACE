#/bin/sh

source ./inc_ims.sh

if [ -f $LOG ]
then
    LOG_SIZE=$(du -m $LOG | cut -f 1)
    if [ $LOG_SIZE -ge 1024 ]
    then
        CURR_TIME=`date`
        echo "$LOG is reset at $CURR_TIME" > $LOG
    fi
fi

echo "Try to pull job remotely from ims at `date`" >> $LOG
CMD_GET_JOB="$SSH_CONN $JOB_SCRIPT getjob"
MSG=`$SSH_CONN "$JOB_SCRIPT getjob"`
if [ -z "$MSG" ]
then
    echo "No job found from ims." >> $LOG
else 
    ID=`echo $MSG | cut -d',' -f1`
	CMD=`echo $MSG | cut -d',' -f2`
    IMS_PATH=`echo $MSG | cut -d',' -f3`
    TACC_PATH=`echo $MSG | cut -d',' -f4`
	echo "Job $ID is found to run $CMD" >> $LOG 
    $SSH_CONN "$JOB_SCRIPT update $ID received"
    PARENT_TACC_PATH=`dirname $TACC_PATH`
    $SCP_CONN:$IMS_PATH $PARENT_TACC_PATH/
    #echo $IMS_PATH
    #echo $TACC_PATH
    #echo $MSG
    #echo $CMD
    eval $CMD $ID
    $SSH_CONN "$JOB_SCRIPT update $ID queued"
	echo "Job $ID is queued at `date`" >> $LOG
fi


