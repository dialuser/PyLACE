#/bin/sh

source inc_ims.sh

GET_JOB_ID=(`$SSH_CONN $JOB_SCRIPT get_running_jobs`)
if [ -z "$GET_JOB_ID" ]
then
    echo "No running job found"
else
    NUM=${#GET_JOB_ID[@]}
    for (( i=0; i<$((NUM/2)); i++));
    do
        JOB_ID=${GET_JOB_ID[i*2]}
        TACC_JOB_ID=${GET_JOB_ID[i*2+1]}
        OUTPUT=`squeue --job $TACC_JOB_ID -o "%T %M" -h 2> /dev/null`
        if [ -n "$OUTPUT" ]
        then
            STATUS=`echo $OUTPUT | cut -d' ' -f1`
            RUNNING_TIME=`echo $OUTPUT | cut -d' ' -f2`
            #echo $STATUS
            #echo $RUNNING_TIME
            $SSH_CONN "$JOB_SCRIPT update_running_jobs $JOB_ID $STATUS $RUNNING_TIME"
        fi
    done
fi

