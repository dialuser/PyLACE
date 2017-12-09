#!/bin/sh

LOG=$WORK/ims/tacc_jobs.log
IMS_SERVER=xxx.xxx.xxx.xxx
IMS_ACCOUNT=xxx
SSH_CONN="ssh -i $HOME/xxx.pem $IMS_ACCOUNT@$IMS_SERVER"
SCP_CONN="scp -rp -i $HOME/xxx.pem $IMS_ACCOUNT@$IMS_SERVER"
SCP_NO_PATH="scp -rp -i $HOME/xxx.pem"
JOB_SCRIPT=/home/$IMS_ACCOUNT/tacc_jobs/job_manager_psql.sh
