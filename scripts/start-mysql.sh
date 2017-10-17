#!/bin/bash
/pipelines/scripts/cromwell.scala mysql 5.7 /pipelines /pipelines
sleep 5
echo "started mysql service for cromwell"


