#!/bin/bash
singularity exec -B /home/users /home/users/pjh/tools/gatk-4.3.0.0.sif gatk "$@"
