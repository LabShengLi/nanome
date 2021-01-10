#!/bin/bash

set -x

export abc="Yang Liu"

ret=$(sbatch --export=ALL,abc1=Tombo --wait test-script.sh)
echo ${ret}

