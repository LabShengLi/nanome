#!/bin/bash

for i in {1..22} X Y; do
	echo "sanity check chr${i}"
	bash sanity_check_deepmod_submit.sh chr${i}
done

exit 0

#bash sanity_check_deepmod_submit.sh chr6
#bash sanity_check_deepmod_submit.sh chr7
#bash sanity_check_deepmod_submit.sh chr9
#bash sanity_check_deepmod_submit.sh chr10
#bash sanity_check_deepmod_submit.sh chr11
#
#exit 0

