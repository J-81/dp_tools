# This script checks all bam files found in the directory and prints a report of checked bams and any issues are logged to bad_bams.txt
find $1 -type f -name "*bam" \
  -execdir bash -c 'echo "Checking bam file: $1"; samtools quickcheck $1 || echo "$1 failed samtools quickcheck" | tee bad_bams.txt' none {} \;