#git init
#git add .
#git commit -m "first commit"
#git remote add origin https://github.com/smorfopoulou/viral_denovo_pipeline.git
#git push -u origin master

# unset SSH_ASKPASS
# git add scripts/Morfopoulou_pipeline_viral_LEGION_v4.sh
# git commit -m "1) samtools -A flag 2) variant calling PathSeek reqs"
# git push -u origin master

 unset SSH_ASKPASS
 git add scripts/*
 git add exec/QUASR_v7.03/*
 git add exec/picard-tools-1.138/*
 git add README
 git add example_quick/
 git add blast_databases/
 git commit -m "version v5: 1) Add Picard step for removing duplicates based on mapping coord. 2) Call consensus with QUASR 3) Removed EBV flag, use continuous instead for EBV data 4) updated README"
 git push -u origin master
