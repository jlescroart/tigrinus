#! /bin/bash

cd /media/labgenoma4/DATAPART4/jonasl/data/tigrinus/generode/subs_dp1.8
for bam in *.bam; do angsd -dofasta 2 -nThreads 3 -doCounts 1 -explode 0 -setMinDepth 5 -minQ 20 -minMapQ 30 -i $bam -r MT: -out /media/labgenoma4/DATAPART4/jonasl/data/tigrinus/generode/mitogenomes/all_from_subs_dp1.8/$bam; done
cd /media/labgenoma4/DATAPART4/jonasl/data/tigrinus/generode/mitogenomes/all_from_subs_dp1.8
gunzip *.gz
rename.ul generode.subs_dp1.8.bam.fa .fasta *.fa

