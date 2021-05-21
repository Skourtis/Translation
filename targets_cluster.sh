r/bin/bash
#$ -q long-sl7
#$ -l virtual_free=64G
#$ -l h_rt=240:00:00

module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load R/4.0.0-foss-2020a
Rscript /nfs/users2/ssdelci/skourtis/Translation/tar_make_clustermq.R
#added gittoken PAB in the Renvironment which allowed piggyback upload. Add piggyback download
