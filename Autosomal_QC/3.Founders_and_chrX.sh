#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G

#########
## make sure the .fam file in all the plink files contains information of the actual pedigrees so plink can identify the founders and the females
ml RPlus 
##input and output variables
inputdir="" ### path to the most QCed set of autosomal chromosome files (family must be completely corrected and samples must match with the ones in X_QC/2_CR_high/chr_X.fam)
GeneralQCDir=""  ### working directory for general QC (required to have ./X_QC/2_CR_high/, ./plots/ and ./X_QC/3_MAF_HWE/ )
codedir="" ### directory with the scripts (requeired to have sub_founders_XchrHW.R )
pedigree_fam="" ### final .fam with pedigree information (parents and sex for each sample)
plinkmod="PLINK/1.9-beta6-20190617"    ### plink module, in case it is different from "PLINK/1.9-beta6-20190617"
##### 
# point for your parameters file.
if [ -z ${parameters_file+x} ]; then
  echo "parameter_file unset. Using default parameters..."
else
  source ${parameters_file}
  echo "parameters file set to: '${parameters_file}'"
fi

if [ ! -z ${parameters_file+x} ]; then
  if [ -f "${GeneralQCDir}/parameters_file_founders.sh" ]; then
    echo "parameter file already present and will be rewritten"
      cp ${parameters_file} "${GeneralQCDir}/parameters_file_founders.sh"
  fi
  cp ${parameters_file} "${GeneralQCDir}/parameters_file_founders.sh"
fi

######### separate founders and create H-W for X chromosome founders stats. will create the folder ${GeneralQCDir}\7_Founder_stats
### create plink files and call_rate stats for individuals and SNPs
Rscript ${codedir}/sub_founders_XchrHW.R -p ${inputdir} \
                                         -px ${GeneralQCDir}/X_QC/2_CR_high/chr_X \
                                         -pex ${plinkmod} \
                                         -ped ${pedigree_fam} \
                                         -out ${GeneralQCDir} \
### take the ${GeneralQCDir}\7_Founder_stats results to finalize cthe QC on both X and Autosomes
##remove SNPs flagged as ouliert for H-W in chromosome X
ml ${plinkmod} 

plink \
--bfile ${GeneralQCDir}/X_QC/2_CR_high/chr_X \
--exclude ${GeneralQCDir}\7_Founder_stats/founder.final.list.excluded.snps \
--make-bed \
--out ${GeneralQCDir}/X_QC/3_MAF_HWE/chr_X

#### create final autosomes
mkdir -p ${GeneralQCDir}/8_final_QCed_autosomes_X
cp ${GeneralQCDir}/X_QC/3_MAF_HWE/chr_X.* ${GeneralQCDir}/ ## copy the QCed X chromosome

for chr in {1..22} "XY"
do
plink \
--bfile ${inputdir}/chr_${chr} \
--exclude ${GeneralQCDir}\7_Founder_stats/founder.final.list.excluded.snps \
--make-bed \
--out ${GeneralQCDir}/8_final_QCed_autosomes_X/chr_${chr}
done

#################
