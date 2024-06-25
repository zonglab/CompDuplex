To run this pipeline, run the following command:  
bash SCRIPTs/s1_CompDup_00_submit_job.sh <working_dir> <fq___prefix> <normal_bulk_path>  
<working_dir>: ideally, the absolute path to where you want to store the analysis results  
<fq___prefix>: path and prefix to the sequenced fastq files, e.g.: /PATH/TO/PREFIX corresponding to /PATH/TO/PREFIX_R1[_001].fastq.gz  
<normal_bulk>: path to the 30X WGS matched normal bulk results, must include: 05_add_RG.bam, 15_fin_hetero.vcf, and 15_fin___homo.vcf  

A test dataset is included surrounding hg19 chr14:100,000,000-108,000,000.

package requirements:
cutadapt 4.8  
     bwa 0.7.13-r1126  
  picard 3.1.1  
samtools 1.12  
bcftools 1.3  
   seqtk 1.4-r122  
 SnpSift 4.1k  

python packages:  
   pysam 0.22.1  
   numpy 1.26.4  
  pandas 2.2.2  
