#!/usr/bin/env nextflow

bams =  Channel.fromPath("../preprocessing/realigned-bams/*.bam")}
bais =  Channel.fromPath("../preprocessing/realigned-bams/*.bai")}

regionTasks = Channel
  .from(file(params.regionTasksFreebayes).readLines())
  .map {line ->
    list = line.split(",")
    task = list[0]
    flag = list[1]
    [ task, flag ]
}


process Freebayes1 {
  publishDir "outputs/regionVCFs"
  tag "$regionTask"
  
  cpus { 1 }
  memory { 32.GB * task.attempt }
  time { 2.d }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'

  input:
  file(bams) from bamFiles.first()
  file(bais) from baiFiles.first()
  set regionTask, regions from regionTasks
  file params.genome
  file params.genomeIndex

  output:
  set regionTask, file("region_${regionTask}.vcf") into region_VCFs

  script:
  input_bams = bams.collect{"-b $it"}.join(' ')
  
  """
  freebayes --use-best-n-alleles 4 --min-coverage 5 \
    -f $params.genome \
    $input_bams \
    $regions > region_${regionTask}.vcf
  """
}
  
all_VCFs = region_VCFs.toList()

process ConcatenateVCFs {
  publishDir "outputs"

  memory { 64.GB }
  time { 24.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 1
  maxErrors '-1'

  input:
  file(vcfs) from region_VCFs

  output:
  set file("freebayes-min5.raw.snps.indels.vcf") into finalVCF

  script:
  input_vcfs = vcfs.collect{"$it"}.join(' ')

  """
  vcf-concat ${input_vcfs} > freebayes-min5.raw.snps.indels.vcf

  """
}
