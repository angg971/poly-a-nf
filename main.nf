// Declare syntax version
nextflow.enable.dsl=2

/* 
 * pipeline input parameters 
 */


params.reads = "$projectDir/data/*_{1,2}.fastq"
params.ref_genome = "$projectDir/data/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa"
params.gtf_file = "$projectDir/data/Glycine_max.Glycine_max_v2.1.55.gtf"
params.outdir = "results"



log.info """\

      =================================================
            R N A S E Q - N F   W O R K F L O W  

            F O R  T H E  D E T E C T I O N   O F

            P O L Y A D E N Y L A T E D  R E A D S  
      =================================================
        reads            : ${params.reads}
        reference genome : ${params.ref_genome}
        gtf              : ${params.gtf_file}
        outdir           : ${params.outdir}
        """
        .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the reference genome file
 */

process INDEX {

  input:
  path ref_genome

  output:
  path "index.*"

  script:
  """
  bowtie2-build $ref_genome index
  """
}

/*
 * define the `GTF2BED` process that creates a BED file
 * given the params.gtf_file as input. 
 */
 
process GTF2BED {

    input:
    path gtf

    output:
    path "${gtf.baseName}.bed"

    script:
    """
    gtf2bed --gtf "${gtf}" --bed "${gtf.baseName}.bed"
    """
}


/*
 * define the `TRIM` process that creates trimmed fastq files
 * given the reads fastq file. 
 */
 
process TRIM {

    tag { pair_id }

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path(trimmed_reads), emit: trimmed_fastqs
    path "${pair_id}/*.atria.log*", emit: logs

    script:
    def (r1, r2) = reads

    trimmed_reads = reads.collect {
        "${pair_id}/${it.baseName}.atria.fastq"
    }

    """
    atria \\
        -r "${r1}" \\
        -R "${r2}" \\
        -o "${pair_id}" \\
        -t 1 \\
        -l 12 \\
        -q 30
    """
}

/*
 * define the `SAM_FOR_STRAND` process that creates a SAM file
 * to be used by `STRAND_CHECK` process to determine 
 * the strandedness of each fastq file pair. 
 */
 
process SAM_FOR_STRAND {

    input:
    path index_files, stageAs: "bt2_index/*"
    tuple val(pair_id), path(reads)

    output:
    path "${pair_id}.sam"

    script:
    
    def idx = index_files[0].getBaseName(2)    

    """
    bowtie2 -x "${idx}" -1 "${reads[0]}" -2 "${reads[1]}" -S "${pair_id}.sam"
    """
}



/*
 * the `CHECK_STRANDEDNESS` process takes SAM file 
 * and a BED file
 * and outputs a a txt file containing information about strandedness. 
 */
 
process CHECK_STRANDEDNESS {

    tag { pair_id }

    input:
    path sam
    path bed
    tuple val(pair_id), path(reads)

    output:
    path "${pair_id}_strand_check.txt"

    script:
    """
    infer_experiment.py -i "${sam}" -r "${bed}" > "${pair_id}_strand_check.txt"
    """
}

/*
 * the `FINDTAIL` process takes "${pair_id}_strand_check.txt" file  
 * and executes commands according to the reads strandedness. 
 * The output generated is FASTA files.  
 */
 
process FINDTAIL {

    tag { pair_id }

    debug true

    input:
    path strand_check_file
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("tail_${pair_id}_{1,2}{,.pr,.sl}.fasta")
    
    script:
    """
#!/usr/bin/env python3

import subprocess
import pandas as pd
import os


result = pd.read_csv("${strand_check_file}", sep="\\r\\n", header=None, engine='python')

failed = float(result.iloc[1,0].replace('Fraction of reads failed to determine: ', ''))
fwd = float(result.iloc[2,0].replace('Fraction of reads explained by "1++,1--,2+-,2-+": ', ''))
rev = float(result.iloc[3,0].replace('Fraction of reads explained by "1+-,1-+,2++,2--": ', ''))
fwd_percent = fwd/(fwd+rev)
rev_percent = rev/(fwd+rev)

if float(result.iloc[1,0].replace('Fraction of reads failed to determine: ', '')) > 0.50:
    cmd1 = "findtail_v1.01 --input_file " + "${reads[0]}" + " --seqlength 500 --endgap 2 --taillength 10 --identity 10 --ptype A --stype T --output_format fasta --output_type sl --d > " + "tail_${pair_id}_1.sl.fasta"
    cmd2 = "findtail_v1.01 --input_file " + "${reads[1]}" + " --seqlength 500 --endgap 2 --taillength 10 --identity 10 --ptype A --stype T --output_format fasta --output_type sl --d > " + "tail_${pair_id}_2.sl.fasta"
    print(cmd1)
    subprocess.call(cmd1, shell=True)
    print(cmd2)
    subprocess.call(cmd2, shell=True)

if fwd_percent > 0.9:
    #Over 90% of reads explained by "1++,1--,2+-,2-+
    #Data is likely FR/fr-secondstrand
    cmd1 = "findtail_v1.01 --input_file " + "${reads[0]}" + " --seqlength 500 --endgap 2 --taillength 10 --identity 10 --ptype A --stype T --output_format fasta --output_type pr --d > " + "tail_${pair_id}_1.pr.fasta"
    cmd2 = "findtail_v1.01 --input_file " + "${reads[1]}" + " --seqlength 500 --endgap 2 --taillength 10 --identity 10 --ptype A --stype T --output_format fasta --output_type pr --d > " + "tail_${pair_id}_2.pr.fasta"
    print(cmd1)
    subprocess.call(cmd1, shell=True)
    print(cmd2)
    subprocess.call(cmd2, shell=True)


elif rev_percent > 0.9:
    # Over 90% of reads explained by "1+-,1-+,2++,2--
    # Data is likely RF/fr-firststrand
    cmd1 = "findtail_v1.01 --input_file " + "${reads[0]}" + " --seqlength 500 --endgap 2 --taillength 10 --identity 10 --ptype A --stype T --output_format fasta --output_type sl --d > " + "tail_${pair_id}_1.sl.fasta"
    cmd2 = "findtail_v1.01 --input_file " + "${reads[1]}" + " --seqlength 500 --endgap 2 --taillength 10 --identity 10 --ptype A --stype T --output_format fasta --output_type sl --d > " + "tail_${pair_id}_2.sl.fasta"
    print(cmd1)
    subprocess.call(cmd1, shell=True)
    print(cmd2)
    subprocess.call(cmd2, shell=True)

elif max(fwd_percent, rev_percent) < 0.6:
    #Under 60% of reads explained by one direction
    #Data is likely unstranded
    cmd1 = "findtail_v1.01 --input_file " + "${reads[0]}" + " --seqlength 500 --endgap 2 --taillength 10 --identity 10 --ptype A --stype T --output_format fasta --output_type sl --d > " + "tail_${pair_id}_1.sl.fasta"
    cmd2 = "findtail_v1.01 --input_file " + "${reads[0]}" + " --seqlength 500 --endgap 2 --taillength 10 --identity 10 --ptype A --stype T --output_format fasta --output_type pr --d > " + "tail_${pair_id}_1.pr.fasta"
    cmd3 = "findtail_v1.01 --input_file " + "${reads[1]}" + " --seqlength 500 --endgap 2 --taillength 10 --identity 10 --ptype A --stype T --output_format fasta --output_type sl --d > " + "tail_${pair_id}_2.sl.fasta"
    cmd4 = "findtail_v1.01 --input_file " + "${reads[1]}" + " --seqlength 500 --endgap 2 --taillength 10 --identity 10 --ptype A --stype T --output_format fasta --output_type pr --d > " + "tail_${pair_id}_2.pr.fasta"
    print(cmd1)
    subprocess.call(cmd1, shell=True)
    print(cmd2)
    subprocess.call(cmd2, shell=True)
    print(cmd3)
    subprocess.call(cmd3, shell=True)
    print(cmd4)
    subprocess.call(cmd4, shell=True)

else:
    #if strand couldnt be detected
    # we assume it is first-stranded
    # since most illumina reads are
    cmd1 = "findtail_v1.01 --input_file " + "${reads[0]}" + " --seqlength 500 --endgap 2 --taillength 10 --identity 10 --ptype A --stype T --output_format fasta --output_type sl --d > " + "tail_${pair_id}_1.sl.fasta"
    cmd2 = "findtail_v1.01 --input_file " + "${reads[1]}" + " --seqlength 500 --endgap 2 --taillength 10 --identity 10 --ptype A --stype T --output_format fasta --output_type sl --d > " + "tail_${pair_id}_2.sl.fasta"
    print(cmd1)
    subprocess.call(cmd1, shell=True)
    print(cmd2)
    subprocess.call(cmd2, shell=True)
    """
}

/*
 * define the `ALIGN` process that outputs unmapped files in FASTA format
 * It will take the FASTA files creates by the `FINDTAIL` process 
 * and try to align them INDIVIDUALLY.
 */

 
process ALIGN {

    tag { pair_id }

    input:
    path index_files, stageAs: "bt2_index/*"
    tuple val(pair_id), path(fasta_file)

    output:
    tuple val(pair_id), path("${fasta_file.baseName}.unmapped.fasta")

    script:
    def idx = index_files[0].getBaseName(2)    

    """
    bowtie2 -x "${idx}" -f "${fasta_file}" --un "${fasta_file.baseName}.unmapped.fasta"
    """
}

/*
 * define the `TRIM_TAILS` process that inputs unmapped files in FASTA format
 * and outputs FASTA file with trimmed tails.
 */
 
process TRIM_TAILS {

    tag { pair_id }

    input:
    tuple val(pair_id), path(unmapped_fasta_file)

    output:
    tuple val(pair_id), path ("${unmapped_fasta_file.baseName}_tail.trimmed.fasta")

    script:
    """
    tailTrimming.py --fasta "${unmapped_fasta_file}" --output "${unmapped_fasta_file.baseName}_tail.trimmed.fasta"
    """
}

/*
 * define the `REALIGN` process takes tail.trimmed.fasta files
 * and outputs a SAM file.
 */
 
process REALIGN {

    tag { pair_id }

    publishDir "${params.outdir}/${pair_id}/", mode: 'copy'

    input:
    path index_files, stageAs: "bt2_index/*"
    tuple val(pair_id), path(fasta_file_trimmed_tails)

    output:
    path "${fasta_file_trimmed_tails.baseName}_polyAtranscripts.sam"

    script:
    def idx = index_files[0].getBaseName(2)    

     """
    bowtie2 -x "${idx}" -f "${fasta_file_trimmed_tails}" -S "${fasta_file_trimmed_tails.baseName}_polyAtranscripts.sam"
    """
}

workflow {

    index_ch = INDEX(params.ref_genome)

    gtf2bed_ch=GTF2BED(params.gtf_file)

    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists:true)

    trim_out = TRIM(read_pairs_ch)

    sam_for_strand_ch = SAM_FOR_STRAND(index_ch, TRIM.out.trimmed_fastqs)

    check_strand_ch = CHECK_STRANDEDNESS (sam_for_strand_ch, gtf2bed_ch, TRIM.out.trimmed_fastqs)

    find_tail_ch = FINDTAIL (check_strand_ch, TRIM.out.trimmed_fastqs)

    align_input_ch = find_tail_ch.transpose()

    align_ch = ALIGN( index_ch, align_input_ch )

    trim_tails_ch = TRIM_TAILS(align_ch)

    realign_ch = REALIGN (index_ch, trim_tails_ch)
}
