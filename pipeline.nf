/*                         
    *** Create a Nextflow pipline that helps to evaluate and trim the bacterial genome sequencing 
    data (paired-end, FASTQ) using fastp, then performs de novo assembly of this data using SPAdes, 
    and finally annotates the assembly results using Prokka.

    This pipeline was run on WGS data of Mycobacterium ulcerans from EMBL-EBI: ERR3335404. Paired-end
    read data was scaled down (take first 4000 line each file) for fast performance and save the 
    computer resources. 
    Link to EMBL database: https://www.ebi.ac.uk/ena/browser/view/ERR3335404 

    *** Raw read data was stored in fastq/. The output of this pipeline will be stored in TEST_OUTPUT/ 
    once executed.
*/

// Create a parameter that store the path of TEST_OUTPUT/ 
params.output = '/home/nhantruong/small_WGS/TEST_OUTPUT'


// 1. Quality control report and data trimming with fastp
process fastp_QC {
    // Copy all files, directories of this process to /home/nhantruong/small_WGS/TEST_OUTPUT/FASTP_OUTPUT
    publishDir "${params.output}/FASTP_OUTPUT", mode: 'copy'

    input:
        // Declare the nescessary input for this process
        path Raw_read1 // Forward strand raw read data
        path Raw_read2 // Reverse strand raw read data
    
    output:
        // Declare the output of this process
        // Use emit to name the output channel so we can get the output files of this 
        // process later for next process
        path "${name1}_trimmed.fq.gz", emit: trimmed_1 // E.g this channel will have the trimmed file of read1  
        path "${name2}_trimmed.fq.gz", emit: trimmed_2
        path "*"

    script:
    // Create name1, name2 and assign it value to get specific name for trimmed output
    // For example: Raw_read1 = '/home/nhantruong/small_WGS/fastq/ERR3335404_1-small.fastq'
    // .getBaseName() method extract the name 'RR3335404_1-small' without dir and extension
    // .split('-')[0] will cut that name into an array devided by '-' then get first element 'ERR3335404_1'
    name1 = Raw_read1.getBaseName().split('-')[0]
    name2 = Raw_read2.getBaseName().split('-')[0]
    """
    fastp -i ${Raw_read1} -I ${Raw_read2} -o '${name1}_trimmed.fq.gz' -O '${name2}_trimmed.fq.gz' 
    """
    // The script above run fastp with 2 input: -i <Forward-strand read data> -I <Reverse-strand read data>
    // and the outputs are: -o <Forwrd-strand trimmed data> -O <Reverse-strand trimmed data> 
    // The quality control report will be created as default: fastp.html, fastp.json
}


// 2. Genome assembly with SPAdes
process spades_assembly {
    // Copy all files, directories of this process to /home/nhantruong/small_WGS/TEST_OUTPUT/SPADES_OUTPUT
    publishDir "${params.output}/SPADES_OUTPUT", mode: 'copy'

    input:
        // declare input of this process: required two paths
        path trimmed_1
        path trimmed_2
    output:
        // get contigs.fasta into the assembled output channel
        path "${spades_out}_assembled/contigs.fasta", emit: assemnled 

    script:
    // get the name ERR3335404
    spades_out = trimmed_1.getBaseName().split('_')[0]

    // run SPAdes --careful: tries to reduce number of mismatches and short indels
    // -1 <trimmed data of Forward strand>; -2 <trimmed data of Reverse strand>
    // -o <name of desired ouput directory>
    """
    spades.py --careful -1 $trimmed_1 -2 $trimmed_2 -o '${spades_out}_assembled'
    """
}


// 3. Genome annotation with Prokka
process prokka_annotate {
    // Copy all files, directories of this process to /home/nhantruong/small_WGS/TEST_OUTPUT/PROKKA_OUTPUT
    publishDir "${params.output}/PROKKA_OUTPUT", mode: 'copy'

    input: 
        path assembled_data
    output:
        // take all output files of this process 
        path "*" 

    script:
    // run prokka --cpus 1 => tell prokka to use just 1 processor of this device
    // --prefix result => set up the name of output directory and files 
    // --kingdom Bacteria => set annotation mode for Bacteria
    // ${assembled_data} => file path of assembled data
    """
    prokka --cpus 1 --prefix result --kingdom Bacteria ${assembled_data}
    """
}


// define parameters that store the path of raw read data
params.read1 = '/home/nhantruong/small_WGS/fastq/ERR3335404_1-small.fastq'
params.read2 = '/home/nhantruong/small_WGS/fastq/ERR3335404_2-small.fastq'

workflow {
    // call process fastp_QC for trimming QC report
    fastp_QC(params.read1, params.read2)

    // call process spades_assembly, this process take up the output file of fastp_QC as input
    spades_assembly(fastp_QC.out.trimmed_1, fastp_QC.out.trimmed_2)

    // call process prokka_annotate, this process take up the output file of spades_assembly as input
    prokka_annotate(spades_assembly.out.assemnled)
}