# Only need to define wildcards agin if you want an if/else statement in defining them (May be useful for combined analysis. # will be useful for plans about host removal

def get_raw_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule unmapped_contigs_nucleotide:
    message:
        """
        DNA classification step 1 (starts after analysis step 3) 
        Blast unassigned contigs {wildcards.sample} blastn
        """
    input:
        contigsunassigned = config["sub_dirs"]["contigs_assigned"] + "/{sample}_matches_unassigned.fa"
    output:
        fasta_clust = config["sub_dirs"]["contigs_assigned_nucl"] + "/{sample}_reads_clust.fasta",
        fasta_clust_ref_file = config["sub_dirs"]["contigs_assigned_nucl"] + "/{sample}_reads_clust.fasta.clstr",
        blastfile = config["sub_dirs"]["contigs_assigned_nucl"] + "/{sample}_matches_nucleotide.m8"
    params:
        blastdb= config["blast_nucleotide_database"],
        blastnfolder = config["sub_dirs"]["contigs_assigned_nucl"]
    log:
        "logs/" + config["sub_dirs"]["contigs_assigned_nucl"] + "/{sample}.log"
    threads: 6
    resources:
        mem_mb=16000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["contigs_assigned_nucl"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.blastnfolder} ]
        then
            mkdir {params.blastnfolder}
        fi && \
        lengthunassigned=(`wc -l {input.contigsunassigned}`) && \
        if [ ${{lengthunassigned}} -ge 3 ]
        then
        cd-hit -i {input.contigsunassigned} -o {output.fasta_clust} -c 0.95 -n 5 -T {threads} -d 0 -M 15000
        blastn -query {output.fasta_clust} \
            -db {params.blastdb} \
            -evalue 0.00001 \
            -max_target_seqs 5 \
            -max_hsps 1 \
            -outfmt '6 qseqid sseqid pident length evalue bitscore staxids stitle qcovhsp' \
            -num_threads {threads} \
            -word_size 33 \
            -out {output.blastfile} \
            2> {log}
        fi && \
        touch {output.blastfile} && \
        touch {output.fasta_clust} && \
        touch {output.fasta_clust_ref_file}
        """

rule analyse_blastn_hits:
    message:
        """
        DNA classification step 2
        generate taxonomies from each Diamond hit for {wildcards.sample} raw reads using R and taxonomizr package
        """
    input:
        blastfile = config["sub_dirs"]["contigs_assigned_nucl"] + "/{sample}_matches_nucleotide.m8",
        unassigned_contigs = config["sub_dirs"]["contigs_assigned"] + "/{sample}_matches_unassigned.fa",
        full_prot_host_id = config["sub_dirs"]["progress_main_pipeline_host"] + "/{sample}_host_id.txt",
        clusterfile = config["sub_dirs"]["contigs_assigned_nucl"] + "/{sample}_reads_clust.fasta.clstr"
    output:
        sampfinished = config["sub_dirs"]["contigs_assigned_nucl_abundances"] + "/{sample}_finished.txt",
        unassignedDNAcontigs = config["sub_dirs"]["contigs_assigned_nucl"] + "/{sample}_matches_unassigned.fa",
        contigsassignedfile = config["sub_dirs"]["contigs_assigned_nucl"] + "/abundances/{sample}_Contigsallinformationassignment.txt",
        contigshostassigned = config["sub_dirs"]["contigs_assigned_nucl"] + "/abundances/{sample}_host_aligned_contigs_list.txt"
    params:
        samplename = "{sample}",
        namenodedatabase = config["Accession_allnamenode"],
        abundances = config["sub_dirs"]["contigs_assigned_nucl"] + "/abundances/",
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["contigs_assigned_nucl"],
        contigname = "{sample}_matches_unassigned.fa"
    log:
        "logs/" + config["sub_dirs"]["contigs_assigned_nucl_abundances"] + "/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["contigs_assigned_nucl_abundances"] + "/{sample}.txt"
    resources:
        mem_mb=3000
    conda: "Rdataplotting"
    threads: 1
    shell:
        """
        length=(`wc -l {input.blastfile}`) && \
        if [ ${{length}} -ge 1 ]
        then
        Rscript {config[program_dir]}scripts/extract_blastn_best_hits_HPC.R \
            --inputblastn {input.blastfile} --name {params.samplename} \
            --threads {threads} --Accnode {params.namenodedatabase} --output {output.unassignedDNAcontigs} --programdir {params.basedir} --savdir {params.wrkdir} \
            --inputcontig {input.unassigned_contigs} --abundances {params.abundances} --savcontig {params.contigname} --hostid {input.full_prot_host_id} --clustfile {input.clusterfile}
        fi && \
        touch {output.sampfinished} && \
        touch {output.contigsassignedfile} && \
        touch {output.unassignedDNAcontigs} && \
        touch {output.contigshostassigned}
        """
