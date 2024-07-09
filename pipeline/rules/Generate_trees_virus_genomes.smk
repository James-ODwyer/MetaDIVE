

# script picks up after reference generated genome and spades assembled genome are made
# the results this picks up from is the multifastas per species per sample.
# Not sure whether it is best to run the mafft as parallel or just as a for loop. external cross sample paralelisation will occur through snakemake but there may be 5+ genomes per sample which will be done in series
# But if I run parallel it will require I devote a large guaranteed number of cores and total memory which may cause wastage whem there is only 1-2 viruses present. I have set it up as linear for loop for the time being

# run with -ge 3 because there are always at least 3 folders if a genome is assembled. hits, references, finalgenomes

rule align_viral_seqs:
    message:
        """
        Advanced viral tree building step 1.
        Creating Mafft alignments for each virus identified for {wildcards.sample}
        """
    input:
        finished = config["sub_dirs"]["Viral_genomes_present_refined_genomes"] + "/{sample}/finished2.txt"
    output:
        finished2 = config["sub_dirs"]["Viral_genomes_present_mafft"] + "/{sample}/finished.txt"
    params:
        programdir = config["program_dir"],
        samplename = "{sample}",
        inresultsdir = config["sub_dirs"]["Viral_genomes_present_refined_genomes"] + "/{sample}",
        ingenomedir = config["sub_dirs"]["Viral_genomes_present_refined_genomes"] + "/{sample}/finished_genomes",
        outconsensusdir = config["sub_dirs"]["Viral_genomes_present_mafft"] + "/{sample}"
    log:
        "logs/" + config["sub_dirs"]["Viral_genomes_present_mafft"] + "/{sample}.txt"
    threads: 4
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Viral_genomes_present_mafft"] + "/{sample}.txt"
    conda: "snakemake7"
    shell:
        """
        if [ ! -d {params.outconsensusdir} ]; then mkdir -p {params.outconsensusdir}; fi && \
        length=$(find {params.inresultsdir} -type d | wc -l) && \
        if [ ${{length}} -ge 3 ]
        then
        multifastas=(`ls {params.ingenomedir}/*_multifasta/*`)
        for file in ${{multifastas[@]}}; do
        speciesname=(`awk '{{sub(/_combined.*/, ""); sub(/.*_multifasta\//, ""); print}}' <<<$file`)
        echo "$speciesname"
        mafft --auto --reorder --maxiterate 50 --retree 20 --adjustdirection --thread {threads} $file > {params.outconsensusdir}/"$speciesname"_mafft.aln
        done
        fi && \
        if [ ${{length}} -le 2 ]
        then
        echo " No genomes were returned for sample {params.samplename} "
        fi && \
        touch {output.finished2}
        """

# awk sub may call just the first \/ but I am hoping it does all. double check if fail
# Trimal has no parallelisation so could be a candidate to run parallel GNU through.

rule trimal_alignments:
    message:
        """
        Advanced viral tree building step 2.
        Creating trimming alignments and sequences for {wildcards.sample} using Trimal
        """
    input:
        finished = config["sub_dirs"]["Viral_genomes_present_mafft"] + "/{sample}/finished.txt"
    output:
        finished2 = config["sub_dirs"]["Viral_genomes_present_trimal"] + "/{sample}/finished.txt"
    params:
        programdir = config["program_dir"],
        samplename = "{sample}",
        inmafftdir = config["sub_dirs"]["Viral_genomes_present_mafft"] + "/{sample}",
        outconsensusdir = config["sub_dirs"]["Viral_genomes_present_trimal"] + "/{sample}"
    priority: 10
    log:
        "logs/" + config["sub_dirs"]["Viral_genomes_present_trimal"] + "/{sample}.txt"
    threads: 2
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Viral_genomes_present_trimal"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.outconsensusdir} ]; then mkdir -p {params.outconsensusdir}; fi && \
        length=$(find {params.inmafftdir} -type f | wc -l) && \
        if [ ${{length}} -ge 2 ]
        then
        maffts=(`ls {params.inmafftdir}/*mafft*`)
        for file in ${{maffts[@]}}; do
        speciesname=(`awk '{{sub(/_mafft.*/, ""); sub(/.*\//, ""); print}}' <<<$file`)
        trimal -in $file -out {params.outconsensusdir}/"$speciesname"_trimal.aln -htmlout {params.outconsensusdir}/"$speciesname"_results_trim.html -automated1 -keepheader 2> {log}
        done
        fi && \
        if [ ${{length}} -le 1 ]
        then
        echo " No genomes were returned for sample {params.samplename} "
        touch {output.finished2} 
        fi && \
        touch {output.finished2}
        """

rule iq_tree:
    message:
        """
        Advanced viral tree building step 3.
        generating ML tree for {wildcards.sample} using IQ-TREE
        """
    input:
        finished = config["sub_dirs"]["Viral_genomes_present_trimal"] + "/{sample}/finished.txt"
    output:
        finished2 = config["sub_dirs"]["Viral_genomes_present_iqtree"] + "/{sample}/finished.txt"
    params:
        programdir = config["program_dir"],
        samplename = "{sample}",
        inmtrimaldir = config["sub_dirs"]["Viral_genomes_present_trimal"] + "/{sample}",
        outconsensusdir = config["sub_dirs"]["Viral_genomes_present_iqtree"] + "/{sample}"
    priority: 10
    log:
        "logs/" + config["sub_dirs"]["Viral_genomes_present_iqtree"] + "/{sample}.txt"
    threads: 6
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Viral_genomes_present_iqtree"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.outconsensusdir} ]; then mkdir -p {params.outconsensusdir}; fi && \
        length=$(find {params.inmtrimaldir} -type f | wc -l) && \
        if [ ${{length}} -ge 2 ]
        then
        maffts=(`ls {params.inmtrimaldir}/*.aln`)
        for file in ${{maffts[@]}}; do
        speciesname=(`awk '{{sub(/_trimal.aln/, ""); sub(/.*{params.samplename}\//, ""); print}}' <<<$file`)
        filename=(`awk '{{sub(/.*{params.samplename}\//, ""); print}}' <<<$file`)
        echo "Virus species to generate tree $speciesname "
        if [ ! -d {params.outconsensusdir}/"$speciesname" ]; then mkdir {params.outconsensusdir}/"$speciesname"; fi
        cp $file {params.outconsensusdir}/"$speciesname"
        samppresence=(`grep "{params.samplename}" {params.outconsensusdir}/"$speciesname"/"$filename" | wc -l`)
        if [ ${{samppresence}} -ge 1 ]; then
        iqtree -s {params.outconsensusdir}/"$speciesname"/"$filename" -B 1000 -T AUTO --threads-max {threads} 2> {log}
        fi
        if [ ${{samppresence}} -eq 0 ]; then
        echo " All viral sequences derrived from the sample were removed due to poor quality "
        touch {params.outconsensusdir}/"$speciesname"/tree_not_generated.txt
        fi
        done
        fi && \
        if [ ${{length}} -le 1 ]
        then
        echo " No genomes were returned for sample {params.samplename} "
        touch {output.finished2} 
        fi && \
        touch {output.finished2}
        """
