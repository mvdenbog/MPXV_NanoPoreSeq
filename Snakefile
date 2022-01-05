
WORKFLOW = 'workflow_mapping_to_assembly'

REF, = glob_wildcards("/data/{ref}.fasta")

rule all:
     input:
         canu_fastmer=expand("/data/{wf}_canu_d/{wf}_canu_fastmer.txt", wf=WORKFLOW),
         canu_homopolish=expand("/data/{wf}_canu_d/homopolish_output/{wf}_canu_homopolished_fastmer.txt", wf=WORKFLOW),
         raconfasta=expand("/data/{wf}_Racon/reads.racon2.fasta", wf=WORKFLOW),
         racon_medaka_homopolish=expand("/data/{wf}_RaconMedakaHomopolish/homopolish_output/consensus_canu_racon_medaka_homopolished_fastmer.txt", wf=WORKFLOW),
         canu_medaka_homopolish=expand("/data/{wf}_canu_d/Medaka/homopolish_output/consensus_canu_medaka_homopolished_fastmer.txt", wf=WORKFLOW),
         dir="/data/FastMer",
         canu_medaka=expand("/data/{wf}_canu_d/Medaka/consensus_canu_medaka_fastmer.txt", wf=WORKFLOW),
         canu_racon_fastmer=expand("/data/{wf}_Racon/reads.racon2.fasta.txt", wf=WORKFLOW)
     shell:
         '''
           cp {input.canu_fastmer} {input.dir}/canu.txt
           cp {input.canu_homopolish} {input.dir}/canu_homopolish.txt
           cp {input.canu_medaka} {input.dir}/canu_medaka.txt
           cp {input.canu_medaka_homopolish} {input.dir}/canu_medaka_homopolish.txt
           cp {input.canu_racon_fastmer} {input.dir}/canu_racon.txt
           cp {input.racon_medaka_homopolish} {input.dir}/racon_medaka_homopolish.txt
           ( echo -n -e "Methods\t"; cut -f 2-8  {input.dir}/canu.txt  | head -n 1 ) > {input.dir}/result.TXT
           for i in {input.dir}/*.txt; do base=`basename \$i .txt`; echo -n ${{base}}; echo -n -e "\t";cut -f 2-8 $i | tail -n +2; done >> {input.dir}/result.TXT
           echo "Pipeline \'MPXV_NanoPoreSeq\' complete."
         '''

rule make_fastmer_dir:
     output:
        fastmerdir=directory("/data/FastMer")
     shell:
         'if [[ ! -d {output.fastmerdir} ]]; then mkdir {output.fastmerdir}; fi'

rule merge_fastq:
     output:
         expand("/data/{wf}_input.fastq", wf=WORKFLOW)
     params:
        dir=expand("/data/{wf}_summary-plots-log-transformed", wf=WORKFLOW)
     shell:
         '''
         cat /data/*.fastq > {output}
         # NanoPlot -t 2 --fastq {output}  --loglength -o {params.dir} --plots hex dot pauvre kde
         '''

rule map:
     input:
         fq=expand("/data/{wf}_input.fastq", wf=WORKFLOW),
         fasta=expand("/data/{ref}.fasta", ref=REF)
     output:
         expand("/data/{wf}_mnmp_bam_mapped.fastq", wf=WORKFLOW)
     params:
        canu=expand("{wf}_canu", wf=WORKFLOW),
        canu_d=expand("/data/{wf}_canu_d", wf=WORKFLOW),
        mapping=expand("/data/{wf}_mnmp", wf=WORKFLOW)
     shell:
        '''
         echo {input.fasta}
         echo minimap ...
         minimap2 -a -x map-ont {input.fasta} {input.fq}  > {params.mapping}
         echo samtools view ... 
         samtools view -bt {input.fasta} -o {params.mapping}.unsorted.bam {params.mapping}
         echo samtools sort ...
         samtools sort -o {params.mapping}.sorted.bam {params.mapping}.unsorted.bam 
         echo samtools depth ...
         samtools depth {params.mapping}.sorted.bam > {params.mapping}.sorted.bam.depth
         echo image ...
         /opt/Conda/bin/python /pipeline/scripts/depth_of_coverage_v2.py {params.mapping}.sorted.bam.depth
         echo mapped reads bam ...
         samtools view -b -F 4 {params.mapping}.sorted.bam > {params.mapping}.mapped.bam
         echo mapped reads fastq ...
         samtools bam2fq {params.mapping}.mapped.bam > {output}
        '''

rule assemble:
     input: expand("/data/{wf}_mnmp_bam_mapped.fastq", wf=WORKFLOW)
     output: 
        assembly=expand("/data/{wf}_canu_d/{wf}_canu.contigs.fasta", wf=WORKFLOW),
        canu_trimmed_reads=expand("/data/{wf}_canu_d/{wf}_canu.trimmedReads.fasta.gz", wf=WORKFLOW)
     params:
        canu=expand("{wf}_canu", wf=WORKFLOW),
        canu_d=expand("/data/{wf}_canu_d", wf=WORKFLOW)
     shell:
        '''
         echo canu ...
         canu -p {params.canu} -d {params.canu_d}  maxThreads=1 useGrid=false maxThreads=4 stopOnReadQuality=false genomeSize=200k minReadLength=500 -nanopore-raw {input}
        '''

rule canu_fastmer:
    input: assembly=expand("/data/{wf}_canu_d/{wf}_canu.contigs.fasta", wf=WORKFLOW),fasta=expand("/data/{ref}.fasta", ref=REF)
    output: fastmer=expand("/data/{wf}_canu_d/{wf}_canu_fastmer.txt", wf=WORKFLOW)
    shell:
       '''
        python2 /opt/FastMer/fastmer.py --reference {input.fasta} --assembly {input.assembly} > {output.fastmer}
       '''


rule canu_homopolish:
    input: expand("/data/{wf}_canu_d/{wf}_canu.contigs.fasta", wf=WORKFLOW)
    output:
        expand("/data/{wf}_canu_d/homopolish_output/{wf}_canu_homopolished.fasta", wf=WORKFLOW),
    params:
       outdir=expand("/data/{wf}_canu_d/homopolish_output", wf=WORKFLOW)
    shell:
       '''
        /opt/Conda/envs/homopolish/bin/python /opt/homopolish/homopolish.py polish -a {input} -s /opt/homopolish/virus.msh -m /opt/homopolish/R9.4.pkl -o {params.outdir}
       '''

rule canu_homopolish_fastmer:
    input:
       assembly=expand("/data/{wf}_canu_d/homopolish_output/{wf}_canu_homopolished.fasta", wf=WORKFLOW),
       fasta=expand("/data/{ref}.fasta", ref=REF),
    output:
       fastmer1=expand("/data/{wf}_canu_d/homopolish_output/{wf}_canu_homopolished_fastmer.txt", wf=WORKFLOW),
    shell:
       '''
        python2 /opt/FastMer/fastmer.py --reference {input.fasta} --assembly {input.assembly} > {output.fastmer1}
       '''

rule canu_medaka:
    input:
       reads=expand("/data/{wf}_canu_d/{wf}_canu.trimmedReads.fasta.gz", wf=WORKFLOW),
       assembly=expand("/data/{wf}_canu_d/{wf}_canu.contigs.fasta", wf=WORKFLOW),
       fasta=expand("/data/{ref}.fasta", ref=REF),
    params:
       dir=expand("/data/{wf}_canu_d/Medaka/", wf=WORKFLOW)
    output:
       medaka=expand("/data/{wf}_canu_d/Medaka/consensus.fasta", wf=WORKFLOW),
       fastmer=expand("/data/{wf}_canu_d/Medaka/consensus_canu_medaka_fastmer.txt", wf=WORKFLOW)
    shell:
       '''
        medaka_consensus -i {input.reads} -d {input.assembly} -o {params.dir}
        python2 /opt/FastMer/fastmer.py --reference {input.fasta} --assembly {output.medaka} > {output.fastmer} 
       '''

rule canu_medaka_homopolish:
    input:
       medakaconsensus=expand("/data/{wf}_canu_d/Medaka/consensus.fasta", wf=WORKFLOW),
       fasta=expand("/data/{ref}.fasta", ref=REF),
    params:
       homopolish_dir=expand("/data/{wf}_canu_d/Medaka/homopolish_output", wf=WORKFLOW)
    output:
       fastmer=expand("/data/{wf}_canu_d/Medaka/homopolish_output/consensus_canu_medaka_homopolished_fastmer.txt", wf=WORKFLOW)
    shell:
       '''
       /opt/Conda/envs/homopolish/bin/python /opt/homopolish/homopolish.py polish -a {input.medakaconsensus} -s /opt/homopolish/virus.msh -m /opt/homopolish/R9.4.pkl -o {params.homopolish_dir}
       python2 /opt/FastMer/fastmer.py --reference {input.fasta} --assembly   {params.homopolish_dir}/consensus_homopolished.fasta  > {output.fastmer}
       '''

rule canu_racon:
    input:
       nmnp_mapped=expand("/data/{wf}_mnmp_bam_mapped.fastq", wf=WORKFLOW),
       fasta=expand("/data/{ref}.fasta", ref=REF),
    params:
       racon_dir=expand("/data/{wf}_Racon", wf=WORKFLOW)
    output:
       expand("/data/{wf}_Racon/reads.racon2.fasta", wf=WORKFLOW),
       expand("/data/{wf}_Racon/reads.racon2.fasta.txt", wf=WORKFLOW)
    shell:
      '''
        minimap2 -x ava-ont -t 4 {input.nmnp_mapped} {input.nmnp_mapped} | gzip -1 > {params.racon_dir}/reads.paf.gz
        miniasm -f {input.nmnp_mapped} {params.racon_dir}/reads.paf.gz > {params.racon_dir}/reads.gfa
        awk '$1 ~/S/ {{print ">"$2"\\n"$3}}' {params.racon_dir}/reads.gfa > {params.racon_dir}/reads.fasta
        minimap2 -t 4 {params.racon_dir}/reads.fasta {input.nmnp_mapped} > {params.racon_dir}/reads.gfa1.paf
        racon -t 4 {input.nmnp_mapped} {params.racon_dir}/reads.gfa1.paf {params.racon_dir}/reads.fasta > {params.racon_dir}/reads.racon1.fasta
        minimap2 -t 4 {params.racon_dir}/reads.racon1.fasta {input.nmnp_mapped} > {params.racon_dir}/reads.gfa2.paf
        racon -t 4 {input.nmnp_mapped} {params.racon_dir}/reads.gfa2.paf {params.racon_dir}/reads.racon1.fasta > {params.racon_dir}/reads.racon2.fasta
        python2 /opt/FastMer/fastmer.py --reference {input.fasta} --assembly {params.racon_dir}/reads.racon2.fasta > {params.racon_dir}/reads.racon2.fasta.txt
      '''

rule canu_racon_medaka_homopolish:
   input:
     canu_trimmed_reads=expand("/data/{wf}_canu_d/{wf}_canu.trimmedReads.fasta.gz", wf=WORKFLOW),
     racon_assembly=expand("/data/{wf}_Racon/reads.racon2.fasta", wf=WORKFLOW),
     fasta=expand("/data/{ref}.fasta", ref=REF),
   params:
     outdir=expand("/data/{wf}_RaconMedakaHomopolish", wf=WORKFLOW)
   output:
     expand("/data/{wf}_RaconMedakaHomopolish/homopolish_output/consensus_canu_racon_medaka_homopolished_fastmer.txt", wf=WORKFLOW)
   shell:
     '''
       medaka_consensus -i {input.canu_trimmed_reads} -d {input.racon_assembly} -o {params.outdir}
       /opt/Conda/envs/homopolish/bin/python /opt/homopolish/homopolish.py polish -a {params.outdir}/consensus.fasta  -s /opt/homopolish/virus.msh -m /opt/homopolish/R9.4.pkl -o {params.outdir}/homopolish_output
       python2 /opt/FastMer/fastmer.py --reference {input.fasta}  --assembly {params.outdir}/homopolish_output/consensus_homopolished.fasta > {output}
     '''
