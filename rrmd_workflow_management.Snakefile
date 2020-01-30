################################################################################
# Project: Archaic genome comparison 
# Part: Quality checks
# Step: Average coverage on chromosome 21
#
# Calculate the average coverage on the chromosome 21 for the three Archaic
# samples.
#
# Alex Huebner, 29/01/2020
################################################################################

workdir: "./coverage"

# Infer sample names from file system
SAMPLES, = glob_wildcards("../data/{sample}.hg19.bam")
print(SAMPLES)

rule all:
    input: 
        expand("{sample}.21_avgcov.txt", sample=SAMPLES[0])


rule sort_by_coordinate:
    output:
        "{sample}.sorted.bam"
    message: "Sort sample {wildcards.sample} by coordinate"
    params: 
        bam = "../data/{sample}.hg19.bam"
    shell:
        """
        samtools sort -o {output} {params.bam}
        """

rule index_sorted_bam:
    input:
        "{sample}.sorted.bam"
    output:
        "{sample}.sorted.bam.bai"
    message: "Create BAM index for sample {wildcards.sample}"
    shell:
        """
        samtools index {input}
        """

rule subset_chrom:
    input:
        bam = "{sample}.sorted.bam",
        bai = "{sample}.sorted.bam.bai"
    output:
        "{sample}.21.bam"
    message: "Subset sample {wildcards.sample} to chromosome 21"
    shell:
        """
        samtools view -bh {input.bam} > {output}
        """

rule depth:
    input:
        "{sample}.21.bam"
    output:
        "{sample}.21_depth.txt.gz"
    message: "Determine the depth on chromosome 21 for sample {wildcards.sample}"
    shell:
        """
        samtools depth -a {input} | bgzip > {output}
        """

rule average_depth:
    input:
        "{sample}.21_depth.txt.gz"
    output:
        "{sample}.21_avgcov.txt"
    message: "Calculate the average coverage for sample {wildcards.sample}"
    shell:
        """
        zcat {input} | awk '{{totalcov += $3}}END{{print totalcov / NR}}' - > {output}
        """
