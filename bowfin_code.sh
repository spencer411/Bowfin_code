# Assign names to each forward and reverse sequence reads file
fwd1="bowfin_1.fq.gz"
rev1="bowfin_2.fq.gz"
name1="bowfin"

# Trim PE reads
fastp \
   -i ${fwd1} \
   -I ${rev1} \
   -o ${name1}_F.trimmed.fq.gz \
   -O ${name1}_R.trimmed.fq.gz \
   --detect_adapter_for_pe \
   --cut_front \
   --cut_tail \
   --cut_window_size=4 \
   --cut_mean_quality=20 \
   --qualified_quality_phred=20 \
   --unqualified_percent_limit=30 \
   --n_base_limit=5 \
   --length_required=50 \
   --low_complexity_filter \
   --complexity_threshold=30 \
   --overrepresentation_analysis \
   --json=${name1}.json \
   --html=${name1}.html \
   --report_title="$name1" \
   --thread=8

# Create MP reads using in silico method
cross-mates AHAT01.fasta bowfin_F.trimmed.fq bowfin_R.trimmed.fq

#Overlap paired-end reads with flash
flash bowfin_F.trimmed.fq  bowfin_R.trimmed.fq -M 180

# Remove identical read pairs with clumpify
#PE reads
echo "Running Clumpify PE mode..."
   clumpify.sh \
      -in=bowfin_F.trimmed.fq \
      -in2=bowfin_R.trimmed.fq \
      -out=bowfin_F.trimmed.deduped.fq.gz \
      -out2=bowfin_R.trimmed.deduped.fq.gz \
      dedupe=t \
      containment=f \
      optical=t \
      dupedist=12000
   echo "Finished Clumpify PE mode..."

#SE reads
echo "Running Clumpify SE mode..."
   clumpify.sh \
      -in=bowfin.extendedFrags.fastq \
      -out=bowfin.extendedFrags.deduped.fq.gz \
      dedupe=t \
      containment=f \
      optical=t \
      dupedist=12000
   echo "Finished Clumpify SE mode..."

#Remove mitochondrial reads with bowtie2
#PE reads
name="bowfin"
bowtie2-build bowfin_mito.fasta bowfin_mito
# Do Mapping for bowfin reads
bowtie2 \
   --phred33 \
   -q \
   --very-sensitive \
   --minins 0 \
   --maxins 1000 \
   --fr \
   --threads 20 \
   --reorder \
   -x bowfin_mito \
   -1 ${name}_F.trimmed.deduped.fq.gz \
   -2 ${name}_R.trimmed.deduped.fq.gz | \
   samtools view -b -F 2 | \
   samtools sort -T ${name}.tmp -n -O bam | \
   bedtools bamtofastq -i - -fq ${name}_F.trimmed.deduped.noMito.fq -fq2 ${name}_R.trimmed.deduped.noMito.fq

#SE reads
name="bowfin"
bowtie2-build bowfin_mito.fasta bowfin_mito
# Do Mapping for bowfin reads
bowtie2 \
   --phred33 \
   -q \
   --very-sensitive \
   --minins 0 \
   --maxins 1000 \
   --threads 20 \
   --reorder \
   -x bowfin_mito \
   -U ${name}.extendedFrags.deduped.fq.gz | \
   samtools view -b -F 2 | \
   samtools sort -T ${name}.tmp -n -O bam | \
   bedtools bamtofastq -i - -fq ${name}.extendedFrags.deduped.noMito.fq

#Error correct with musket
#PE reads
musket -p 64 bowfin_F.trimmed.deduped.noMito.fq bowfin_R.trimmed.deduped.noMito.fq -omulti trout.trimmed.deduped.noMito.corrected -inorder

#SE reads
musket -p 64 bowfin.extendedFrags.deduped.noMito.fq -omulti bowfin.extendedFrags.deduped.noMito.corrected -inorder

#Assemble with abyss
abyss-pe k=104 name=bowfin lib='pea peb' mp='mpc mpd mpe mpf mpg mph mpi' \
	pea='bowfin.trimmed.deduped.noMito.corrected_1.fq bowfin.trimmed.deduped.noMito.corrected_2.fq' \
	peb='bowfin_F.trimmed-pe-500_1.fq bowfin_F.trimmed-pe-500_2.fq' \
	se='bowfin.extendedFrags.deduped.noMito.corrected.fq'\
	mpc='bowfin_F.trimmed-mp-1000_1.fq bowfin_F.trimmed-mp-1000_2.fq' \
	mpd='bowfin_F.trimmed-mp-1500_1.fq bowfin_F.trimmed-mp-1500_2.fq' \
	mpe='bowfin_F.trimmed-mp-2000_1.fq bowfin_F.trimmed-mp-2000_2.fq' \
	mpf='bowfin_F.trimmed-mp-5000_1.fq bowfin_F.trimmed-mp-5000_2.fq' \
	mpg='bowfin_F.trimmed-mp-10000_1.fq bowfin_F.trimmed-mp-10000_2.fq' \
	mph='bowfin_F.trimmed-mp-20000_1.fq bowfin_F.trimmed-mp-20000_2.fq' \
	mpi='bowfin_F.trimmed-mp-50000_1.fq bowfin_F.trimmed-mp-50000_2.fq'

#Mapping with bwa, fastp, samtools and seqtk
# Set the input file:
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p samples)
fq=$(grep "$sample" seq_files.list)

# Run SE mapping
./Mapper_SE \
   -r bowfin-8.fa \
   -i $fq \
   -t 4 \
   -o $sample

#Mapper_SE executable script
# V3.2 4/02/2020 by Bob Fitak - SE reads
# This script will take an input fastq file of raw
# reads and clean, map, and process them. The
# reference genome will be indexed using BWA if
# no index is found.

# You may need to adjust other parameters for
# fastp, bwa and samtools as you see fit.
# The program assumes paired-end sequencing
# libraries

# REQUIREMENTS (must be in your PATH):
# bwa
# fastp
# samtools (v1.3 or greater)
# seqtk (Heng Li's github page)

scriptName="${0##*/}"

function printUsage() {
   cat<<EOF

Synopsis

   $scriptName [-h | --help] [-i forward.fastq{.gz}] [-j reverse.fastq{.gz}] [-r reference.fa] [-X integer] [-t integer] [-o output_prefix] [-u]

Description

   Clean, trim, and map a sequencing dataset to
   a reference genome, keeping the mapped or
   unmapped reads file (SAM format).

   -i
      The forward reads file in fastq format (gzipped accepted)

   -j
      The reverse reads file in fastq format (gzipped accepted)

   -r
      the reference genome fasta file for mapping

   -X
      include only X number of reads (for testing purposes)
      Default: include all reads

   -t
      the number of threads (default = 1)

   -o
      prefix for output files (default = "I-should-have-named-these-files")

   -u
      if set (T), then only unmapped reads are kept and the resulting BAM files deleted to save space
      Default: F, don't keep unmapped reads in separate files

   -h, --help
      Brings up this help page

Author
   Robert Fitak, Department of Biology; Genomics and Bioinformatics Cluster, University of Central Florida
   robert.fitak@ucf.edu

EOF
}

if [ $# -lt 1 ]
   then printUsage
   exit 1
fi

# Set default number of threads
threads=1

# Set default for not keeping unmapped reads
unmapped=F

# Set default output files name
out="I-should-have-named-these-files"

# Read input parameters
while getopts ':hi:j:r:X:t:o:u' option
   do
   case "${option}"
   in
   i) fwd=${OPTARG};;
   j) rev=${OPTARG};;
   r) ref=${OPTARG};;
   X) reads=${OPTARG};;
   t) threads=${OPTARG};;
   o) out=${OPTARG};;
   u) unmapped=T;;
   h) printUsage
     exit
     ;;
   help) printUsage
     exit
     ;;
   esac
done

# Check for reference genome index, or index if none is found
if [[ -f $ref.bwt ]]
then
   echo "Reference index exists, skipping reference index build"
else
   echo "Making bwa index"
   bwa index $ref
fi

# Process all sequences
if [[ -z $reads ]]
then
fastp \
   --in1 $fwd \
   --stdout \
   --detect_adapter_for_pe \
   --adapter_fasta /network/rit/home/bobfitak/.bin/adapters.fa \
   --cut_front \
   --cut_tail \
   --cut_window_size=4 \
   --cut_mean_quality=20 \
   --qualified_quality_phred=20 \
   --unqualified_percent_limit=30 \
   --n_base_limit=5 \
   --length_required=50 \
   --low_complexity_filter \
   --complexity_threshold=30 \
   --overrepresentation_analysis \
   --html=${out}.html \
   --json=${out}.json \
   --report_title="$out" \
   --thread=${threads} | \
bwa mem \
   -M \
   -p \
   -t ${threads} \
   ${ref} \
   - | \
samtools sort \
   -T ${out}.tmp \
   -O bam \
   -o ${out}.sorted.bam -

#${bowtie2} \
#   -x ${ref} \
#   --interleaved - \
#   -q \
#   --very-sensitive-local \
#   --no-unal \
#   --threads ${threads} \
#   -S ${SRR}.sam \
#   --un-conc-gz ${SRR}.unmapped.fq.gz

# Or process only X reads
elif [[ ! -z $reads ]]
then
# multiply 4 by the number of reads desired
n=$(( $reads * 4 ))
seqtk seq -l0 $fwd | head -n ${n} | \
fastp \
   --stdin \
   --stdout \
   --detect_adapter_for_pe \
   --adapter_fasta /network/rit/home/bobfitak/.bin/adapters.fa \
   --cut_front \
   --cut_tail \
   --cut_window_size=4 \
   --cut_mean_quality=20 \
   --qualified_quality_phred=20 \
   --unqualified_percent_limit=30 \
   --n_base_limit=5 \
   --length_required=50 \
   --low_complexity_filter \
   --complexity_threshold=30 \
   --overrepresentation_analysis \
   --html=${out}.html \
   --json=${out}.json \
   --report_title="$out" \
   --thread=${threads} | \
bwa mem \
   -M \
   -p \
   -t ${threads} \
   ${ref} \
   - | \
samtools sort \
   -T ${out}.tmp \
   -O bam \
   -o ${out}.sorted.bam -
fi

# Delete extra files
rm -rf ${out}.tmp* ${out}.json

# BAM file, get stats, index.
samtools index ${out}.sorted.bam
samtools stats ${out}.sorted.bam > ${out}.bamstats
rm -rf ${out}.tmp*

# If -u is set, get unmapped reads
if [ "$unmapped" = "T" ]
   then
   echo "Keeping unmapped reads and deleting BAM files"
   samtools fastq \
      -1 ${out}.unmapped.F.fq.gz \
      -2 ${out}.unmapped.R.fq.gz \
      -0 ${out}.unmapped.0.fq.gz \
      -s ${out}.unmapped.s.fq.gz \
      -c 6 \
      -n \
      -f 4 \
      ${out}.sorted.bam
   rm ${out}.sorted.bam ${out}.sorted.bam.bai
fi

#echo "input = $SRR, X = $reads, reference = $ref, threads = $threads"
echo "Mapping analysis complete.  All results are saved with the prefix: $out"

#SNP calling with STACKS

/network/rit/home/sb939359/turnerlab/bin/stacks-2.53/bin/ref_map.pl -T 16 \
--popmap /network/rit/home/sb939359/turnerlab/spencer/stacks/names.txt \
-X "populations:--vcf" \
-X "populations:--structure" \
-X "populations:--phylip" \
-o /network/rit/home/sb939359/turnerlab/spencer/stacks/output \
--samples /network/rit/home/sb939359/turnerlab/spencer/stacks \

vcftools --vcf filteredSNPs.vcf --remove-filtered-all --maf 0.10 --remove-indels --max-missing 0.5 --recode --recode-INFO-all --out SNPs_filtered_maf0.10_no_indels_50percent
