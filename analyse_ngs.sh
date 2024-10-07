#!/bin/bash

ORIGINAL_DIR=$(pwd)
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
SCRIPT_NAME="$(basename "$0")"
# DIR_SRC=/mnt/chaelab/rachelle/src

HOSTNAME=$(hostname)

PREFIX_DEFAULT="analyseNGS"
GZ_DEFAULT=1
# FASTQ_TO_SAMPLE_DEFAULT="echo '%f_%f' | sed 's/_S0/_S/'"
SAMPLE_ID_FORMAT_DEFAULT=1
SKIP_TO_STR_DEFAULT="start"
STOP_AT_STR_DEFAULT="end"
## derep options
TRUNCATION_LENGTH_DEFAULT=70
FASTQ_SUFFIX_PATTERN_DEFAULT='_L001_R\d_001'
## derep dada2 options
maxEE_DEFAULT=2
RM_PHIX_DEFAULT=TRUE
MULTITHREAD_DEFAULT=TRUE
## blast options
if [[ "${HOSTNAME}" == 'chaelab-ws.nus.edu.sg' ]]; then
    FA_REF_DEFAULT=/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta
    GFF_REF_DEFAULT=/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.gff
elif [[ "${HOSTNAME}" == 'chaelab2' ]]; then
    FA_REF_DEFAULT=/media/HDD3/rachelle/data/genomes/Athaliana/TAIR10/a_thaliana_all.fasta
    GFF_REF_DEFAULT=/media/HDD3/rachelle/data/genomes/Athaliana/TAIR10/TAIR10_GFF3_genes.gff
fi
BLAST_SHORT_THRESHOLD_DEFAULT=50
## read identity options
EXCEL_GENE_PATTERN_DEFAULT='([^_]+)$'
## misc paths
if [[ "${HOSTNAME}" == 'chaelab-ws.nus.edu.sg' ]]; then
    CONDA_DEFAULT=/home/rachelle/miniconda3/condabin/conda
elif [[ "${HOSTNAME}" == 'chaelab2' ]]; then
    CONDA_DEFAULT=/media/HDD3/rachelle/programmes/Conda/miniconda3/bin/conda
fi
## crispresso2 wrapper options
CRISPRESSO_WRAPPER_VERSION_DEFAULT='4-3'
EDITOR_DEFAULT=Cas9

params=${@}

while (( "$#" )); do
    case "$1" in
        ## general parameters
        -d|--dir|--output_dir) DIR=$(realpath "${2}");; ## output directory
        -p|--prefix) PREFIX="${2}";; ## output prefix
        --skip-to|--start-at) SKIP_TO_STR="${2}";; ## valid inputs: derep, blast, intersect, id-seq, id-read, demultiplex
        --stop-at) STOP_AT_STR="${2}";; ## valid inputs: derep, blast, intersect, id-seq, id-read, demultiplex
        --sample-id-format) SAMPLE_ID_FORMAT="${2}";; ## id of sample ID format (default=1); see scripts/fastq2sample.sh for what each id's format is
        --fastq-suffix) FASTQ_SUFFIX_PATTERN="${2}";; ## fastq file suffix regex; use 'R\d' pattern for R1/R2 matching; do not include extensions '.fastq' and/or '.gz'
        ## input parameters
        -f|--fastq|--fastq_dir) DIR_FASTQ=$(realpath "${2}");; ## path to directory containing fastq files
        --gz) GZ=1;; ## raise if fastq files are gzipped (extension fastq.gz)
        --unzipped) GZ=0;; ## raise if fastq files are not gzipped (extension fastq)
        ## derep options
        -t|--truncation) TRUNCATION_LENGTH="${2}";; ## trim reads to first X bp for analysis
        ## derep dada2 parameters
        --maxEE) maxEE="${2}";; ## maxEE parameter for dada2::filterAndTrim
        --rm-phix) RM_PHIX="${2}";; ## rm.phix parameter for dada2::filterAndTrim
        -m|--multithread) MULTITHREAD="${2}";; ## multithread parameter for dada2::learnErrors
        ## blast options
        --assembly) FA_REF="${2}";; ## path to reference assembly FASTA file
        --gff) GFF_REF="${2}";; ## path to reference GFF file
        --blast-short-threshold) BLAST_SHORT_THRESHOLD="${2}";; ## maximum length for blastn-short task
        ## read identity options
        -e|--excel|-b|--booking-sheet|--amplicon_book) EXCEL="${2}";; ## path to Excel sheet of ngs metadata
        -n|-s|--sheetname|--sheet) EXCEL_SHEETNAME="${2}";; ## name of Excel sheet containing metadata
        --metadata|--excel-tsv) EXCEL_TSV="${2}";; ## path to tsv of Excel sheet of ngs metadata
        --excel-gene-pattern) EXCEL_GENE_PATTERN="${2}";; ## gene regex; used with str_extract to get gene ID from Excel file
        ## demultiplexing options (none lol)
        ## crispresso2 options
        --crispresso-wrapper-version) CRISPRESSO_WRAPPER_VERSION="${2}";; ## Crispresso2_wrapper version
        --amplicon|--amplicon_fasta_dir) DIR_AMPLICON=$(realpath "${2}");; ## directory of amplicon FASTA files
        --editor) EDITOR="${2}";; ## valid values: ABE, Cas9, unknown
        ## misc paths
        --conda) CONDA="${2}";; ## path to conda executable
        -h|--help) man -l ${SCRIPT_DIR}/MANUAL_analyse_ngs.1; exit 0;;
        --readme) cat ${SCRIPT_DIR}/README.md; exit 0;;
    esac
    shift
done

## check required arguments
if [ -z "${DIR}" ]; then
    echo "Directory required. Please provide directory using '-d <path to directory>'"
    exit 1
fi
## throw error if Excel inputs not provided in correct combination
if ( [[ ! -z "${EXCEL_TSV}" ]] && ( [[ ! -z "${EXCEL}" ]] || [[ ! -z "${EXCEL_SHEETNAME}" ]] ) ) ||
       ( [[ ! -z "${EXCEL}" ]] && [[ -z "${EXCEL_SHEETNAME}" ]] ) ||
       ( [[ -z "${EXCEL}" ]] && [[ ! -z "${EXCEL_SHEETNAME}" ]] ); then
    echo "Use only '--excel-tsv <path to file>' OR '--excel <path to file> --sheetname <sheet name>'."
    exit 1
fi

## apply defaults
## use sheetname as output prefix if --sheetname is used and --prefix is not used
PREFIX=$(echo "${PREFIX:-${EXCEL_SHEETNAME:-${PREFIX_DEFAULT}}}" | sed 's/ /_/g')
SAMPLE_ID_FORMAT="${SAMPLE_ID_FORMAT:-${SAMPLE_ID_FORMAT_DEFAULT}}"
SKIP_TO_STR="${SKIP_TO_STR:-${SKIP_TO_STR_DEFAULT}}"
STOP_AT_STR="${STOP_AT_STR:-${STOP_AT_STR_DEFAULT}}"
GZ=${GZ:-${GZ_DEFAULT}}
## derep options
TRUNCATION_LENGTH="${TRUNCATION_LENGTH:-${TRUNCATION_LENGTH_DEFAULT}}"
FASTQ_SUFFIX_PATTERN="${FASTQ_SUFFIX_PATTERN:-${FASTQ_SUFFIX_PATTERN_DEFAULT}}"
## derep dada2 options
maxEE="${maxEE:-${maxEE_DEFAULT}}"
RM_PHIX="${RM_PHIX:-${RM_PHIX_DEFAULT}}"
MULTITHREAD="${MULTITHREAD:-${MULTITHREAD_DEFAULT}}"
## blast options
FA_REF="${FA_REF:-${FA_REF_DEFAULT}}"
GFF_REF="${GFF_REF:-${GFF_REF_DEFAULT}}"
BLAST_SHORT_THRESHOLD="${BLAST_SHORT_THRESHOLD:-${BLAST_SHORT_THRESHOLD_DEFAULT}}"
## read identity options
EXCEL_GENE_PATTERN="${EXCEL_GENE_PATTERN:-${EXCEL_GENE_PATTERN_DEFAULT}}"
## crispresso2 wrapper options
CRISPRESSO_WRAPPER_VERSION="${CRISPRESSO_WRAPPER_VERSION:-${CRISPRESSO_WRAPPER_VERSION_DEFAULT}}"
EDITOR="${EDITOR:-${EDITOR_DEFAULT}}"
## misc paths
CONDA="${CONDA:-${CONDA_DEFAULT}}"

## convert EXCEL to abspath
if [[ ! -z "${EXCEL}"  ]]; then
    EXCEL=$(realpath "${EXCEL}")
fi

## write log
script_path="$SCRIPT_DIR/${SCRIPT_NAME}"
mkdir -p ${DIR}/logfile
printf -- "${params}\n\n${script_path}\n\n## general\n-d|--dir|--output_dir:\t${DIR}\n-p|--prefix:\t${PREFIX}\n--skip-to|--start-at:\t${SKIP_TO_STR}\n--stop-at:\t${STOP_AT_STR}\n--sample-id-format:\t${SAMPLE_ID_FORMAT}\n\n## input\n-f|--fastq|--fastq_dir:\t${DIR_FASTQ}\n--gz:\t${GZ}\n\n## derep\n-t|--truncation:\t${TRUNCATION_LENGTH}\n-s|--fastq-suffix:\t${FASTQ_SUFFIX_PATTERN}\n\n## derep dada2\n--maxEE:\t${maxEE}\n--rm-phix:\t${RM_PHIX}\n-m|--multithread:\t${MULTITHREAD}\n\n## blast\n--assembly:\t${FA_REF}\n--gff:\t${GFF}\n--blast-short-threshold:\t${BLAST_SHORT_THRESHOLD}\n\n## read identity\n-e|--excel|-b|--booking-sheet|--amplicon_book:\t${EXCEL}\n-n|-s|--sheetname|--sheet:\t${EXCEL_SHEETNAME}\n--metadata|--excel-tsv:\t${EXCEL_TSV}\n--excel-gene-pattern:\t${EXCEL_GENE_PATTERN}\n\n## crispresso2 wrapper options\n--crispresso-wrapper-version:\t${CRISPRESSO_WRAPPER_VERSION}\n--editor:\t${EDITOR}\n--amplicon_fasta_dir:\t${DIR_AMPLICON}\n" > ${DIR}/logfile/${PREFIX}_analyseNGS.log

## prep conda
source "$(dirname ${CONDA%*/conda})/etc/profile.d/conda.sh"

## set SKIP_TO and STOP_AT (i.e. manually implement some kind of goto/jump functionality)
declare -A step_map=(
    [start]=0
    [parse-excel]=1
    [derep]=2
    [blast]=3
    [intersect]=4
    [id-seq]=5
    [id-read]=6
    [demultiplex]=7
    [crispresso]=8
    [mutation-table]=9
    [end]=100
)
step_keys=( ${!step_map[@]} )
# if [[ ! ${!step_map[@]} =~ ${SKIP_TO_STR} ]]; then
if [[ ! ${step_keys[@]} =~ ${SKIP_TO_STR} ]]; then
    echo "Invalid --skip-to input. Valid inputs: derep, blast, intersect, id, demultiplex."
    exit 1
elif [[ ! ${step_keys[@]} =~ ${STOP_AT_STR} ]]; then
    echo "Invalid --stop-at input. Valid inputs: derep, blast, intersect, id, demultiplex."
    exit 1
fi
SKIP_TO=${step_map[${SKIP_TO_STR}]}
STOP_AT=${step_map[${STOP_AT_STR}]}
## throw error if STOP_AT is earlier than SKIP_TO
if [[ ${SKIP_TO} -gt ${STOP_AT} ]]; then
    echo "Invalid --skip-to/--stop-at combination. --skip-to step must be same as or earlier than --stop-at step."
    exit 1
fi
## throw error if crispresso step is to be run but --amplicon_fasta_dir is not used
if ( [[ ${SKIP_TO} -le ${step_map[crispresso]} ]] && [[ ${STOP_AT} -ge ${step_map[crispresso]} ]] &&
         [[ -z "${DIR_AMPLICON}" ]] ); then
    echo "--amplicon_fasta_dir is required if Crispresso2 is to be executed."
    exit 1
fi

## set paths
grp_id=${PREFIX}_trunc${TRUNCATION_LENGTH} ## new prefix after derep step
dir_metadata=${DIR}/metadata
dir_derep=${DIR}/derep/${grp_id} ## generated by derep step
dir_blast=${DIR}/blast/${grp_id} ## first used by blast step
dir_identity=${DIR}/identity ## generated by id-read step
dir_demultiplex=${DIR}/demultiplex/${grp_id} ## first used by demultiplex
dir_crispresso=${DIR}/crispresso/${grp_id} ## first used by crispresso
dir_mutation=${DIR}/mutation/${grp_id} ## first used by mutation
f_tmp=${DIR}/tmp.txt
f_tmp2=${DIR}/tmp2.txt

## set common variables
blast_header="sseqid sstart send qseqid qstart qend length sframe qframe mismatch gaps bitscore"

## parse excel
step="parse-excel"
if [[ -z "${EXCEL_TSV}" ]]; then
    echo "--Parsing Excel file--"
    mkdir -p ${dir_metadata}
    EXCEL_TSV=${dir_metadata}/${PREFIX}.tsv
    conda activate analyse_ngs
    # xlsx2csv -n "${EXCEL_SHEETNAME}" -d tab --ignoreempty --escape "${EXCEL}" | sed 's/\\n//g' > ${EXCEL_TSV}
    ${SCRIPT_DIR}/scripts/booking_sheet_wide2long.py \
                 --booking-excel ${EXCEL} --sheet ${EXCEL_SHEETNAME} --out ${EXCEL_TSV}
    conda deactivate
fi

## derep
step="derep"
if [[ ${SKIP_TO} -le ${step_map[${step}]} ]] && [[ ${STOP_AT} -ge ${step_map[${step}]} ]]; then
    echo "--Dereplicating fastq files--"
    conda activate analyse_ngs
    if [[ -z ${GZ} ]]; then
        ${SCRIPT_DIR}/scripts/derep_fastq.R --fastq ${DIR_FASTQ} --out ${DIR} --prefix ${PREFIX} \
                     --truncation ${TRUNCATION_LENGTH} --fastq-suffix ${FASTQ_SUFFIX_PATTERN} \
                     --maxEE ${maxEE} --rm-phix ${RM_PHIX} --multithread ${MULTITHREAD} \
                     --metadata ${EXCEL_TSV} --fastq2sample ${SCRIPT_DIR}/scripts/fastq2sample.R \
                     --sample-id-format ${SAMPLE_ID_FORMAT}
    else
        ${SCRIPT_DIR}/scripts/derep_fastq.R --fastq ${DIR_FASTQ} --out ${DIR} --prefix ${PREFIX} \
                     --truncation ${TRUNCATION_LENGTH} --fastq-suffix ${FASTQ_SUFFIX_PATTERN} \
                     --maxEE ${maxEE} --rm-phix ${RM_PHIX} --multithread ${MULTITHREAD} \
                     --metadata ${EXCEL_TSV} --fastq2sample ${SCRIPT_DIR}/scripts/fastq2sample.R \
                     --sample-id-format ${SAMPLE_ID_FORMAT} \
                     --gz
    fi
    conda deactivate
fi

## blast
step="blast"
if [[ ${SKIP_TO} -le ${step_map[${step}]} ]] && [[ ${STOP_AT} -ge ${step_map[${step}]} ]]; then
    echo "--Executing BLASTN--"
    # ## execute blast + identify top hit's gene
    # source /mnt/chaelab/rachelle/src/run_blast6.sh
    ## set input & output directories
    dir_blast_inpt=${dir_derep}
    ## make output directory
    mkdir -p ${dir_blast}
    ## prepare parameters
    if [[ ${TRUNCATION_LENGTH} < ${BLAST_SHORT_THRESHOLD} ]]; then
        task=blastn-short
    else
        task=blastn
    fi
    ## blast
    source ${SCRIPT_DIR}/scripts/fastq2sample.sh
    fwd_rvs=( fwd rvs )
    while read -r fastq_id; do
        sample_id=$(fastq2sample ${SAMPLE_ID_FORMAT} ${fastq_id})
        for read_type in ${fwd_rvs[@]}; do
            fasta=${dir_blast_inpt}/${sample_id}.derep.${read_type}.fasta
            base=$(basename ${fasta%*.fasta})
            echo ${base}
            tsv_blast=${dir_blast}/${base}.blastn.tsv
            blastn -query ${fasta} -subject ${FA_REF} \
                   -out ${f_tmp} -task ${task} \
                   -outfmt "6 ${blast_header}"
            ## modify blast output so that smaller coord is always in 'start' col (and also 0-index it)
            python3 -c "f = open('${f_tmp}', 'r'); dat = [x.split('\t') for x in f.readlines()]; f.close(); dat = [x[:1] + sorted(map(int, x[1:3])) + x[3:] for x in dat]; dat = [x[:1] + [str(x[1]-1), str(x[2])] + x[3:] for x in dat]; f = open('${tsv_blast}', 'w+'); f.write(''.join('\t'.join(x) for x in dat)); f.close()"
            ## remove tmp files
            rm ${f_tmp}
        done
    done <<< "$(awk 'NR>1 {print $2}' ${EXCEL_TSV} | cut -f2 | sort | uniq)"
fi

## intersect (overwrites original tsv_blast)
step="intersect"
if [[ ${SKIP_TO} -le ${step_map[${step}]} ]] && [[ ${STOP_AT} -ge ${step_map[${step}]} ]]; then
    echo "--Intersecting BLAST hits with annotations--"
    conda activate analyse_ngs
    source ${SCRIPT_DIR}/scripts/fastq2sample.sh
    fwd_rvs=( fwd rvs )
    while read -r fastq_id; do
    # for tsv_blast in ${dir_blast}/*.blastn.tsv; do
        sample_id=$(fastq2sample ${SAMPLE_ID_FORMAT} ${fastq_id})
        for read_type in ${fwd_rvs[@]}; do
            base=${sample_id}.derep.${read_type}
            tsv_blast=${dir_blast}/${base}.blastn.tsv
            ## get coordinates of overlapping genes in BED format
            bedtools intersect -u -b ${tsv_blast} -a ${GFF_REF} |
                awk '$3=="gene"' | gff2bed | cut -f-4 > ${f_tmp}
            ## generate final intersect output
            printf "${blast_header} gchrom gstart gend gene\n" | tr ' ' '\t' > ${f_tmp2}
            bedtools intersect -loj -a ${tsv_blast} -b ${f_tmp} >> ${f_tmp2}
            ## move output to final destination
            mv ${f_tmp2} ${dir_blast}/${base}.blastn.intersectGene.tsv
            ## remove tmp files
            rm ${f_tmp}
        done
    done <<< "$(awk 'NR>1 {print $2}' ${EXCEL_TSV} | cut -f2 | sort | uniq)"
    conda deactivate
fi

## get identity of truncated derep sequences
step="id-seq"
if [[ ${SKIP_TO} -le ${step_map[${step}]} ]] && [[ ${STOP_AT} -ge ${step_map[${step}]} ]]; then
    echo "--Identifying gene(s) most similar to each dereplicated sequence--"
    conda activate analyse_ngs
    source ${SCRIPT_DIR}/scripts/fastq2sample.sh
    fwd_rvs=( fwd rvs )
    while read -r fastq_id; do
    # for tsv_blast in ${dir_blast}/*.blastn.intersectGene.tsv; do
        sample_id=$(fastq2sample ${SAMPLE_ID_FORMAT} ${fastq_id})
        for read_type in ${fwd_rvs[@]}; do
            base=${sample_id}.derep.${read_type}
            tsv_blast=${dir_blast}/${base}.blastn.intersectGene.tsv
            ${SCRIPT_DIR}/scripts/get_derep_seq_identity.R --collapse \
                         --in ${tsv_blast} \
                         --out ${dir_derep}/${base}.identity.txt
        done
    done <<< "$(awk 'NR>1 {print $2}' ${EXCEL_TSV} | cut -f2 | sort | uniq)"
    conda deactivate
fi

## get identity of reads
step="id-read"
if [[ ${SKIP_TO} -le ${step_map[${step}]} ]] && [[ ${STOP_AT} -ge ${step_map[${step}]} ]]; then
    echo "--Tabulating read identities--"
    conda activate analyse_ngs
    ${SCRIPT_DIR}/scripts/get_read_identity.R --prefix ${grp_id} \
                 --out ${DIR} --derep ${dir_derep} --blast ${dir_blast} \
                 --metadata ${EXCEL_TSV} --gene-pattern ${EXCEL_GENE_PATTERN}
    conda deactivate
fi

## demultiplex, use XXX.identity-ontarget.txt from id-read to determine which genes to keep per fastq file
step="demultiplex"
if [[ ${SKIP_TO} -le ${step_map[${step}]} ]] && [[ ${STOP_AT} -ge ${step_map[${step}]} ]]; then
    echo "--Demultiplexing fastq files--"
    conda activate analyse_ngs
    declare -A fastq_suffixes=(
        [fwd]=$(printf ${FASTQ_SUFFIX_PATTERN} | sed 's/R\\d/R1/g').fastq
        [rvs]=$(printf ${FASTQ_SUFFIX_PATTERN} | sed 's/R\\d/R2/g').fastq
    )
    dir_tmp=${DIR}/tmp
    mkdir -p ${dir_demultiplex} ${dir_tmp}
    ## import functions
    source ${SCRIPT_DIR}/scripts/fastq2sample.sh
    ## set common variables
    f_id=${dir_identity}/${grp_id}.identity-ontarget.txt
    ## iterate through fastq-id/gene-id combos
    while read -r fastq_id gene_id; do
        echo "${fastq_id} ${gene_id}"
        ## set variables shared by given fastq-id/gene-id combo
        sample_id=$(fastq2sample ${SAMPLE_ID_FORMAT} ${fastq_id})
        declare -A fastq_filtered=(
            [fwd]=${dir_tmp}/${sample_id}.${gene_id}.fwd.fastq.gz
            [rvs]=${dir_tmp}/${sample_id}.${gene_id}.rvs.fastq.gz
        )
        ## process forward/fwd/R1 and reverse/rvs/R2 reads
        for read_type in ${!fastq_suffixes[@]}; do
            tmp_prefix=${sample_id}.derep.${read_type}
            fastq=${DIR_FASTQ}/${sample_id}${fastq_suffixes[${read_type}]}
            fasta=${dir_derep}/${tmp_prefix}.fasta
            f_map=${dir_derep}/${tmp_prefix}.identity.txt
            f_patterns_for_seqkit=${f_tmp}
            ## generate file for sequence patterns for seqkit
            echo $fastq_id
            ${SCRIPT_DIR}/scripts/filter_derep_by_gene.py --fasta ${fasta} --map ${f_map} \
                         --gene ${gene_id} --out ${f_patterns_for_seqkit}
            ## filter for sequences that best match the given gene with seqkit
            if [[ -z ${GZ} ]]; then
                cat ${fastq} |
                    seqkit grep --by-seq --only-positive-strand \
                           --region 1:${TRUNCATION_LENGTH} \
                           --pattern-file ${f_patterns_for_seqkit} \
                           --out-file ${fastq_filtered[${read_type}]}
            else
                zcat ${fastq}.gz |
                    seqkit grep --by-seq --only-positive-strand \
                           --region 1:${TRUNCATION_LENGTH} \
                           --pattern-file ${f_patterns_for_seqkit} \
                           --out-file ${fastq_filtered[${read_type}]}
            fi
            ## remove tmp files
            rm ${f_patterns_for_seqkit}
        done
        ## filter seqkit output for sequences present in both fwd and rvs files
        seqkit pair --read1 ${fastq_filtered[fwd]} --read2 ${fastq_filtered[rvs]} \
               --out-dir ${dir_demultiplex}
        ## remove tmp files
        rm ${fastq_filtered[fwd]} ${fastq_filtered[rvs]}
    done <<< "$(awk -F '\t' 'NR>1 && $3=="fwd_rvs" {print $1,$2}' ${f_id})"
    ## remove tmp directory
    rmdir ${dir_tmp}
    conda deactivate
fi

step="crispresso"
if [[ ${SKIP_TO} -le ${step_map[${step}]} ]] && [[ ${STOP_AT} -ge ${step_map[${step}]} ]]; then
    echo "--Executing Crispresso2 on demultiplexed files--"
    if [[ "${EDITOR}" == 'ABE' ]]; then
        preset_args='--cleavage_offset -12 --quantification_window_size 6'
    elif [[ "${EDITOR}" == 'unknown' ]]; then
        preset_args='--screen-grna'
    fi
    conda activate crispresso2_env
    mkdir -p ${dir_crispresso}
    ${SCRIPT_DIR}/scripts/Crispresso2_wrapper/Crispresso2_wrapper_v${CRISPRESSO_WRAPPER_VERSION}.py \
                 --excel-tsv ${EXCEL_TSV} --PCR_Product_fa_gene_pattern "${EXCEL_GENE_PATTERN_DEFAULT}" \
                 --fastq_dir ${dir_demultiplex} --amplicon_fasta_dir ${DIR_AMPLICON} \
                 --output_dir ${dir_crispresso} ${preset_args}
    conda deactivate
fi

step="mutation-table"
if [[ ${SKIP_TO} -le ${step_map[${step}]} ]] && [[ ${STOP_AT} -ge ${step_map[${step}]} ]]; then
    echo "--Making mutation table--"
    mkdir -p ${dir_mutation}
    fout=${dir_mutation}/${PREFIX}.Quantification_window_nucleotide_percentage_table.long.txt
    EXCEL_TSV=${dir_metadata}/${PREFIX}.tsv
    conda activate analyse_ngs
    ${SCRIPT_DIR}/scripts/make_mut_table.py \
                 --metadata ${EXCEL_TSV} \
                 --crispresso-dir ${dir_crispresso} \
                 --out ${fout} \
                 --fastq2sample ${SCRIPT_DIR}/scripts/fastq2sample.py \
                 --sample-id-format ${SAMPLE_ID_FORMAT} \
                 --gene-pattern ${EXCEL_GENE_PATTERN}
    conda deactivate
fi

exit 0
