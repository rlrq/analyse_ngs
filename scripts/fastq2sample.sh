#!/usr/bin

fastq2sample() {
    local fmt=${1}
    local fastq_id=${2}
    if [[ ${fmt} -eq 1 ]]; then
        output=$(echo "${fastq_id}_${fastq_id}" | sed 's/_S0/_S/')
    else
        output=${fastq_id}
    fi
    echo "${output}"
}
