{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
START_TIME=$(date +%s)
save START_TIME

NUM_THREADS={{ opt.parallel }}
save NUM_THREADS

# Add masurca to $PATH
export PATH="$(readlinkf "$(which masurca)" | xargs dirname):$PATH"

#----------------------------#
# Read stats of PE reads
#----------------------------#
log_info Symlink/copy input files
if [ ! -e pe.cor.fa ]; then
    ln -s {{ args.0 }} pe.cor.fa
fi
cp {{ args.1 }} env.json

log_info Read stats of PE reads

SUM_COR=$( faops n50 -H -N 0 -S pe.cor.fa )
save SUM_COR

KMER='{{ opt.kmer | replace(from=" ", to=",") }}'
save KMER
log_debug "You set kmer size of $KMER for the graph"

{% if opt.estsize == "auto" -%}
ESTIMATED_GENOME_SIZE=$( cat env.json | jq '.ESTIMATED_GENOME_SIZE | tonumber' )
{% else -%}
ESTIMATED_GENOME_SIZE={{ opt.estsize }}
save ESTIMATED_GENOME_SIZE
{% endif -%}
log_debug "ESTIMATED_GENOME_SIZE: $ESTIMATED_GENOME_SIZE"

#----------------------------#
# Build unitigs
#----------------------------#
if [ ! -e unitigs.fasta ]; then
log_info Creating unitigs

{% set kmers = opt.kmer | split(pat=" ") %}
{% for kmer in kmers -%}
    log_debug with k={{ kmer }}
{% if opt.unitigger == "tadpole" -%}
    tadpole.sh \
        in=pe.cor.fa \
        out=unitigs_K{{ kmer }}.fasta \
        threads={{ opt.parallel }} \
        k={{ kmer }} \
        overwrite
{% elif opt.unitigger == "bcalm" -%}
    bcalm -in pe.cor.fa \
        -kmer-size {{ kmer }} -abundance-min 3 -verbose 0 \
        -nb-cores {{ opt.parallel }} -out K{{ kmer }}
    mv K{{ kmer }}.unitigs.fa unitigs_K{{ kmer }}.fasta
{% elif opt.unitigger == "bifrost" -%}
    Bifrost build --input-seq-file pe.cor.fa  \
        --kmer-length {{ kmer }} --clip-tips --del-isolated \
        --threads {{ opt.parallel }} --fasta --output-file unitigs_K{{ kmer }}
{% else -%}
    create_k_unitigs_large_k -c $(({{ kmer }}-1)) -t {{ opt.parallel }} \
        -m {{ kmer }} -n $ESTIMATED_GENOME_SIZE -l {{ kmer }} -f 0.000001 \
        pe.cor.fa \
        > unitigs_K{{ kmer }}.fasta
{% endif -%}
{% endfor -%}
{# Keep a blank line #}
log_info Creating non-contained unitigs
    anchr contained \
{% for kmer in kmers -%}
        unitigs_K{{ kmer }}.fasta \
{% endfor -%}
        --len {{ opt.min }} --idt 0.9999 --ratio 0.99999 --parallel {{ opt.parallel }} \
        -o unitigs.non-contained.fasta

if [ -s unitigs.non-contained.fasta ]; then
{% if opt.merge == "1" -%}
    log_info Merging unitigs
    anchr orient unitigs.non-contained.fasta \
        --len {{ opt.min }} --idt 0.999 --parallel {{ opt.parallel }} \
        -o unitigs.orient.fasta
    anchr merge unitigs.orient.fasta \
        --len {{ opt.min }} --idt 0.9999 --parallel {{ opt.parallel }} \
        -o unitigs.fasta
{% else -%}
    mv unitigs.non-contained.fasta unitigs.fasta
{% endif -%}
else
    touch unitigs.fasta
fi

rm unitigs_K*.fasta

fi

#----------------------------#
# Done.
#----------------------------#
END_TIME=$(date +%s)
save END_TIME

RUNTIME=$((END_TIME-START_TIME))
save RUNTIME

log_info Done.

exit 0
