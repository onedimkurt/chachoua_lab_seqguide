#!/bin/bash
# Checks the genome-browser files written by SeqGuide (pure Python) with the reference tool, samtools.
# Runs inside the SeqGuide container:  docker compose run --rm seqguide bash /data/check_exports.sh
cd /data/runs/validation || { echo "no validation results"; exit 1; }
fail=0
for d in */genome_browser; do
  s=${d%%/*}
  echo "== $s"
  fa=$(ls "$d"/*.fasta | head -1)
  name=$(cut -f1 "$fa.fai" | head -1)
  got=$(samtools faidx "$fa" "$name:1-30" | tail -n +2 | tr -d '\n')
  want=$(grep -v '^>' "$fa" | tr -d '\n' | cut -c1-30)
  [ "$got" = "$want" ] && echo "  fasta .fai index OK" || { echo "  fasta .fai index MISMATCH"; fail=1; }
  if [ -f "$d/reads.bam" ]; then
    if samtools quickcheck -v "$d/reads.bam" && samtools view "$d/reads.bam" > /dev/null; then
      echo "  BAM readable (quickcheck + full decode) OK"
    else
      echo "  BAM FAILED"; fail=1
    fi
    total=$(samtools view -c "$d/reads.bam")
    idx=$(samtools idxstats "$d/reads.bam" | awk '{s+=$3} END {print s}')
    [ "$total" = "$idx" ] && echo "  .bai index counts OK ($total reads)" || { echo "  .bai counts MISMATCH ($total vs $idx)"; fail=1; }
    a=$(samtools depth -a "$d/reads.bam" | awk -v n="$name" '$1==n && $2>=1000 && $2<=1020 {print $3}' | tr '\n' ' ')
    b=$(samtools depth -a -r "$name:1000-1020" "$d/reads.bam" | awk '{print $3}' | tr '\n' ' ')
    [ "$a" = "$b" ] && echo "  indexed region query OK (depth at 1000-1020: ${a%% *}...)" || { echo "  indexed query MISMATCH"; fail=1; }
  else
    echo "  (no BAM in this folder)"
  fi
done
[ $fail = 0 ] && echo "ALL GENOME-BROWSER FILES OK" || echo "SOME GENOME-BROWSER FILES FAILED"
