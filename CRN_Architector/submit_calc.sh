
for VAR in $(seq 0 10)
do
    sed -i '/python build_structures_arch.py/c\python build_structures_arch.py -n '"$VAR"' -b 1' batch.sub
    sbatch batch.sub
done
