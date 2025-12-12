# BF528_Final

## SRA Toolkit Download Workflow

Follow these steps to download sequencing data from the SRA and convert it to compressed FASTQ files.

1. **Prepare the environment**: Install the SRA Toolkit for your platform (see the [NCBI download page](https://github.com/ncbi/sra-tools/wiki/Downloads)). Make sure the `prefetch`, `fastq-dump`, and `vdb-config` executables are on your `PATH`.
2. **Configure cache location (important)**: Run the interactive configuration tool so you can place the cache somewhere with sufficient space instead of the default home directory location:
   ```bash
   vdb-config --interactive
   ```
3. **Download SRA runs**:
   * Pre-download `.sra` files listed in `SRR_Acc_List.txt`:
     ```bash
     prefetch --option-file SRR_Acc_List.txt
     ```
   * Convert each run to compressed FASTQ. The `--split-files` flag is essential for paired-end data:
     ```bash
     cat SRR_Acc_List.txt | xargs -n 1 fastq-dump --split-files --gzip
     ```

Keep `SRR_Acc_List.txt` in the working directory and list one accession per line.
