region_counter
==

A program to quickly count the number of reads in an RNAseq sam/bam/cram file.

Counts are 
  - Mapped (total mapped reads)
  - Mapped, exon (mapped within an exon region)
  - Unmapped

Usage:
```
region_counter [OPTIONS] --bamfile <BAMFILE> --gtf <GTF>

Options:
  -b, --bamfile <BAMFILE>
  -g, --gtf <GTF>
  -q, --minmapqual <MINMAPQUAL>        [default: 35]
  -f, --required-flag <REQUIRED_FLAG>  [default: 3]
  -F, --filtered-flag <FILTERED_FLAG>  [default: 2816]
  -h, --help                           Print help
  -V, --version                        Print version
```

### Filtering reads

Reads contribute to the count if they are greater than or equal to a minimum mapping threshold, if they satisfy all `required-flag` flags (default=3 - include only if read is paired and mapped in proper pair) and have no `filtered-flag` flags (default=2816 - exclude if read is secondary, read fails vendor quality checks, or read is supplementary).
