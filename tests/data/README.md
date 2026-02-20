# Test Data

## chr2_1-2333275.fa.gz

This file contains the first 2,333,275 bases of chromosome 2 from the hg38
reference genome, bgzipped for efficient storage and random access.

**Source**: UCSC hg38 reference genome **Region**: chr2:1-2,333,275 **Size**:
~737 KB (bgzipped)

This region is used for testing the `spliceai_utils` module, specifically for
the test variant at chr2:1,860,161 and surrounding splicing analysis.

### How it was created

```bash
# Extract region from reference
samtools faidx /path/to/Reference.fa chr2:1-2333275 > chr2_1-2333275.fa

# Fix header to be just ">chr2" (required for pysam)
sed -i '1s/>chr2:1-2333275/>chr2/' chr2_1-2333275.fa

# Compress with bgzip for efficient random access
bgzip chr2_1-2333275.fa

# Create index
samtools faidx chr2_1-2333275.fa.gz
```

### Files

- `chr2_1-2333275.fa.gz` - Bgzipped FASTA file with sequence data
- `chr2_1-2333275.fa.gz.fai` - FASTA index file (for samtools/pysam)
- `chr2_1-2333275.fa.gz.gzi` - Bgzip index file (for random access)

### Why bgzip?

- **Much smaller**: 737 KB vs 2.3 MB (68% smaller)
- **Git-friendly**: Smaller diffs when tracking in version control
- **Transparent**: pysam/samtools can read bgzipped files directly
- **Random access**: With `.gzi` index, can efficiently access any region
