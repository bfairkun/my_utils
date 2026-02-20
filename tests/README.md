# Testing Guide

## Running Tests

### Run all tests:

```bash
pytest
```

### Run with verbose output:

```bash
pytest -v
```

### Run with output from print statements:

```bash
pytest -v -s
```

### Run a specific test file:

```bash
pytest tests/test_spliceai_utils.py
```

### Run a specific test function:

```bash
pytest tests/test_spliceai_utils.py::test_variant_creation
```

### Run with coverage report:

```bash
pytest --cov=my_utils --cov-report=term-missing
```

### Run tests in parallel (if you have pytest-xdist installed):

```bash
pytest -n auto
```

## Reference Genome Setup

The `test_spliceai_utils.py` tests require a reference genome FASTA file.

### Default Path

By default, tests look for:

```
/project2/yangili1/bjf79/ReferenceGenomes/Human_UCSC.hg38_GencodeComprehensive46/Reference.fa
```

### Custom Path

Set a custom path using an environment variable:

```bash
export SPLICEAI_TEST_FASTA=/path/to/your/reference.fa
pytest tests/test_spliceai_utils.py
```

### Download Reference Genome

If you don't have a reference genome, download one:

**Option 1: UCSC hg38**

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
export SPLICEAI_TEST_FASTA=./hg38.fa
```

**Option 2: GENCODE (includes gene annotations)**

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
export SPLICEAI_TEST_FASTA=./GRCh38.primary_assembly.genome.fa
```

### Checksum Validation

The tests compute a checksum of the first 10MB of the FASTA file for validation.
On first run, the checksum will be printed. You can add it to
`EXPECTED_FASTA_SHA256` in the test file for future validation.

If tests are skipped due to missing FASTA, you'll see helpful error messages
with download instructions.

## Installing Test Dependencies

Make sure you have the test dependencies installed:

```bash
pip install -e ".[dev]"
```

## Test Structure

- `tests/test_package.py` - Basic package tests
- `tests/test_spliceai_utils.py` - Tests for spliceai_utils module

## Notes

- Tests use pytest fixtures to load SpliceAI models once per session for
  efficiency
- FASTA file is validated with checksum for consistency
- Tests are automatically skipped if reference genome is not available
- Use `-s` flag to see print statements during test execution
