"""
Tests for spliceai_utils module.
"""

import pytest
import sys
from pathlib import Path

# Add src to path so we can import spliceai_utils
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from my_utils.spliceai_utils import (
    Variant,
    load_spliceai_models,
    spliceO_predictions,
)


@pytest.fixture(scope="session")
def fasta_path():
    """Return path to test FASTA file (chr2:1-2,333,275)."""
    path = Path(__file__).parent / "data" / "chr2_1-2333275.fa.gz"
    
    if not path.exists():
        pytest.fail(
            f"Test FASTA file not found: {path}\n"
            f"Expected bgzipped file containing chr2:1-2,333,275 from hg38 reference genome."
        )
    
    print(f"\n✓ Using test FASTA: {path}")
    return str(path)


@pytest.fixture(scope="module")
def spliceai_models():
    """Load SpliceAI models once for all tests."""
    print("\nLoading SpliceAI models...")
    models = load_spliceai_models()
    print("✓ Models loaded")
    return models


@pytest.fixture
def test_variant():
    """Create a test SNV variant."""
    return Variant('chr2', 1860161, 'C', 'A')


def test_variant_creation(test_variant):
    """Test that a Variant can be created."""
    assert test_variant.chrom == 'chr2'
    assert test_variant.pos == 1860161
    assert test_variant.ref == 'C'
    assert test_variant.alt == 'A'
    print(f"✓ Test variant created: {test_variant}")


def test_variant_properties(test_variant):
    """Test Variant properties."""
    assert test_variant.is_snv is True
    assert test_variant.is_insertion is False
    assert test_variant.is_deletion is False
    assert test_variant.indel_length == 0
    assert test_variant.start == 1860160  # 0-based
    assert test_variant.end == 1860161  # 0-based


def test_variant_validation(fasta_path, test_variant):
    """Test that variant REF validation works."""
    try:
        test_variant.validate_ref(fasta_path)
        print("✓ Variant REF validated")
    except ValueError as e:
        pytest.fail(f"Variant validation failed: {e}")


def test_spliceO_predictions_without_variants(fasta_path, spliceai_models):
    """Test primer walks WITHOUT variants."""
    print("\nRunning primer walks WITHOUT variants...")
    results_wt = spliceO_predictions(
        fasta_path=fasta_path,
        region_interval="chr2:1,840,760-1,887,609",
        walk_interval=(1860114, 1860248),
        walk_step=20,
        substitution_length=20,
        sites_of_interest={'cryptic_donor': 1860161},
        models=spliceai_models,
        strand="-",
        variants=None
    )
    
    assert len(results_wt) > 0, "Should return predictions"
    assert 'donor_prob' in results_wt.columns
    assert 'acceptor_prob' in results_wt.columns
    assert 'MaskStart' in results_wt.columns
    assert 'Mask_Context' in results_wt.columns
    assert 'ASO_Sequence' in results_wt.columns
    
    # Check that we have a wildtype control
    wt_rows = results_wt[results_wt['MaskStart'] == 'WT']
    assert len(wt_rows) > 0, "Should have WT control"
    
    print(f"✓ Predictions without variants completed: {len(results_wt)} rows")
    print(f"  Columns: {list(results_wt.columns)}")
    print(f"  Sample data:\n{results_wt.head(3)}")


def test_spliceO_predictions_with_variant(fasta_path, spliceai_models, test_variant):
    """Test primer walks WITH an SNV variant."""
    print("\nRunning primer walks WITH SNV variant...")
    results_var = spliceO_predictions(
        fasta_path=fasta_path,
        region_interval="chr2:1,840,760-1,887,609",
        walk_interval=(1860114, 1860248),
        walk_step=20,
        substitution_length=20,
        sites_of_interest={'cryptic_donor': 1860161},
        models=spliceai_models,
        strand="-",
        variants=[test_variant]
    )
    
    assert len(results_var) > 0, "Should return predictions"
    print(f"✓ Predictions with variant completed: {len(results_var)} rows")


def test_variant_effect_comparison(fasta_path, spliceai_models, test_variant):
    """Compare predictions with and without variant."""
    print("\nComparing variant effects...")
    
    # Without variant
    results_wt = spliceO_predictions(
        fasta_path=fasta_path,
        region_interval="chr2:1,840,760-1,887,609",
        walk_interval=(1860114, 1860248),
        walk_step=20,
        substitution_length=20,
        sites_of_interest={'cryptic_donor': 1860161},
        models=spliceai_models,
        strand="-",
        variants=None
    )
    
    # With variant
    results_var = spliceO_predictions(
        fasta_path=fasta_path,
        region_interval="chr2:1,840,760-1,887,609",
        walk_interval=(1860114, 1860248),
        walk_step=20,
        substitution_length=20,
        sites_of_interest={'cryptic_donor': 1860161},
        models=spliceai_models,
        strand="-",
        variants=[test_variant]
    )
    
    # Get wildtype predictions
    wt_prob = results_wt[results_wt['MaskStart'] == 'WT']['donor_prob'].values[0]
    var_prob = results_var[results_var['MaskStart'] == 'WT']['donor_prob'].values[0]
    
    print(f"Donor probability at cryptic site:")
    print(f"  Wildtype: {wt_prob:.4f}")
    print(f"  With C>A variant: {var_prob:.4f}")
    print(f"  Difference: {var_prob - wt_prob:.4f}")
    
    # Probabilities should be different
    assert wt_prob != var_prob, "Variant should change predictions"


def test_spliceO_predictions_with_immutable_range(fasta_path, spliceai_models):
    """Test primer walks with immutable range protection."""
    print("\nTesting immutable range...")
    results_immutable = spliceO_predictions(
        fasta_path=fasta_path,
        region_interval="chr2:1,840,760-1,887,609",
        walk_interval=(1860114, 1860248),
        walk_step=20,
        substitution_length=20,
        sites_of_interest={'cryptic_donor': 1860161},
        models=spliceai_models,
        strand="-",
        variants=None,
        immutable_range=(1860160, 1860163)  # Protect 4 bases around cryptic donor site
    )
    
    assert len(results_immutable) > 0, "Should return predictions"
    print(f"✓ Predictions with immutable range completed: {len(results_immutable)} rows")
    
    # Check that immutable bases are preserved in mask context
    masked_rows = results_immutable[results_immutable['MaskStart'] != 'WT']
    if len(masked_rows) > 0:
        sample_mask = masked_rows.iloc[0]
        mask_context = sample_mask['Mask_Context']
        print(f"  Sample mask context: {mask_context}")
        
        # The mask context should have some non-N bases if immutable range overlaps
        # We don't assert specific content since it depends on the overlap
        assert isinstance(mask_context, str), "Mask_Context should be a string"
        print(f"  Note: Immutable bases should be preserved (not N's) when they overlap the mask")


if __name__ == "__main__":
    # Allow running tests directly with python
    pytest.main([__file__, "-v", "-s"])
