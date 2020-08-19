from genmod.annotate_models.models import check_compounds
from genmod.vcf_tools import Genotype

from ped_parser import FamilyParser

FAMILY_FILE = "tests/fixtures/recessive_trio.ped"

def get_family(family_file = None, family_lines = None):
    """Return a family object
    
    """
    family = None
    if family_file:
        family = FamilyParser(open(family_file, 'r'))
    elif family_lines:
        family = FamilyParser(family_lines)
        
    return family

# Affected
### Test family with het healthy parents and het compound proband affected
def test_2_het_parents_affected():
    """Test healthy parents and an affected proband
    """
    # family = get_family(family_file=FAMILY_FILE)
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\tfather\tmother\t1\t2\n",
        "1\tmother\t0\t0\t2\t1\n",
        "1\tfather\t0\t0\t1\t1\n"
    ]
    family = get_family(family_lines=family_lines)

    recessive_variant = {'genotypes': {}}
    recessive_variant['genotypes']['proband'] = Genotype(**{'GT':'0/1'})
    recessive_variant['genotypes']['father'] = Genotype(**{'GT':'0/1'})
    recessive_variant['genotypes']['mother'] = Genotype(**{'GT':'0/0'})

    recessive_variant_2 = {'genotypes': {}}
    recessive_variant_2['genotypes']['proband'] = Genotype(**{'GT':'0/1'})
    recessive_variant_2['genotypes']['father'] = Genotype(**{'GT':'0/0'})
    recessive_variant_2['genotypes']['mother'] = Genotype(**{'GT':'0/1'})
    
    assert check_compounds(
        variant_1 = recessive_variant,
        variant_2 = recessive_variant_2,
        family = family,
        intervals = dict(),
        phased = False
    ) == True

### Test family with one het healthy parent and het compound proband affected, with dn variant
def test_1_het_parents_affected():
    """Test healthy parents and an affected proband
    """
    # family = get_family(family_file=FAMILY_FILE)
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\tfather\tmother\t1\t2\n",
        "1\tmother\t0\t0\t2\t1\n",
        "1\tfather\t0\t0\t1\t1\n"
    ]
    family = get_family(family_lines=family_lines)

    recessive_variant = {'genotypes': {}}
    recessive_variant['genotypes']['proband'] = Genotype(**{'GT':'0/1'})
    recessive_variant['genotypes']['father'] = Genotype(**{'GT':'0/0'})
    recessive_variant['genotypes']['mother'] = Genotype(**{'GT':'0/0'})

    recessive_variant_2 = {'genotypes': {}}
    recessive_variant_2['genotypes']['proband'] = Genotype(**{'GT':'0/1'})
    recessive_variant_2['genotypes']['father'] = Genotype(**{'GT':'0/0'})
    recessive_variant_2['genotypes']['mother'] = Genotype(**{'GT':'0/1'})
    
    assert check_compounds(
        variant_1 = recessive_variant,
        variant_2 = recessive_variant_2,
        family = family,
        intervals = dict(),
        phased = False
    ) == True
    
# Healthy
### Test family with one het healthy parent for both variants
def test_1_het_parents_affected():
    """Test healthy parents and an affected proband
    """
    # family = get_family(family_file=FAMILY_FILE)
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\tfather\tmother\t1\t2\n",
        "1\tmother\t0\t0\t2\t1\n",
        "1\tfather\t0\t0\t1\t1\n"
    ]
    family = get_family(family_lines=family_lines)

    recessive_variant = {'genotypes': {}}
    recessive_variant['genotypes']['proband'] = Genotype(**{'GT':'0/1'})
    recessive_variant['genotypes']['father'] = Genotype(**{'GT':'0/1'})
    recessive_variant['genotypes']['mother'] = Genotype(**{'GT':'0/0'})

    recessive_variant_2 = {'genotypes': {}}
    recessive_variant_2['genotypes']['proband'] = Genotype(**{'GT':'0/1'})
    recessive_variant_2['genotypes']['father'] = Genotype(**{'GT':'0/1'})
    recessive_variant_2['genotypes']['mother'] = Genotype(**{'GT':'0/1'})
    
    assert check_compounds(
        variant_1 = recessive_variant,
        variant_2 = recessive_variant_2,
        family = family,
        intervals = dict(),
        phased = False
    ) == False