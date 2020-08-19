from genmod.annotate_models.models import check_recessive
from genmod.vcf_tools import Genotype

from ped_parser import FamilyParser

FAMILY_FILE = "../fixtures/recessive_trio.ped"

def get_family(family_file = None, family_lines = None):
    """Return a family object
    
    """
    family = None
    if family_file:
        family = FamilyParser(open(family_file, 'r'))
    elif family_lines:
        family = FamilyParser(family_lines)
        
    return family

################# Test affected ###############

def test_recessive_affected_homozygote_male():
    """Test an affected homozygote male"""
    
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\t0\t0\t1\t2\n"
    ]
    
    family = get_family(family_lines=family_lines)
    
    homozygote_variant = {'genotypes': {}}
    homozygote_variant['genotypes']['proband'] = Genotype(**{'GT':'1/1'})
    
    assert check_recessive(
        variant = homozygote_variant,
        family = family,
        strict = False
    ) == True

def test_recessive_affected_het_male():
    """Test a sick male
    """
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\t0\t0\t1\t2\n"
    ]
    
    family = get_family(family_lines=family_lines)
    
    recessive_variant = {'genotypes': {}}
    recessive_variant['genotypes']['proband'] = Genotype(**{'GT':'0/1'})
    
    assert check_recessive(
        variant = recessive_variant,
        family = family,
        strict = False
    ) == False

def test_recessive_affected_male_ref_call():
    """Test an affected ref call male"""
    
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\t0\t0\t1\t2\n"
    ]
    
    family = get_family(family_lines=family_lines)
    
    homozygote_variant = {'genotypes': {}}
    homozygote_variant['genotypes']['proband'] = Genotype(**{'GT':'0/0'})
    
    assert check_recessive(
        variant = homozygote_variant,
        family = family,
        strict = False
    ) == False

############### Test family ##############

def test_family_one_affected():
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
    recessive_variant['genotypes']['proband'] = Genotype(**{'GT':'1/1'})
    recessive_variant['genotypes']['father'] = Genotype(**{'GT':'0/1'})
    recessive_variant['genotypes']['mother'] = Genotype(**{'GT':'0/0'})
    
    assert check_recessive(
        variant = recessive_variant,
        family = family,
        strict = False
    ) == True