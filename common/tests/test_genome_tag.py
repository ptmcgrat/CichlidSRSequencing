"""cichlid_common.common_light.genome_tag: the suffix of the files and work folders of a build on a genome other than the default."""
import unittest
from cichlid_common.common_light import DEFAULT_GENOME, genome_tag


class TestGenomeTag(unittest.TestCase):
    def test_the_default_genome_keeps_every_name_it_has(self):
        self.assertEqual(DEFAULT_GENOME, 'Mzebra_GT3_NCBI')
        for g in ('Mzebra_GT3_NCBI', None, ''):
            self.assertEqual(genome_tag(g), '')

    def test_another_genome_gets_a_dot_suffix(self):
        self.assertEqual(genome_tag('Mconophoros_GT1_NCBI'), '.Mconophoros_GT1_NCBI')
        self.assertEqual(f"bc_family{genome_tag('Mconophoros_GT1_NCBI')}.snps.vcf.gz", 'bc_family.Mconophoros_GT1_NCBI.snps.vcf.gz')
        self.assertEqual(f"bc_family{genome_tag('Mzebra_GT3_NCBI')}.snps.vcf.gz", 'bc_family.snps.vcf.gz')

    def test_the_default_can_be_given(self):
        self.assertEqual(genome_tag('X', default='X'), '')
        self.assertEqual(genome_tag('Mzebra_GT3_NCBI', default='X'), '.Mzebra_GT3_NCBI')


if __name__ == '__main__':
    unittest.main()
