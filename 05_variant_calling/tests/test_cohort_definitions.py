"""Cohort definitions: loading, validation, and the members a definition selects."""
import os, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_fixtures import *
from variant_calling import cohort as C

YAML = """
defaults:
  gatk: {heterozygosity: 0.00175}
  cohort:
    genome: Other_Genome
    caller: bcftools
    groups: per_sample
callsets:
  one:
    kind: family
    broods: [B1]
  two:
    kind: sex_set
    species: X y
    caller: gatk
    groups: species
"""


class TestLoading(unittest.TestCase):
    def test_defaults_are_filled_and_an_entry_can_override_them(self):
        p = tempfile.mktemp(suffix='.yaml'); open(p, 'w').write(YAML)
        c = C.load_cohorts(p)
        self.assertEqual(c['one']['genome'], 'Other_Genome'); self.assertEqual(c['one']['caller'], 'bcftools'); self.assertEqual(c['one']['groups'], 'per_sample')
        self.assertEqual((c['two']['caller'], c['two']['groups']), ('gatk', 'species')); self.assertEqual(c['one']['name'], 'one')
        self.assertEqual(c['one']['exclude_ecogroups'], ['Error', 'Exclude'], 'the module default stays when the file does not set it')

    def test_the_shipped_file_loads_and_every_entry_has_a_kind(self):
        c = C.load_cohorts()
        self.assertTrue({'bc_family', 'yh_broods', 'mc_broods', 'gertrudae_orange_collar', 'saulosi_green_face_peacock'} <= set(c))
        self.assertFalse([n for n in c if 'plus_adults' in n], 'the adult arms were dropped: Category is not used')
        self.assertEqual(c['gertrudae_orange_collar']['kind'], 'sex_set'); self.assertEqual(c['saulosi_green_face_peacock']['species'], 'Aulonocara saulosi green face peacock')
        self.assertFalse([n for n, v in c.items() if 'categories' in v], 'no shipped cohort selects by Category')
        self.assertEqual(c['bc_family']['kind'], 'bc_parents'); self.assertEqual(c['yh_broods']['kind'], 'family')
        for name, spec in c.items():
            self.assertEqual(C.validate(spec), [], name)
            self.assertEqual(spec['groups'], 'per_sample', 'lab and farm stock: no shared Hardy-Weinberg prior (decision 63)')
            self.assertEqual(spec['caller'], 'bcftools')

    def test_category_is_not_a_selector_and_the_message_says_so(self):
        msg = ' '.join(C.validate(family_spec(categories=['OtherYHs'])))
        self.assertIn('categories is not a selector', msg); self.assertNotIn('typo', msg)

    def test_validation_names_the_problem(self):
        self.assertEqual(C.validate(family_spec()), [])
        self.assertIn('kind', ' '.join(C.validate(family_spec(kind='trio'))))
        self.assertIn('caller', ' '.join(C.validate(family_spec(caller='freebayes'))))
        self.assertIn('groups', ' '.join(C.validate(family_spec(groups='clade'))))
        self.assertIn('typo', ' '.join(C.validate(family_spec(bloods=['B1']))))
        self.assertIn('selects nothing', ' '.join(C.validate(dict(kind='family', caller='bcftools', groups='per_sample'))))
        self.assertIn('selects nothing', ' '.join(C.validate(dict(kind='family', caller='bcftools', groups='per_sample', categories=['X']))))
        self.assertIn('must be a list', ' '.join(C.validate(family_spec(broods='B1'))))


class TestMembers(unittest.TestCase):
    def setUp(self):
        extra = [fish('ADULT1', 'M', eco='Deep_Benthic'), fish('ADULT2', 'F', eco='Deep_Benthic'), fish('UNSEXED', 'U'), fish('BAD', 'M', eco='Error'), fish('OTHERSP', 'M', species='Other species')]
        self.dt = sheet(brood() + extra)

    def test_a_brood_gives_offspring_then_parents(self):
        ids, notes = C.resolve_members(self.dt, family_spec())
        self.assertEqual(len(ids), 26); self.assertEqual(ids[0], 'B1_m0'); self.assertEqual(ids[-2:], ['S_m', 'D_f'], 'sire then dam, in the order of the sheet'); self.assertEqual(notes, [])

    def test_parents_can_be_left_out(self):
        ids, _ = C.resolve_members(self.dt, family_spec(include_parents=False))
        self.assertEqual(len(ids), 24); self.assertNotIn('S_m', ids)

    def test_species_next_to_broods_is_only_the_declared_species(self):
        ids, _ = C.resolve_members(self.dt, family_spec(species=YH))
        self.assertEqual(len(ids), 26, 'it must not pull in the whole species')

    def test_species_alone_selects_the_species_sexed_and_not_flagged(self):
        ids, notes = C.resolve_members(self.dt, sexset_spec())
        self.assertIn('ADULT1', ids); self.assertIn('B1_m0', ids); self.assertNotIn('UNSEXED', ids); self.assertNotIn('BAD', ids); self.assertNotIn('OTHERSP', ids)
        self.assertEqual(len(notes), 2); self.assertTrue(any('Ecogroup' in n for n in notes) and any('Sex is not M or F' in n for n in notes))

    def test_unsexed_fish_are_kept_when_sexing_is_not_required(self):
        ids, _ = C.resolve_members(self.dt, sexset_spec(require_sexed=False))
        self.assertIn('UNSEXED', ids)

    def test_ecogroups_narrow_a_species_selection(self):
        ids, _ = C.resolve_members(self.dt, sexset_spec(ecogroups=['Deep_Benthic']))
        self.assertEqual(sorted(ids), ['ADULT1', 'ADULT2'])

    def test_species_names_match_whatever_the_case_and_spacing(self):
        dt = sheet([fish('A', 'M', species='Aulonocara  saulosi Green face PEACOCK'), fish('B', 'F', species='aulonocara saulosi green face peacock'), fish('C', 'F', species='Aulonocara saulosi')])
        ids, _ = C.resolve_members(dt, sexset_spec(species='Aulonocara saulosi green face peacock'))
        self.assertEqual(sorted(ids), ['A', 'B'], 'a shorter species name is another species')

    def test_explicit_samples_come_first_are_kept_even_when_unknown_and_are_not_repeated(self):
        ids, _ = C.resolve_members(self.dt, family_spec(samples=['NOT_IN_SHEET', 'B1_f3'], include_parents=False))
        self.assertEqual(ids[:2], ['NOT_IN_SHEET', 'B1_f3']); self.assertEqual(ids.count('B1_f3'), 1); self.assertEqual(len(ids), 25)

    def test_exclude_samples_removes_even_explicit_ones(self):
        ids, _ = C.resolve_members(self.dt, family_spec(exclude_samples=['S_m', 'B1_m0']))
        self.assertNotIn('S_m', ids); self.assertNotIn('B1_m0', ids)

    def test_blank_or_partial_parents_do_not_become_members(self):
        dt = sheet([fish('A', 'M', brood='B9', parents='P1,'), fish('B', 'F', brood='B9', parents=',P2'), fish('C', 'F', brood='B9', parents=''), fish('P1', 'M'), fish('P2', 'F')])
        ids, _ = C.resolve_members(dt, family_spec(broods=['B9']))
        self.assertEqual(sorted(ids), ['A', 'B', 'C', 'P1', 'P2'])

    def test_a_sheet_without_optional_columns_still_resolves(self):
        dt = sheet(brood())[['SampleID', 'Sex', 'Species', 'BroodID', 'Parents']]
        ids, _ = C.resolve_members(dt, family_spec())
        self.assertEqual(len(ids), 26)


class TestParsing(unittest.TestCase):
    def test_sex_is_strict(self):
        self.assertEqual([C.parse_sex(x) for x in ('M', 'f', ' F ', 'U', 'F?', 'M?', '', None, float('nan'))], [1, 0, 0, None, None, None, None, None, None])

    def test_parents_follow_sexloci(self):
        self.assertEqual(C.parse_parents('S,D'), ('ok', 'S', 'D')); self.assertEqual(C.parse_parents(''), ('none', None, None)); self.assertEqual(C.parse_parents('S,?'), ('partial', 'S', None))
        self.assertEqual(C.parse_parents(',D'), ('partial', None, 'D')); self.assertEqual(C.parse_parents('S'), ('malformed', None, None)); self.assertEqual(C.parse_parents('S', known={'S'}), ('partial_single', 'S', None))
        self.assertEqual(C.parse_parents('a,b,c'), ('malformed', None, None))


if __name__ == '__main__':
    unittest.main()
