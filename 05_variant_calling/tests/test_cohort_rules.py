"""The cohort check: every rule, by kind of cohort. A clean cohort must pass; each rule has a case that triggers it."""
import json, os, sys, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_fixtures import *
from variant_calling import rules as R

RULES = dict(R.DEFAULT_RULES)


def run(spec, rows, led=None, rules=None):
    dt = sheet(rows)
    from variant_calling.cohort import resolve_members
    members, _ = resolve_members(dt, spec)
    return R.check_cohort('t', spec, dt, led if led is not None else ledger(members), rules or RULES)


class Base(unittest.TestCase):
    def has(self, rep, level, code):
        self.assertIn(code, rep.codes(level), f'expected {level} {code}; got {[(f.level, f.code) for f in rep.findings]}')

    def lacks(self, rep, code):
        self.assertNotIn(code, rep.codes(), f'unexpected {code}: {rep.text()}')


class TestFamily(Base):
    def test_a_clean_family_passes_without_warnings(self):
        rep = run(family_spec(), brood())
        self.assertFalse(rep.blocked, rep.text()); self.assertEqual(rep.count('WARN'), 0, rep.text()); self.assertEqual(len(rep.members), 26)

    def test_an_invalid_sex_blocks_and_parents_are_included(self):
        rows = brood(); rows[0]['Sex'] = 'U'; rows[5]['Sex'] = 'F?'
        rep = run(family_spec(), rows); self.has(rep, 'BLOCK', 'sex_not_mf')
        self.assertIn('S_m', [s for f in rep.findings if f.code == 'sex_not_mf' for s in f.samples], 'a parent with sex U blocks too')

    def test_two_species_block_and_a_wrong_declared_species_blocks(self):
        rows = brood(); rows[3]['Species'] = 'Other sp.'
        self.has(run(family_spec(), rows), 'BLOCK', 'species_mismatch')
        self.has(run(family_spec(species='Not this one'), brood()), 'BLOCK', 'species_not_as_declared')
        rows = brood(); rows[2]['Species'] = ''
        self.has(run(family_spec(), rows), 'BLOCK', 'species_missing')

    def test_species_spelling_differences_in_case_and_spaces_do_not_block(self):
        rows = brood(); rows[3]['Species'] = YH.upper(); rows[4]['Species'] = '  ' + YH.replace(' ', '  ')
        self.assertFalse(run(family_spec(), rows).blocked)

    def test_a_listed_parent_missing_from_the_cohort_blocks_unlike_sexloci(self):
        rep = run(family_spec(include_parents=False), brood())
        self.has(rep, 'BLOCK', 'parent_not_in_cohort'); self.assertEqual(rep.codes('BLOCK').count('parent_not_in_cohort'), 2)

    def test_no_member_listing_parents_blocks(self):
        rows = brood(); [r.update(Parents='') for r in rows]
        self.has(run(family_spec(), rows), 'BLOCK', 'no_pedigree')

    def test_one_unknown_parent_and_a_malformed_field_warn(self):
        rows = brood(); rows[2]['Parents'] = 'S_m,'; rows[3]['Parents'] = 'S_m D_f'
        rep = run(family_spec(), rows)
        self.has(rep, 'WARN', 'parent_unknown'); self.has(rep, 'WARN', 'parents_malformed'); self.assertFalse(rep.blocked, rep.text())

    def test_swapped_parent_sexes_warn(self):
        rows = brood(); rows[0]['Sex'] = 'F'; rows[1]['Sex'] = 'M'
        self.has(run(family_spec(), rows), 'WARN', 'parent_sex_mismatch')

    def test_a_small_brood_warns_and_few_of_one_sex_block(self):
        rep = run(family_spec(), brood(n_m=5, n_f=5)); self.has(rep, 'WARN', 'small_family'); self.has(rep, 'WARN', 'low_per_sex'); self.assertFalse(rep.blocked)
        self.has(run(family_spec(), brood(n_m=2, n_f=14)), 'BLOCK', 'too_few_per_sex')

    def test_parents_are_not_counted_towards_the_per_sex_minimum(self):
        rep = run(family_spec(), brood(n_m=2, n_f=2)); self.has(rep, 'BLOCK', 'too_few_per_sex')
        self.assertIn('2 males and 2 females', rep.text())

    def test_coverage_warnings_by_role(self):
        led = ledger([r['SampleID'] for r in brood()], **{'B1_m0': 1.0, 'S_m': 6.0})
        rep = run(family_spec(), brood(), led)
        self.has(rep, 'WARN', 'low_coverage_offspring'); self.has(rep, 'WARN', 'low_coverage_parent'); self.assertFalse(rep.blocked)
        led = ledger([r['SampleID'] for r in brood()], **{'B1_m0': 3.0})
        self.lacks(run(family_spec(), brood(), led), 'low_coverage_offspring')      # offspring only warn below 2x

    def test_a_member_without_an_alignment_row_blocks(self):
        led = ledger([r['SampleID'] for r in brood() if r['SampleID'] != 'B1_f1'])
        rep = run(family_spec(), brood(), led); self.has(rep, 'BLOCK', 'no_alignment'); self.assertIn('B1_f1', rep.findings[[f.code for f in rep.findings].index('no_alignment')].samples)

    def test_a_different_genome_blocks_and_a_blank_genome_version_does_not(self):
        ids = [r['SampleID'] for r in brood()]
        led = ledger(ids); led.loc[led.SampleID == 'B1_m1', 'GenomeVersion'] = 'Mconophoros_GT1'
        self.has(run(family_spec(), brood(), led), 'BLOCK', 'genome_version_mismatch')
        led.loc[led.SampleID == 'B1_m1', 'GenomeVersion'] = ''
        self.lacks(run(family_spec(), brood(), led), 'genome_version_mismatch')

    def test_a_sample_not_in_the_sheet_blocks_with_a_likely_match(self):
        rep = run(family_spec(samples=['B1-m0'], include_parents=True), brood())
        self.has(rep, 'BLOCK', 'sample_not_in_table'); self.assertIn('B1_m0', rep.text())

    def test_duplicate_sheet_rows_block(self):
        rows = brood(); rows.append(dict(rows[5]))
        self.has(run(family_spec(), rows), 'BLOCK', 'duplicate_rows')

    def test_a_flagged_ecogroup_member_included_by_name_warns(self):
        rows = brood(); rows[6]['Ecogroup'] = 'Error'
        self.has(run(family_spec(), rows), 'WARN', 'flagged_ecogroup')

    def test_an_alignment_row_without_a_coverage_value_is_noted(self):
        led = ledger([r['SampleID'] for r in brood()], **{'B1_m1': None})
        rep = run(family_spec(), brood(), led); self.has(rep, 'INFO', 'no_coverage_value'); self.assertFalse(rep.blocked)

    def test_a_bad_definition_blocks_before_anything_is_read(self):
        rep = R.check_cohort('t', family_spec(kind='nope'), sheet(brood()), ledger([]), RULES)
        self.has(rep, 'BLOCK', 'bad_definition'); self.assertEqual(rep.members, [])

    def test_selectors_that_match_nothing_block(self):
        self.has(run(family_spec(broods=['NOPE']), brood()), 'BLOCK', 'no_members')

    def test_two_broods_are_judged_separately(self):
        rows = brood('B1') + brood('B2', sire='S2_m', dam='D2_f', n_m=4, n_f=4)
        rep = run(family_spec(broods=['B1', 'B2']), rows)
        small = [f for f in rep.findings if f.code == 'small_family']
        self.assertEqual(len(small), 1); self.assertIn('S2_m', small[0].message)


class TestSexSet(Base):
    def pool(self, n=12):
        return [fish(f'U{i}', 'M' if i % 2 else 'F', brood=f'FAM{i}') for i in range(n)]

    def test_a_clean_pool_passes(self):
        rep = run(sexset_spec(), self.pool()); self.assertFalse(rep.blocked, rep.text()); self.assertEqual(rep.count('WARN'), 0, rep.text())

    def test_parents_are_not_required_and_no_pedigree_is_fine(self):
        self.lacks(run(sexset_spec(), self.pool()), 'no_pedigree')

    def test_too_few_of_one_sex_blocks(self):
        rows = [fish(f'M{i}', 'M') for i in range(8)] + [fish('F0', 'F'), fish('F1', 'F')]
        self.has(run(sexset_spec(), rows), 'BLOCK', 'too_few_per_sex')

    def test_siblings_and_parent_with_offspring_are_flagged(self):
        rows = self.pool(12); rows[0]['BroodID'] = 'SAME'; rows[1]['BroodID'] = 'SAME'
        rep = run(sexset_spec(), rows); self.has(rep, 'WARN', 'related_members'); self.assertIn('U0/U1', rep.text())
        rows = self.pool(12); rows[4]['Parents'] = 'U1,U0'
        self.has(run(sexset_spec(), rows), 'WARN', 'related_members')

    def test_sex_confounded_with_a_batch_column_warns(self):
        rows = [fish(f'M{i}', 'M', project='PA', brood=f'A{i}') for i in range(6)] + [fish(f'F{i}', 'F', project='PB', brood=f'B{i}') for i in range(6)]
        rep = run(sexset_spec(), rows); self.has(rep, 'WARN', 'sex_confounded'); self.assertIn('Project', rep.text())
        mixed = [fish(f'M{i}', 'M', project='PA' if i < 3 else 'PB', brood=f'A{i}') for i in range(6)] + [fish(f'F{i}', 'F', project='PA' if i < 3 else 'PB', brood=f'B{i}') for i in range(6)]
        self.lacks(run(sexset_spec(), mixed), 'sex_confounded')

    def test_shallow_fish_are_listed_and_the_unusable_ones_counted(self):
        rows = self.pool(12); led = ledger([r['SampleID'] for r in rows], U1=4.0, U2=1.5)
        rep = run(sexset_spec(), rows, led); self.has(rep, 'WARN', 'low_coverage'); self.assertIn('1 below 3x, unusable', rep.text())

    def test_unsexed_fish_are_skipped_by_the_selector_not_blocked(self):
        rows = self.pool(12) + [fish('UNK', 'U')]
        rep = run(sexset_spec(), rows); self.assertFalse(rep.blocked); self.assertNotIn('UNK', rep.members); self.has(rep, 'INFO', 'selection')


class TestBcParents(Base):
    def cross(self):
        return [fish('YH_008_m', 'M', YH, eco='Deep_Benthic'), fish('MC-1R6R-f', 'F', 'Mchenga conophoros'), fish('F1_002', 'U', 'hybrid', parents='YH_008_m,MC-1R6R-f'), fish('F1_003', 'F', 'hybrid', parents='YH_008_m,MC-1R6R-f')]

    def test_two_species_unknown_sex_and_no_per_sex_minimum_pass(self):
        rep = run(bc_spec(['YH_008_m', 'MC-1R6R-f', 'F1_002', 'F1_003']), self.cross())
        self.assertFalse(rep.blocked, rep.text()); self.has(rep, 'WARN', 'sex_unknown'); self.has(rep, 'INFO', 'cohort_summary')
        for code in ('species_mismatch', 'too_few_per_sex', 'sex_not_mf', 'no_pedigree'):
            self.lacks(rep, code)

    def test_the_same_cohort_as_a_family_would_be_refused(self):
        rows = self.cross(); [r.update(BroodID='X') for r in rows]
        rep = run(family_spec(broods=['X']), rows); self.assertTrue(rep.blocked); self.has(rep, 'BLOCK', 'species_mismatch')

    def test_an_f1_whose_parent_is_not_in_the_cohort_blocks(self):
        rep = run(bc_spec(['YH_008_m', 'F1_002', 'F1_003']), self.cross())
        self.has(rep, 'BLOCK', 'pedigree_incomplete'); self.assertIn('MC-1R6R-f', rep.text())

    def test_founders_need_no_parents(self):
        self.assertFalse(run(bc_spec(['YH_008_m', 'MC-1R6R-f']), self.cross()).blocked)

    def test_shallow_members_block_and_moderate_ones_warn(self):
        ids = ['YH_008_m', 'MC-1R6R-f', 'F1_002', 'F1_003']
        rep = run(bc_spec(ids), self.cross(), ledger(ids, **{'F1_002': 2.0, 'F1_003': 7.0}))
        self.has(rep, 'BLOCK', 'coverage_too_low'); self.has(rep, 'WARN', 'low_coverage')
        self.assertIn('provisional', rep.text())

    def test_the_threshold_is_configurable(self):
        ids = ['YH_008_m', 'MC-1R6R-f']
        rules = dict(RULES, bc_parents_min_coverage=15.0)
        self.has(run(bc_spec(ids), self.cross(), ledger(ids, cov=10.0), rules), 'BLOCK', 'coverage_too_low')

    def test_a_member_without_an_alignment_blocks_here_too(self):
        self.has(run(bc_spec(['YH_008_m', 'MC-1R6R-f']), self.cross(), ledger(['YH_008_m'])), 'BLOCK', 'no_alignment')


class TestEmptySelectors(Base):
    """Found on the live sheet (2026-11-10): a cohort that selected adults by a sheet label gave exactly the members of its brood and the check said ok, because a selector that matches
    nothing is not an error. It is a WARN, so a cohort cannot silently come out smaller than meant."""
    def test_an_ecogroup_that_matches_nobody_warns(self):
        rows = brood() + [fish('A1', 'M', eco='Shallow_Benthic'), fish('A2', 'F', eco='Shallow_Benthic')]
        rep = run(sexset_spec(ecogroups=['Deep_Benthic']), rows)
        self.has(rep, 'WARN', 'selector_empty'); self.has(rep, 'BLOCK', 'no_members')
        msgs = ' '.join(f.message for f in rep.findings if f.code == 'selector_empty')
        self.assertIn('no sample in the sheet matches species', msgs); self.assertIn("ecogroup 'Deep_Benthic' appears nowhere in the sheet column Ecogroup", msgs)

    def test_an_ecogroup_that_exists_but_not_for_this_species_warns_once(self):
        rows = brood() + [fish('OTHER1', 'M', species='Other species', eco='Deep_Benthic'), fish('OTHER2', 'F', species='Other species', eco='Deep_Benthic')]
        rep = run(sexset_spec(ecogroups=['Deep_Benthic']), rows)
        msgs = [f.message for f in rep.findings if f.code == 'selector_empty']
        self.assertEqual(len(msgs), 1); self.assertIn('no sample in the sheet matches', msgs[0])

    def test_a_species_name_that_matches_nobody_blocks_and_explains(self):
        rep = run(sexset_spec(species='Aulonocara nowhere'), brood())
        self.has(rep, 'BLOCK', 'no_members'); self.has(rep, 'WARN', 'selector_empty')

    def test_a_typo_among_several_broods_warns_and_the_rest_still_select(self):
        rep = run(family_spec(broods=['B1', 'B1_TYPO']), brood())
        self.has(rep, 'WARN', 'selector_empty'); self.assertIn("brood 'B1_TYPO'", rep.text()); self.assertEqual(len(rep.members), 26)

    def test_selectors_that_match_do_not_warn(self):
        self.lacks(run(family_spec(), brood()), 'selector_empty')
        self.lacks(run(family_spec(species=YH), brood()), 'selector_empty')             # species next to broods only declares the species
        self.lacks(run(sexset_spec(), brood()), 'selector_empty')

    def test_an_empty_selection_both_blocks_and_explains(self):
        rep = run(family_spec(broods=['NOPE']), brood()); self.has(rep, 'BLOCK', 'no_members'); self.has(rep, 'WARN', 'selector_empty')


class TestBroodConsistency(Base):
    """A BroodID is one family of two parents: the children must all name the same sire and dam. (On the older sheet export 3 of 11 broods, the BC broods, had children with a blank Parents field.)"""
    def test_a_consistent_brood_has_no_such_warning(self):
        rep = run(family_spec(), brood())
        self.lacks(rep, 'brood_parents_inconsistent'); self.lacks(rep, 'brood_member_without_parents')

    def test_children_naming_different_parent_pairs_warn(self):
        rows = brood(); rows[4]['Parents'] = 'S_m,OTHER_f'; rows[5]['Parents'] = 'S_m,OTHER_f'
        rows.append(fish('OTHER_f', 'F'))
        rep = run(family_spec(), rows)
        self.has(rep, 'WARN', 'brood_parents_inconsistent'); self.assertIn('2 different parent pairs', rep.text()); self.assertIn('S_m x OTHER_f: 2', rep.text())

    def test_a_child_with_a_blank_parents_field_warns(self):
        rows = brood(); rows[6]['Parents'] = ''; rows[7]['Parents'] = ''
        rep = run(family_spec(), rows)
        self.has(rep, 'WARN', 'brood_member_without_parents'); self.assertIn('2 of 24 children', rep.text()); self.lacks(rep, 'brood_parents_inconsistent')

    def test_each_brood_is_judged_on_its_own(self):
        rows = brood('B1') + brood('B2', sire='S2_m', dam='D2_f'); rows[-1]['Parents'] = ''
        rep = run(family_spec(broods=['B1', 'B2']), rows)
        msgs = [f.message for f in rep.findings if f.code == 'brood_member_without_parents']
        self.assertEqual(len(msgs), 1); self.assertIn('brood B2', msgs[0])

    def test_it_does_not_apply_to_a_sex_set_or_a_cross(self):
        rows = [fish(f'U{i}', 'M' if i % 2 else 'F', brood='SAME') for i in range(12)]
        self.lacks(run(sexset_spec(), rows), 'brood_member_without_parents')


class TestRelatednessUnknown(Base):
    def test_lab_reared_fish_with_no_pedigree_are_noted(self):
        rows = [fish(f'U{i}', 'M' if i % 2 else 'F') for i in range(12)]
        rep = run(sexset_spec(), rows)
        self.has(rep, 'INFO', 'relatedness_unknown'); self.assertIn('12 of 12 analysed fish are lab-reared', rep.text()); self.assertFalse(rep.blocked); self.assertEqual(rep.count('WARN'), 0, rep.text())

    def test_wild_caught_fish_and_fish_with_a_brood_are_not_noted(self):
        rows = [fish(f'W{i}', 'M' if i % 2 else 'F', lab='0') for i in range(12)]
        self.lacks(run(sexset_spec(), rows), 'relatedness_unknown')
        rows = [fish(f'U{i}', 'M' if i % 2 else 'F', brood=f'FAM{i}') for i in range(12)]
        self.lacks(run(sexset_spec(), rows), 'relatedness_unknown')


class TestReport(Base):
    def test_text_dict_and_json(self):
        rep = run(family_spec(), brood(n_m=2, n_f=2))
        d = rep.to_dict(); json.dumps(d)
        self.assertTrue(d['blocked']); self.assertEqual(d['kind'], 'family'); self.assertEqual(d['n_members'], 6); self.assertGreaterEqual(d['n_block'], 1)
        t = rep.text(); self.assertTrue(t.splitlines()[0].startswith('cohort t (family, Mzebra_GT3_NCBI): 6 samples')); self.assertIn('RESULT: BLOCKED', t)
        self.assertIn('RESULT: ok to build', run(family_spec(), brood()).text())

    def test_blocks_are_listed_before_warnings(self):
        text = run(family_spec(), brood(n_m=2, n_f=2)).text()
        self.assertLess(text.index('[BLOCK]'), text.index('[WARN]'))

    def test_the_shipped_rules_file_matches_the_defaults(self):
        self.assertEqual(R.load_rules(), RULES)


if __name__ == '__main__':
    unittest.main()
