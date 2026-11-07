"""unit_scripts/sampledb.py with fixtures that reproduce the cases measured on the real export: a curated Exclude that conflicts with Blumer's clade, a different sex, species named differently by SRA and Blumer,
a Blumer-project sample that is not in Table S1, a sample with two SRA organisms, blank trailing columns, new samples with blank fields."""
import os, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'unit_scripts'))
import numpy as np
import pandas as pd
import sampledb as S

NAN = float('nan')
COLS = ['SampleID', 'Sex', 'Species', 'DoB', 'BroodID', 'Parents', 'Ecogroup', 'Subgroup', 'LabReared', 'Inversion10', 'Category', 'gsdf number', 'TE1', 'TE2', 'Unnamed: 14', 'Unnamed: 15']


def sheet():
    rows = [
        ['SAMEA1', 'M', 'Cynotilapia afra', '', '', '', 'Mbuna', 'Rocky', 0, 0.0, 'x', 2.0, '0/0', '0/0', NAN, NAN],                     # existing, in S1, SRA species differs from Blumer's
        ['SAMEA2', NAN, "Maylandia sp. 'pearly'", '', '', '', NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN],                      # new row of the campaign: Sex, Ecogroup, LabReared blank
        ['SAMEA3', 'U', 'Tropheops sp. olive', '', '', '', 'Exclude', NAN, 0, 2.0, NAN, NAN, NAN, NAN, NAN, NAN],                       # curated Exclude, Blumer says Shallow
        ['SAMEA4', 'F', 'Aulonocara minutus', '', '', '', 'Deep_Benthic', NAN, 0, 0.0, NAN, 3.0, '1/1', '0/1', NAN, NAN],                # sheet says F, Blumer says M
        ['SAMEA5', NAN, 'Labeotropheus fuelleborni', '', '', '', NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN],                     # Malinsky sample, not in Blumer's S1
        ['LAB1', 'F', 'Maylandia zebra', 'd', 'b', 'p', 'Mbuna', NAN, 1, NAN, NAN, NAN, NAN, NAN, NAN, NAN],                            # lab fish, no runs anywhere
        ['SAMEA6', NAN, 'Copadichromis sp.', '', '', '', NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN],                             # in S1 but Blumer's species is NA
        ['SAMEA7', NAN, 'Rhamphochromis longiceps', '', '', '', NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN],                      # two different SRA organisms among its runs
    ]
    return pd.DataFrame(rows, columns=COLS)


S1_ROWS = [  # sequence_id, specimen_id, seq_depth, genus, species, sex, clade, location, sublocation, Biosample, Bioproject
    ['cichlid1', 'D1', '15', 'Cynotilapia', 'zebroides', 'M', 'Mbuna', 'L', 'S', 'SAMEA1', 'PRJ'],
    ['cichlid2', 'D2', '15', 'Maylandia', 'pearly', 'F', 'Mbuna', 'L', 'S', 'SAMEA2', 'PRJ'],
    ['cichlid3', 'D3', '15', 'Tropheops', 'olive', 'U', 'Shallow', 'L', 'S', 'SAMEA3', 'PRJ'],
    ['cichlid4', 'D4', '15', 'Aulonocara', 'minutus', 'M', 'Deep', 'L', 'S', 'SAMEA4', 'PRJ'],
    ['cichlid6', 'D6', '15', 'Copadichromis', 'NA', 'U', 'Utaka', 'L', 'S', 'SAMEA6', 'PRJ'],
    ['cichlid7', 'D7', '15', 'Rhamphochromis', 'greyback', 'F?', 'Rhampho', 'L', 'S', 'SAMEA7', 'PRJ'],
    ['cichlid8', 'D8', '15', 'Diplotaxodon', 'bigeye', 'M', 'Diplo', 'L', 'S', 'SAMEA8', 'PRJ'],       # not in the sheet
    ['cichlid9', 'NA', '15', 'Astatotilapia', 'calliptera', 'M', 'AstCal', '', '', 'NA', 'NA'],        # no BioSample
]
S7_ROWS = [['cichlid1', '0', '0', '0', '0', '0'], ['cichlid2', '2', '1', '0', '0', '0'], ['cichlid3', '0', '0', '2', '0', '1'], ['cichlid4', '0', '0', '0', '1', '0'], ['cichlid6', '1', '1', '1', '1', '1'],
           ['cichlid7', '0', '0', '0', '0', '0'], ['cichlid8', '0', '0', '0', '0', '0'], ['cichlid9', '0', '0', '0', '0', '0']]


def tables():
    d = tempfile.mkdtemp() + '/'
    pd.DataFrame(S1_ROWS, columns=['sequence_id', 'specimen_id', 'seq_depth', 'genus', 'species', 'sex', 'clade', 'location', 'sublocation', 'Biosample', 'Bioproject']).to_csv(d + 'S1.tsv', sep='\t', index=False)
    pd.DataFrame(S7_ROWS, columns=['sample_id', 'chr2', 'chr9', 'chr10', 'chr11', 'chr13']).to_csv(d + 'S7.tsv', sep='\t', index=False)
    return d + 'S1.tsv', d + 'S7.tsv'


def reads():
    return pd.DataFrame([['SAMEA1', 'ERR1'], ['SAMEA1', 'ERR1b'], ['SAMEA2', 'ERR2'], ['SAMEA3', 'ERR3'], ['SAMEA4', 'ERR4'], ['SAMEA5', 'ERR5'], ['SAMEA6', 'ERR6'], ['SAMEA7', 'ERR7a'], ['SAMEA7', 'ERR7b'], ['ORPHAN', 'ERR9']],
                        columns=['SampleID', 'RunID'])


def tabs():
    blumer = pd.DataFrame([['ERR1', 'Cynotilapia afra', 'in_database', 'SAMEA1', ''], ['ERR1b', 'Cynotilapia afra', 'in_database', 'SAMEA1', ''], ['ERR2', "Maylandia sp. 'pearly'", 'in_database', 'SAMEA2', ''],
                           ['ERR3', 'Tropheops sp. olive', 'in_database', 'SAMEA3', ''], ['ERR4', 'Aulonocara minutus', 'in_database', 'SAMEA4', ''], ['ERR6', 'Copadichromis sp.', 'in_database', 'SAMEA6', ''],
                           ['ERR7a', 'Rhamphochromis longiceps', 'in_database', 'SAMEA7', ''], ['ERR7b', 'Rhamphochromis sp. chilingali', 'in_database', 'SAMEA7', ''], ['ERR77', 'Xx yy', 'in_database', 'SAMEA77', '']],
                          columns=['Run', 'Organism', 'DownloadStatus', 'BioSample', 'Library Name'])
    mal = pd.DataFrame([['ERR5', 'Labeotropheus fuelleborni', 'in_database', 'SAMEA5', '']], columns=['Run', 'Organism', 'DownloadStatus', 'BioSample', 'Library Name'])
    return {'Blumer': blumer, 'Malinsky': mal}


class WS:
    def __init__(self, fm): self.fm = fm
    row_count = property(lambda self: self.fm.grid[0]); col_count = property(lambda self: self.fm.grid[1])
    def resize(self, rows, cols): self.fm.grid = (rows, cols)


class FM:
    def __init__(self, corrupt=False, grid=(10, 16)):
        self.remote = sheet(); self.localTempDir = tempfile.mkdtemp() + '/'; self.s_ID, self.gc, self.grid, self.corrupt = 'S', self, grid, corrupt
        self.reads_n = self.writes = 0
    def open_by_key(self, k): return self
    def worksheet(self, name): return WS(self)
    def readSampleDatabase(self): self.reads_n += 1; self.sample_dt = self.remote.copy(); self.reads_dt = reads()
    def readDownloadTab(self, t): return tabs()[t]
    def writer(self, fm, df, log):
        assert df.shape[1] <= self.grid[1], 'the grid would have to be enlarged first'
        self.writes += 1; self.remote = df.copy()
        if self.corrupt:
            self.remote.loc[2, 'Species'] = 'DAMAGED'


def run(argv, fm=None):
    fm = fm or FM(); out = []
    res = S.main(argv, fm=fm, writer=fm.writer, out=out.append)
    return fm, res, '\n'.join(out)


class TestTheTablesLiveInTheRepository(unittest.TestCase):
    """The tool must run from a checkout with no flags: Blumer's tables are part of the repository, and they must be what the mapping and the join expect."""
    def test_the_default_tables_exist_and_join(self):
        self.assertTrue(os.path.exists(S.DEFAULT_S1) and os.path.exists(S.DEFAULT_S7), 'data/blumer/Table_S1.tsv and Table_S7.tsv must be in the repository')
        b, problems = S.read_blumer(S.DEFAULT_S1, S.DEFAULT_S7)
        self.assertEqual(problems, []); self.assertEqual(len(b), 1375)

    def test_every_clade_maps_to_an_ecogroup_and_the_inversion_codes_are_0_1_2(self):
        b, _ = S.read_blumer(S.DEFAULT_S1, S.DEFAULT_S7)
        self.assertEqual(sorted(set(b.clade) - set(S.ECO_MAP)), [])                                    # a new clade in the table would otherwise be skipped silently
        for c in S.INVERSIONS:
            self.assertLessEqual(set(b[c]) - {''}, {'0', '1', '2'}, c)
        self.assertTrue(b.sequence_id.is_unique and (b.Biosample[~b.Biosample.isin(['', 'NA'])]).is_unique)

    def test_commands_run_with_no_table_flags(self):
        fm = FM(); out = []
        S.main(['blumer', '--out-dir', tempfile.mkdtemp() + '/'], fm=fm, writer=fm.writer, out=out.append)
        self.assertIn('samples in Table S1: 1375, matched to the SampleDatabase: 0', '\n'.join(out))        # none of the fixture's sample IDs is a real Blumer BioSample: the default tables were read
        out = []; S.main(['check'], fm=FM(), out=out.append); self.assertIn('Blumer samples (Table S1) that are not in the SampleDatabase: 1375 of 1375', '\n'.join(out))
        out = []; S.main(['check', '--no-blumer'], fm=FM(), out=out.append); self.assertNotIn('Blumer samples', '\n'.join(out))

    def test_a_missing_table_says_where_it_should_be(self):
        with self.assertRaises(SystemExit) as cm:
            S.main(['blumer', '--s1', '/nonexistent/S1.tsv'], fm=FM(), out=lambda *a: None)
        self.assertIn('data/blumer', str(cm.exception))


class TestPlan(unittest.TestCase):
    def setUp(self):
        s1, s7 = tables(); self.b, self.problems = S.read_blumer(s1, s7)
        self.sra, self.multi = S.sra_species(reads(), tabs())
        self.new, self.changes, self.summary = S.plan(sheet(), self.b, self.sra)
        self.row = lambda sid: self.new[self.new.SampleID == sid].iloc[0]

    def test_the_two_blumer_tables_join_and_problems_are_reported(self):
        self.assertEqual(self.problems, []); self.assertEqual(len(self.b), 8)
        s1, s7 = tables(); pd.read_csv(s7, sep='\t', dtype=str).iloc[:-1].to_csv(s7, sep='\t', index=False)
        _, problems = S.read_blumer(s1, s7); self.assertTrue(any('not in S7' in p for p in problems))

    def test_species_is_overwritten_with_blumers_and_the_sra_name_is_kept_elsewhere(self):
        self.assertEqual((self.row('SAMEA1').Species, self.row('SAMEA2').Species, self.row('SAMEA4').Species), ('Cynotilapia zebroides', 'Maylandia pearly', 'Aulonocara minutus'))
        self.assertEqual((self.row('SAMEA1').SRASpecies, self.row('SAMEA2').SRASpecies), ('Cynotilapia afra', "Maylandia sp. 'pearly'"))
        self.assertEqual(self.row('SAMEA5').Species, 'Labeotropheus fuelleborni')                                                  # not in S1: unchanged
        self.assertEqual(self.row('SAMEA6').Species, 'Copadichromis sp.')                                                           # Blumer's species is NA: not overwritten
        self.assertTrue(any(c['SampleID'] == 'SAMEA6' and c['action'].startswith('skipped (Blumer species missing') for c in self.changes))
        self.assertIn((('Cynotilapia afra', 'Cynotilapia zebroides'), 1), S.species_pairs(self.changes))

    def test_blanks_are_filled_and_differing_values_are_kept_and_listed(self):
        r = self.row('SAMEA2')
        self.assertEqual((r.Sex, r.Ecogroup, r.LabReared, r.BlumerSequenceID, r.BlumerSpecimenID), ('F', 'Mbuna', 0, 'cichlid2', 'D2'))
        self.assertEqual((r.Inversion2, r.Inversion9, r.Inversion10, r.Inversion11, r.Inversion13), (2.0, 1.0, 0.0, 0.0, 0.0))
        self.assertEqual(self.row('SAMEA3').Ecogroup, 'Exclude'); self.assertEqual(self.row('SAMEA4').Sex, 'F')               # curated / different: not overwritten
        kept = {(c['SampleID'], c['column']): (c['old'], c['new']) for c in self.changes if c['action'] == 'kept (differs)'}
        self.assertEqual(kept[('SAMEA3', 'Ecogroup')], ('Exclude', 'Shallow_Benthic')); self.assertEqual(kept[('SAMEA4', 'Sex')], ('F', 'M'))
        self.assertNotIn(('SAMEA1', 'Sex'), kept)                                                                                    # identical values are not changes at all

    def test_inversion10_that_already_agrees_is_untouched_and_a_conflict_is_kept(self):
        self.assertEqual(self.row('SAMEA3').Inversion10, 2.0)                                                                        # sheet 2, Blumer chr10 2: same
        sh = sheet(); sh.loc[sh.SampleID == 'SAMEA4', 'Inversion10'] = 1.0                                                          # sheet 1, Blumer chr10 0
        new, changes, _ = S.plan(sh, self.b, self.sra)
        self.assertEqual(new[new.SampleID == 'SAMEA4'].iloc[0].Inversion10, 1.0); self.assertIn(('SAMEA4', 'Inversion10', 'kept (differs)'), [(c['SampleID'], c['column'], c['action']) for c in changes])

    def test_sra_species_covers_every_sample_with_runs_and_joins_several_organisms(self):
        self.assertEqual(self.row('SAMEA5').SRASpecies, 'Labeotropheus fuelleborni')                                                # a Malinsky sample, not in Blumer's tables
        self.assertEqual(self.row('SAMEA7').SRASpecies, 'Rhamphochromis longiceps | Rhamphochromis sp. chilingali'); self.assertEqual(self.multi, ['SAMEA7'])
        self.assertEqual(self.row('LAB1').SRASpecies, ''); self.assertNotIn('SAMEA77', self.sra)                                   # no runs / a run that is not in DNAReadsDatabase

    def test_columns_untouched_values_and_order(self):
        self.assertEqual(list(self.new.columns), COLS[:14] + ['SRASpecies', 'BlumerSequenceID', 'BlumerSpecimenID', 'Inversion2', 'Inversion9', 'Inversion11', 'Inversion13'] + COLS[14:])
        old = sheet(); self.assertEqual(list(self.new.SampleID), list(old.SampleID))
        for c in ('gsdf number', 'TE1', 'TE2', 'DoB', 'Subgroup', 'Category'):                                                      # nothing outside the planned columns moves
            self.assertTrue(all(S.same_value(a, b) for a, b in zip(old[c], self.new[c])), c)
        self.assertEqual((self.row('LAB1').Sex, self.row('LAB1').LabReared, self.row('LAB1').Species), ('F', 1, 'Maylandia zebra'))

    def test_the_summary_counts(self):
        c = self.summary['counts']
        self.assertEqual(c[('Species', 'overwritten')], 4)                                                                          # SAMEA1, 2, 3, 7; SAMEA4 is spelled the same, SAMEA6 has no Blumer species
        self.assertEqual(c[('Ecogroup', 'kept')], 1); self.assertEqual(c[('Sex', 'kept')], 1)
        self.assertEqual(self.summary['matched'], {'SAMEA1', 'SAMEA2', 'SAMEA3', 'SAMEA4', 'SAMEA6', 'SAMEA7'})


class TestCommand(unittest.TestCase):
    def test_report_only_changes_nothing_and_writes_the_reports(self):
        s1, s7 = tables(); od = tempfile.mkdtemp() + '/'
        fm, res, text = run(['blumer', '--s1', s1, '--s7', s7, '--out-dir', od])
        self.assertEqual(fm.writes, 0); self.assertIn('report only: nothing written', text)
        for f in ('sampledb_changes.tsv', 'sampledb_species_pairs.tsv', 'sampledb_blumer_unmatched.tsv', 'sampledb_blumer_runs_not_in_S1.tsv', 'sampledb_multiple_sra_species.tsv'):
            self.assertTrue(os.path.exists(od + f), f)
        un = open(od + 'sampledb_blumer_unmatched.tsv').read().split('\n'); self.assertIn('SAMEA8', un); self.assertIn('NA', un)
        self.assertEqual(open(od + 'sampledb_blumer_runs_not_in_S1.tsv').read().strip().split('\n'), ['SampleID_with_Blumer_runs_not_in_S1'])                                   # every Blumer-project sample of the fixture is in S1
        self.assertIn('samples in Table S1: 8, matched to the SampleDatabase: 6, not in the sheet: 2', text)
        pairs = open(od + 'sampledb_species_pairs.tsv').read(); self.assertIn("Cynotilapia afra\tCynotilapia zebroides", pairs)

    def test_write_backs_up_enlarges_the_grid_writes_once_and_verifies_cell_by_cell(self):
        s1, s7 = tables(); fm = FM(grid=(10, 16))
        fm, res, text = run(['blumer', '--s1', s1, '--s7', s7, '--out-dir', tempfile.mkdtemp() + '/', '--write'], fm=fm)
        self.assertEqual(fm.writes, 1); self.assertEqual(fm.grid[1], 23); self.assertIn('grid enlarged', text); self.assertIn('verified cell by cell: 8 rows', text)       # 23 = 14 + 7 new + 2 empty
        self.assertEqual(fm.remote[fm.remote.SampleID == 'SAMEA1'].iloc[0].Species, 'Cynotilapia zebroides')
        b = [d for d in os.listdir(fm.localTempDir) if d.startswith('database_backup_')]; self.assertEqual(len(b), 1)
        saved = pd.read_csv(fm.localTempDir + b[0] + '/SampleDatabase.csv'); self.assertEqual(saved[saved.SampleID == 'SAMEA1'].iloc[0].Species, 'Cynotilapia afra')      # the backup holds the old values

    def test_a_sheet_that_does_not_hold_what_was_written_is_reported_with_the_backup(self):
        s1, s7 = tables()
        with self.assertRaises(SystemExit) as cm:
            run(['blumer', '--s1', s1, '--s7', s7, '--out-dir', tempfile.mkdtemp() + '/', '--write'], fm=FM(corrupt=True))
        self.assertIn('Species: 1 cells differ', str(cm.exception)); self.assertIn('restore from', str(cm.exception))

    def test_readback_comparison_ignores_empty_versus_nan_and_2_versus_2_0(self):
        a = pd.DataFrame(dict(SampleID=['A', 'B'], X=[2, NAN], Y=['', 'v'])); b = pd.DataFrame(dict(SampleID=['A', 'B'], X=[2.0, ''], Y=[NAN, 'v']))
        self.assertEqual(S.readback_problems(a, b), []); b.loc[1, 'Y'] = 'w'; self.assertEqual(len(S.readback_problems(a, b)), 1)
        self.assertEqual(S.readback_problems(a, b.iloc[:1]), ['row count: wrote 2, found 1'])


class TestChecks(unittest.TestCase):
    def test_check_reports_the_problems_of_the_fixture(self):
        s1, s7 = tables(); out = []
        sh = sheet(); sh = pd.concat([sh, sh.iloc[[0]]], ignore_index=True); sh.loc[len(sh) - 2, 'Sex'] = 'male'; sh.loc[1, 'LabReared'] = 7
        fm = FM(); fm.remote = sh
        n = S.main(['check', '--s1', s1, '--s7', s7], fm=fm, out=out.append); text = '\n'.join(out)
        self.assertRegex(text, r'ATTENTION\s+duplicate SampleIDs: 1'); self.assertRegex(text, r'ATTENTION\s+samples with reads but no SampleDatabase row: 1\s+e.g. \[\'ORPHAN\'\]')
        self.assertRegex(text, r'ATTENTION\s+Sex values'); self.assertRegex(text, r'ATTENTION\s+LabReared values')
        self.assertRegex(text, r'ATTENTION\s+DownloadDatabase Blumer: 9 runs say in_database, 1 of them are not in DNAReadsDatabase')
        self.assertRegex(text, r'ATTENTION\s+Blumer samples \(Table S1\) that are not in the SampleDatabase: 2 of 8')
        self.assertIn('matched samples: 6; Species differs from Blumer for', text); self.assertGreater(n, 4)

    def test_a_clean_sheet_has_no_attention_lines_for_the_structure(self):
        out = []; fm = FM()
        S.main(['check'], fm=fm, out=out.append)
        self.assertNotRegex('\n'.join(out), r'ATTENTION\s+duplicate'); self.assertNotRegex('\n'.join(out), r'ATTENTION\s+Sex values')

    def test_species_written_in_several_ways_are_grouped(self):
        self.assertEqual(S.spelling_key("Tropheops sp. 'yellow gular'"), S.spelling_key('Tropheops yellow-gular')); self.assertNotEqual(S.spelling_key('Maylandia zebra'), S.spelling_key('Maylandia callainos'))
        sh = sheet(); sh.loc[0, 'Species'] = "Maylandia sp. 'pearly'"; sh.loc[1, 'Species'] = 'Maylandia pearly'; fm = FM(); fm.remote = sh; out = []
        S.main(['check'], fm=fm, out=out.append); self.assertRegex('\n'.join(out), r"ATTENTION\s+species names: \d+ distinct; 1 taxa written in more than one way")


if __name__ == '__main__':
    unittest.main()
