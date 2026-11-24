"""클라우드 명단 전체 선택·VM별 선택·내부 worker와 기존 실행 호환을 검사한다."""
from contextlib import redirect_stdout
from functools import lru_cache
import copy
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from handoff import cli, local_queue, planner

ROOT = Path(planner.__file__).resolve().parents[1]


class CloudTargetTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name)
        self.config = {'data_root': str(self.base/'input'), 'output_root': str(self.base/'output'),
                       'state_root': str(self.base/'state'), 'output_layout': 'run_root',
                       'model_target_bytes': 1, 'memory_limit_bytes': 10**12,
                       'sample_cpus': 16, 'local_jobs': 1, 'cpu_pool_limit': 16,
                       'default_selection': 'manifests/리드별/medai/001/long_read.tsv',
                       'slurm': {'cpus': 16, 'concurrency': 1}}
        reader = lru_cache(maxsize=None)(planner.read_tsv)
        patcher = patch.object(planner, 'read_tsv', reader)
        patcher.start(); self.addCleanup(patcher.stop)
        patcher = patch.object(cli, 'ROOT', ROOT)
        patcher.start(); self.addCleanup(patcher.stop)

    def args(self, **changes):
        from types import SimpleNamespace
        fields = dict(target=None, selection=None, lab=None, batch=None,
                      sample=None, effective_config=None, command='dry-run')
        fields.update(changes)
        return SimpleNamespace(**fields)

    def selected(self, target='azure_poc', **extra):
        return cli.selection(self.args(target=target, **extra))

    def invoke(self, command, *extra, config=None):
        output = io.StringIO()
        with patch.object(sys, 'argv', ['handoff', command, *extra]), \
             patch.object(cli, 'load_config', return_value=copy.deepcopy(config or self.config)) as loader, \
             redirect_stdout(output):
            result = cli.main()
        return result, output.getvalue(), loader.call_args

    def snapshot(self, samples=None):
        rows = self.selected(sample=samples)
        c = dict(self.config, selection_target='azure_poc')
        if samples: c['selection_scope_samples'] = samples
        c = cli.batch_config(c, rows)
        effective = self.base/'effective.json'; effective.write_text(json.dumps(c))
        selection = self.base/'selection.tsv'; planner.write_tsv(selection, rows)
        return rows, c, effective, selection

    def test_four_targets_cover_all_samples_without_overlap(self):
        scopes = {name: self.selected(name) for name in cli.CLOUD_TARGETS}
        self.assertEqual({k: len(v) for k, v in scopes.items()},
                         {'azure_poc': 228, 'azure_main': 1257, 'aws_main': 1257, 'server': 470})
        ids = [r['sample_id'] for rows in scopes.values() for r in rows]
        self.assertEqual(len(ids), 3212)
        self.assertEqual(len(set(ids)), 3212)

    def test_short_pairs_are_split_between_poc_and_server(self):
        groups = {name: self.selected(name) for name in cli.CLOUD_TARGETS}
        counts = {name: sum(r['read_type'] == 'short_read' and r['data_group'] == 'wgs_rna_pair' for r in rows)
                  for name, rows in groups.items()}
        self.assertEqual(counts, {'azure_poc': 228, 'azure_main': 0, 'aws_main': 0, 'server': 188})
        for name, expected in [('azure_main', 0), ('aws_main', 0), ('server', 221)]:
            self.assertEqual(sum(r['data_group'] == 'wgs_rna_pair' for r in groups[name]), expected)

    def test_all_long_wgs_are_assigned_to_cloud(self):
        groups = {name: self.selected(name) for name in cli.CLOUD_TARGETS}
        self.assertEqual(sum(r['read_type'] == 'long_read' and r['data_group'] == 'wgs'
                             for name in ['azure_main', 'aws_main'] for r in groups[name]), 261)
        self.assertFalse(any(r['read_type'] == 'long_read' and r['data_group'] == 'wgs' for r in groups['server']))
        self.assertEqual(sum(r['read_type'] == 'long_read' and r['data_group'] == 'wgs_rna_pair'
                             for r in groups['server']), 33)

    def test_target_automatically_uses_its_profile_and_ignores_old_default_selection(self):
        result, output, call = self.invoke('dry-run', '--target', 'azure_poc')
        self.assertEqual(result, 0)
        payload = json.loads(output)
        self.assertEqual(payload['selected'], 228)
        self.assertEqual(payload['completion_policy'], 'all_selected')
        self.assertEqual(payload['config']['batch_key'], 'cloud_azure_poc')
        self.assertEqual(Path(call.args[0]), ROOT/'manifests/클라우드_분담/설정/azure_poc_pileup.json')
        self.assertFalse((self.base/'state').exists())

    def test_single_sample_does_not_need_a_lab(self):
        rows = self.selected(sample=['HG00308'])
        self.assertEqual([r['sample_id'] for r in rows], ['HG00308'])
        result, output, _ = self.invoke('dry-run', '--target', 'azure_poc', '--sample', 'HG00308')
        self.assertEqual(result, 0)
        self.assertEqual(json.loads(output)['config']['batch_key'], 'cloud_azure_poc_HG00308')

    def test_unknown_or_duplicate_samples_are_rejected(self):
        for samples in [['NOT_IN_TARGET'], ['HG00308', 'HG00308']]:
            with self.subTest(samples=samples), self.assertRaises(ValueError):
                self.selected(sample=samples)

    def test_legacy_selectors_cannot_be_mixed_with_target(self):
        for extra in [dict(lab='medai'), dict(batch='002'), dict(selection='anything.tsv')]:
            with self.subTest(extra=extra), self.assertRaises(ValueError):
                self.selected(**extra)

    def test_full_target_and_individual_jobs_have_distinct_states(self):
        full = cli.batch_config(dict(self.config, selection_target='azure_poc'), self.selected())
        keys = []
        for sample in ['HG00308', 'HG00260']:
            config = dict(self.config, selection_target='azure_poc', selection_scope_samples=[sample])
            keys.append(cli.batch_config(config, self.selected(sample=[sample]))['state_root'])
        self.assertEqual(len({full['state_root'], *keys}), 3)
        self.assertTrue(full['state_root'].endswith('cloud_azure_poc'))

    def test_internal_snapshot_accepts_full_target_above_old_batch_limit(self):
        rows, _, effective, selection = self.snapshot()
        got = cli.selection(self.args(effective_config=str(effective), selection=str(selection), command='worker'))
        self.assertEqual(got, rows)

    def test_internal_snapshot_order_change_is_rejected(self):
        rows, _, effective, selection = self.snapshot()
        planner.write_tsv(selection, list(reversed(rows)))
        with self.assertRaises(ValueError):
            cli.selection(self.args(effective_config=str(effective), selection=str(selection)))

    def test_saved_target_requires_snapshot_and_matching_target(self):
        _, _, effective, selection = self.snapshot()
        with self.assertRaises(ValueError):
            cli.selection(self.args(effective_config=str(effective)))
        with self.assertRaises(ValueError):
            cli.selection(self.args(target='aws_main', effective_config=str(effective), selection=str(selection)))

    def test_internal_subset_keeps_its_single_sample_scope(self):
        rows, _, effective, selection = self.snapshot(['HG00308'])
        self.assertEqual(cli.selection(self.args(effective_config=str(effective), selection=str(selection))), rows)

    def test_parent_run_passes_entire_target_to_existing_queue(self):
        with patch.object(local_queue, 'run', return_value=0) as run:
            result, _, _ = self.invoke('run', '--target', 'azure_poc')
        self.assertEqual(result, 0)
        self.assertEqual(len(run.call_args.args[1]), 228)
        self.assertEqual(run.call_args.args[0]['selection_target'], 'azure_poc')
        self.assertFalse((self.base/'state').exists())

    def test_internal_worker_can_use_index_above_220(self):
        rows, c, effective, selection = self.snapshot()
        with patch.object(cli, 'worker', return_value=0) as worker:
            result, _, _ = self.invoke('worker', '--effective-config', str(effective),
                                       '--selection', str(selection), '--task-index', '227', config=c)
        self.assertEqual(result, 0)
        self.assertEqual(worker.call_args.args[1], rows[227])

    def test_internal_download_worker_keeps_cloud_selection(self):
        from handoff import download_only
        rows, c, effective, selection = self.snapshot()
        with patch.object(download_only, 'worker', return_value=0) as worker:
            result, _, _ = self.invoke('download-worker', '--effective-config', str(effective),
                                       '--selection', str(selection), '--task-index', '227', config=c)
        self.assertEqual(result, 0)
        self.assertEqual(worker.call_args.args[1], rows[227])

    def test_cloud_claim_and_queue_are_not_cut_off_by_old_byte_goal(self):
        state = self.base/'state'; (state/'samples').mkdir(parents=True)
        (state/'samples/S.json').write_text(json.dumps({'sample': 'S', 'status': 'success', 'model_bytes': 2}))
        self.assertFalse(cli.claim(self.config, 'S2'))
        self.assertTrue(cli.claim(dict(self.config, selection_target='azure_poc'), 'S2'))
        adapter = object.__new__(local_queue.Adapter)
        adapter.config = dict(self.config)
        adapter.rows = [{'sample_id': 'S'}]
        adapter.verified = lambda row: True
        adapter.sample_state = lambda sample: state
        adapter.completed_sizes = {'S': 2}
        self.assertTrue(adapter.target_reached())
        adapter.config['selection_target'] = 'azure_poc'
        self.assertFalse(adapter.target_reached())

    def test_legacy_selection_and_state_key_remain_unchanged(self):
        row = self.selected(sample=['HG00308'])[0]
        selected = cli.selection(self.args(lab=row['assigned_lab'].lower(), sample=[row['sample_id']]))
        config = cli.batch_config(self.config, selected)
        self.assertEqual(config['batch_key'], row['assigned_lab'].lower()+'_'+row['batch_id'])
        with self.assertRaises(ValueError):
            cli.selection(self.args(selection=str(ROOT/'manifests/클라우드_분담/azure_poc.tsv')))

    def test_resume_samples_do_not_shrink_full_target_scope(self):
        self.assertEqual(len(self.selected(sample=['HG00308'], command='resume')), 228)

    def test_old_assignment_plan_cannot_be_invoked_as_cloud_preview(self):
        with patch.object(planner, 'build') as build:
            with self.assertRaises(ValueError): self.invoke('plan', '--target', 'azure_poc')
        build.assert_not_called()


if __name__ == '__main__':
    unittest.main()
