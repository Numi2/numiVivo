"""Corrupt continuation boundaries deliberately; fixtures are not MD evidence."""
import copy
import math
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from importlib import metadata
from ensemble_campaign import read, write, digest, analysis_environment, campaign_grid
from continue_ensemble import CONTRACT, SOURCES, continuation_schedule, verify_chunk, verify_files, file_manifest, load_parent, analyze
from test_ensemble_campaign import fixture

TOOLS = Path(__file__).parent


def chunk_fixture(steps=100):
    request, original, policy = fixture()
    start = copy.deepcopy(original["dynamics"]["finalCheckpoint"])
    start["numericalContract"] = CONTRACT
    final = copy.deepcopy(start)
    final["acceptedStep"] += steps
    final["timePS"] += steps * request["configuration"]["timeStepPS"]
    samples = []
    for offset in range(50, steps+1, 50):
        obs = copy.deepcopy(original["dynamics"]["end"])
        obs.update(stepIndex=start["acceptedStep"]+offset, timePS=start["timePS"]+offset*.001)
        samples.append(dict(stepIndex=obs["stepIndex"], timePS=obs["timePS"],
            positionsNM=copy.deepcopy(final["positionsNM"]), periodicCell=copy.deepcopy(final["periodicCell"]), observables=obs))
    report = dict(schema="numivivo.org/md-run-report/v1", requestedSteps=steps, committedSteps=steps,
        startStep=start["acceptedStep"], endStep=final["acceptedStep"], startTimePS=start["timePS"], endTimePS=final["timePS"],
        systemFingerprint=start["systemFingerprint"], configurationFingerprint=start["configurationFingerprint"],
        finalCheckpoint=final, samples=samples)
    return request, start, report, policy


class ContinuationTests(unittest.TestCase):
    def test_fixed_six_replica_schedule_preserves_acceptance_policy(self):
        base = read(TOOLS/"ensemble_policy.json")
        original = copy.deepcopy(base)
        followup = read(TOOLS/"ensemble_continuation_policy.json")
        chunks, every = continuation_schedule(base, followup)
        self.assertEqual(chunks, [20000]*16)
        self.assertEqual(every, 100)
        self.assertEqual(base, original)
        self.assertEqual(len([r for r in campaign_grid(base) if r[0] == followup["timeStepPS"]]), 6)

    def test_unbounded_or_incomplete_schedules_refused(self):
        base = read(TOOLS/"ensemble_policy.json")
        for key, value in [("totalDurationPS", 50), ("totalDurationPS", 210.025),
            ("totalDurationPS", 10000), ("chunkDurationPS", 100), ("chunkDurationPS", .025),
            ("bootstrapBlockPS", 7), ("bootstrapBlockPS", 100), ("bootstrapBlockPS", math.nan),
            ("timeStepPS", .001), ("chunkDurationPS", True)]:
            policy = read(TOOLS/"ensemble_continuation_policy.json")
            policy[key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                continuation_schedule(base, policy)
        policy = read(TOOLS/"ensemble_continuation_policy.json")
        policy["maximumStatisticalZ"] = 100
        with self.assertRaisesRegex(ValueError, "duration/block changes only"):
            continuation_schedule(base, policy)

    def test_zero_step_preserves_exact_words_and_random_counter(self):
        request, start, report, policy = chunk_fixture(0)
        self.assertEqual(verify_chunk(request, start, report, 0, 50, policy)["kinetic"], [])
        report["finalCheckpoint"]["velocitiesNMPerPS"] = [dict(x=2, y=0, z=0) for _ in start["positionsNM"]]
        with self.assertRaisesRegex(ValueError, "exact checkpoint"):
            verify_chunk(request, start, report, 0, 50, policy)

    def test_continuation_excludes_duplicate_boundary_without_new_equilibration(self):
        request, start, report, policy = chunk_fixture()
        result = verify_chunk(request, start, report, 100, 50, policy)
        self.assertEqual(result["kinetic"], [1.5, 1.5])
        self.assertTrue(result["endpointKinetic"]["passed"])

    def test_changed_clocks_missing_observations_and_state_identity_refused(self):
        mutations = [lambda r:r.update(startStep=0), lambda r:r.update(startTimePS=0),
            lambda r:r.update(committedSteps=99), lambda r:r.update(rejected={"committed":False}),
            lambda r:r["samples"].pop(), lambda r:r["samples"].insert(0, r["samples"][0]),
            lambda r:r["samples"][0].update(stepIndex=50000),
            lambda r:r["samples"][0].update(timePS=math.nan),
            lambda r:r["samples"][0]["observables"].update(configurationFingerprint="other"),
            lambda r:r["samples"][0]["observables"].update(degreesOfFreedom=9),
            lambda r:r["samples"][0]["observables"].update(totalEnergyKJPerMol=8),
            lambda r:r["samples"][0]["observables"].update(temperatureK=300),
            lambda r:r["samples"][0]["positionsNM"][0].update(x=math.nan),
            lambda r:r["samples"][-1]["positionsNM"][0].update(x=.001),
            lambda r:r["finalCheckpoint"].update(numericalContract="other"),
            lambda r:r["finalCheckpoint"]["positionCorrectionsNM"][0].update(x=2**-30),
            lambda r:r["finalCheckpoint"]["velocitiesNMPerPS"][0].update(x=10)]
        for i, mutate in enumerate(mutations):
            request, start, report, policy = chunk_fixture()
            mutate(report)
            with self.subTest(mutation=i), self.assertRaises((ValueError, AssertionError)):
                verify_chunk(request, start, report, 100, 50, policy)

    def test_matching_energy_accounting_does_not_hide_wrong_kinetic_energy(self):
        request, start, report, policy = chunk_fixture()
        for sample in report["samples"]:
            sample["observables"].update(kineticEnergyKJPerMol=3, totalEnergyKJPerMol=3,
                temperatureK=6/(6*policy["boltzmannKJPerMolK"]))
        with self.assertRaisesRegex(ValueError, "endpoint kinetic"):
            verify_chunk(request, start, report, 100, 50, policy)

    def test_artifact_changes_and_escaping_paths_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root/"checkpoint.json", dict(acceptedStep=10))
            files = {"checkpoint.json":digest(root/"checkpoint.json")}
            verify_files(root, files)
            (root/"checkpoint.json").write_text('{"acceptedStep":0}')
            with self.assertRaises(ValueError):verify_files(root, files)
            with self.assertRaises(ValueError):verify_files(root, {"../checkpoint.json":"invalid"})

    def test_parent_qualification_is_required_and_failure_is_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = read(TOOLS/"ensemble_policy.json")
            write(root/"policy.json", base)
            write(root/"binary-manifest.json", {})
            write(root/"parent-manifest.json", {})
            identities = dict(analysisEnvironment=analysis_environment(), policySHA256=digest(root/"policy.json"),
                binaryManifestSHA256=digest(root/"binary-manifest.json"), parentManifestSHA256=digest(root/"parent-manifest.json"),
                tools={n:digest(TOOLS/n) for n in ("ensemble_campaign.py", "ensemble_statistics.py", "run_campaign.py")})
            campaign = dict(schema="numivivo.org/md-ensemble-campaign/v1", policy=base, identities=identities,
                runs=[dict(timeStepPS=d, temperatureK=t, seed=s) for d,t,s in campaign_grid(base)])
            write(root/"campaign.json", campaign)
            qualification = root/"qualification.json"
            with self.assertRaises(FileNotFoundError):load_parent(root, qualification)
            write(qualification, dict(campaignSHA256=digest(root/"campaign.json"), preparedOutcome="failed"))
            _, result = load_parent(root, qualification)
            self.assertEqual(result["preparedOutcome"], "failed")
            (root/"policy.json").write_text('{}')
            with self.assertRaises(ValueError):load_parent(root, qualification)

    def test_offline_analysis_retains_failed_replicas_and_rejects_altered_identity(self):
        for corruption in (None, "duplicate", "policy", "source", "parent"):
            with self.subTest(corruption=corruption), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                base = read(TOOLS/"ensemble_policy.json")
                policy = read(TOOLS/"ensemble_continuation_policy.json")
                write(root/"policy.json", policy)
                write(root/"parent-campaign.json", dict(policy=base))
                write(root/"parent-qualification.json", dict(preparedOutcome="inconclusive"))
                identities = dict(parentCampaignSHA256=digest(root/"parent-campaign.json"),
                    parentQualificationSHA256=digest(root/"parent-qualification.json"),
                    environment=dict(**analysis_environment(), openmm=metadata.version("openmm")),
                    tools={n:digest(TOOLS/n) for n in SOURCES})
                write(root/"identities.json", identities)
                snapshots = file_manifest(root)
                rows = []
                for dt, t, seed in campaign_grid(base):
                    if dt != policy["timeStepPS"]:continue
                    name = f"fixture-{t}-{seed}"
                    (root/name).mkdir()
                    rows.append(dict(identifier=base["cases"][0], directory=name, timeStepPS=dt,
                        temperatureK=t, seed=seed, files={}, chunks=[], outcome="failed"))
                campaign = dict(schema="numivivo.org/md-ensemble-continuation/v1", policy=policy,
                    parentPolicy=base, identities=identities, snapshots=snapshots, runs=rows,
                    parentPreparedOutcome="inconclusive", inputFailures=[dict(identifier="original-water", outcome="preparation-failed")])
                if corruption == "duplicate":rows.append(rows[0])
                elif corruption == "policy":policy["bootstrapBlockPS"] = 10
                elif corruption == "source":identities["tools"]["continue_ensemble.py"] = "changed"
                elif corruption == "parent":identities["parentQualificationSHA256"] = "changed"
                write(root/"campaign.json", campaign)
                args = SimpleNamespace(campaign=root, out=root/"qualification.json")
                if corruption:
                    with self.assertRaises(ValueError):analyze(args)
                    self.assertFalse(args.out.exists())
                else:
                    self.assertEqual(analyze(args), 1)
                    result = read(args.out)
                    self.assertEqual(result["preparedOutcome"], "failed")
                    self.assertEqual(len(result["validationErrors"]), 6)
                    self.assertEqual(result["inputFailures"], campaign["inputFailures"])
                    self.assertEqual(result["parentPreparedOutcome"], "inconclusive")


if __name__ == "__main__":
    unittest.main()
