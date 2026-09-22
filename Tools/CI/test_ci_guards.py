"""Portable tests of CI guards. These do not execute Apple or biological workloads."""
from __future__ import annotations
import copy
import unittest
from verify_apple_toolchain import swift_version, validate
from verify_swift_suite_log import validate_log


def good_record():
    def result(output):
        return {"returncode": 0, "output": output}
    return {
        "system": "Darwin", "machine": "arm64", "developer_directory_exists": True,
        "commands": {
            "swift": result("Apple Swift version 6.3 (swiftlang-6.3.0.100 clang-1700.1)"),
            "xcrun_swift": result("Apple Swift version 6.3 (swiftlang-6.3.0.100 clang-1700.1)"),
            "xcode": result("Xcode 26.6\nBuild version 17F113"),
            "metal": result("/Applications/Xcode_26.6.app/Contents/Developer/Toolchains/Metal.xctoolchain/usr/bin/metal\n"),
        },
    }


GOOD_LOG = '''Test Suite 'Selected tests' passed at 2026-09-22 00:00:00.
Executed 0 tests, with 0 failures (0 unexpected) in 0.000 seconds
◇ Suite "Cell response learning" started.
✔ Test "rejects a failed update without replacing the model" passed after 0.2 seconds.
✔ Suite "Cell response learning" passed after 0.3 seconds.
✔ Test run with 1 test passed after 0.3 seconds.
'''


class ToolchainTests(unittest.TestCase):
    def test_apple_version_without_patch(self):
        self.assertEqual(swift_version("Apple Swift version 6.3 (swiftlang-x)"), (6,3,0))

    def test_full_version(self):
        self.assertEqual(swift_version("Swift version 6.3.1 (swift-6.3.1-RELEASE)"), (6,3,1))

    def test_driver_is_not_compiler_version(self):
        with self.assertRaises(ValueError):
            swift_version("swift-driver version: 1.120.5")

    def test_good_build_prerequisites(self):
        self.assertEqual(validate(good_record()), [])

    def test_original_6_1_toolchain_rejected(self):
        record=good_record()
        for key in ("swift", "xcrun_swift"):
            record["commands"][key]["output"]="Apple Swift version 6.1.2"
        self.assertEqual(len(validate(record)), 2)

    def test_6_2_also_rejected(self):
        record=good_record()
        for key in ("swift", "xcrun_swift"):
            record["commands"][key]["output"]="Swift version 6.2.9"
        self.assertEqual(len(validate(record)), 2)

    def test_numeric_not_decimal_version_order(self):
        record=good_record()
        for key in ("swift", "xcrun_swift"):
            record["commands"][key]["output"]="Swift version 6.10.1"
        self.assertEqual(validate(record), [])

    def test_no_linux_fallback(self):
        record=good_record(); record["system"]="Linux"
        self.assertTrue(validate(record))

    def test_arm_host_required(self):
        record=good_record(); record["machine"]="x86_64"
        self.assertTrue(validate(record))

    def test_explicit_xcode_required(self):
        record=good_record(); record["developer_directory_exists"]=False
        self.assertTrue(validate(record))

    def test_path_compiler_must_agree(self):
        record=good_record(); record["commands"]["swift"]["output"]="Swift version 6.4"
        self.assertTrue(any("differ" in e for e in validate(record)))

    def test_missing_each_tool_is_rejected(self):
        for key in good_record()["commands"]:
            with self.subTest(key=key):
                record=good_record(); record["commands"][key]={"returncode":127, "output":"unavailable"}
                self.assertTrue(validate(record))

    def test_malformed_version_is_rejected(self):
        record=good_record(); record["commands"]["swift"]["output"]="unrecognized output"
        self.assertTrue(validate(record))

    def test_empty_metal_resolution_rejected(self):
        record=good_record(); record["commands"]["metal"]["output"]="\n"
        self.assertTrue(validate(record))

    def test_validator_does_not_modify_its_input(self):
        record=good_record(); original=copy.deepcopy(record)
        validate(record)
        self.assertEqual(record, original)


class ExecutedSuiteTests(unittest.TestCase):
    def test_successful_swift_testing_not_empty_xctest(self):
        self.assertEqual(validate_log(GOOD_LOG, "Cell response learning"), 1)

    def test_no_matched_test_rejected(self):
        with self.assertRaises(ValueError):
            validate_log("Test run with 0 tests passed after 0.001 seconds", "Cell response learning")

    def test_other_suite_not_sufficient(self):
        with self.assertRaises(ValueError):
            validate_log(GOOD_LOG, "Other suite")

    def test_zero_even_with_suite_success_rejected(self):
        with self.assertRaises(ValueError):
            validate_log(GOOD_LOG.replace("with 1 test", "with 0 tests"), "Cell response learning")

    def test_partial_log_rejected(self):
        with self.assertRaises(ValueError):
            validate_log(GOOD_LOG.split("✔ Test run")[0], "Cell response learning")

    def test_skipped_test_rejected(self):
        with self.assertRaises(ValueError):
            validate_log(GOOD_LOG + '◇ Test "training" skipped after 0 seconds.\n', "Cell response learning")

    def test_skipped_reason_rejected(self):
        with self.assertRaises(ValueError):
            validate_log(GOOD_LOG + '◇ Test "training" skipped: no GPU.\n', "Cell response learning")

    def test_quoted_skip_reason_rejected(self):
        with self.assertRaises(ValueError):
            validate_log(GOOD_LOG + '◇ Test "training" skipped: "no GPU".\n', "Cell response learning")

    def test_quoted_failure_reason_rejected(self):
        with self.assertRaises(ValueError):
            validate_log(GOOD_LOG + '✘ Suite "Cell response learning" failed: "device lost".\n', "Cell response learning")

    def test_escaped_quote_in_skipped_test_name(self):
        with self.assertRaises(ValueError):
            validate_log(GOOD_LOG + r'◇ Test "the \"training\" case" skipped: "no GPU".' + "\n", "Cell response learning")

    def test_skipped_word_inside_passed_name_is_not_skip(self):
        text = GOOD_LOG.replace('rejects a failed update without replacing the model',
                                'handles a skipped batch without replacing the model')
        self.assertEqual(validate_log(text, "Cell response learning"), 1)

    def test_failed_test_rejected(self):
        with self.assertRaises(ValueError):
            validate_log(GOOD_LOG + '✘ Test "training" failed after 0 seconds.\n', "Cell response learning")

    def test_failed_name_is_not_failure_status(self):
        self.assertEqual(validate_log(GOOD_LOG, "Cell response learning"), 1)

    def test_repeated_log_is_rejected(self):
        with self.assertRaises(ValueError):
            validate_log(GOOD_LOG+GOOD_LOG, "Cell response learning")

    def test_ansi_and_suite_count(self):
        text=GOOD_LOG.replace('with 1 test passed', 'with 12 tests in 1 suite passed')
        self.assertEqual(validate_log("\x1b[32m"+text+"\x1b[0m", "Cell response learning"), 12)

    def test_empty_suite_rejected(self):
        with self.assertRaises(ValueError):
            validate_log(GOOD_LOG, "")

    def test_compiler_error_rejected(self):
        with self.assertRaises(ValueError):
            validate_log(GOOD_LOG+'error: load failed\n', "Cell response learning")


if __name__ == "__main__":
    unittest.main()
