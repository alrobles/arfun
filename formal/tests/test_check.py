import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location("arfun_formal_check", Path(__file__).resolve().parents[1] / "check.py")
CHECK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECK)


class AxiomAuditTests(unittest.TestCase):
    def test_accepts_standard_axioms(self):
        name = "ArfunFormal.Geometry.valid"
        result = CHECK.parse_axioms(
            f"'{name}' depends on axioms: [propext, Classical.choice, Quot.sound]", [name]
        )
        self.assertEqual(set(result[name]), CHECK.ALLOWED_AXIOMS)

    def test_accepts_axiom_free_theorem(self):
        name = "ArfunFormal.Geometry.valid"
        self.assertEqual(CHECK.parse_axioms(f"'{name}' does not depend on any axioms", [name]), {name: []})

    def test_accepts_wrapped_output(self):
        name = "ArfunFormal.Geometry.valid"
        result = CHECK.parse_axioms(f"'{name}' depends on axioms:\n[propext,\n Quot.sound]", [name])
        self.assertEqual(result[name], ["Quot.sound", "propext"])

    def test_rejects_unproved_dependencies(self):
        name = "ArfunFormal.Geometry.invalid"
        with self.assertRaises(CHECK.VerificationError):
            CHECK.parse_axioms(f"'{name}' depends on axioms: [propext, sorryAx]", [name])

    def test_rejects_custom_axiom(self):
        name = "ArfunFormal.Geometry.invalid"
        with self.assertRaises(CHECK.VerificationError):
            CHECK.parse_axioms(f"'{name}' depends on axioms: [Unproved.assumption]", [name])

    def test_rejects_native_trust_extension(self):
        name = "ArfunFormal.Geometry.invalid"
        with self.assertRaises(CHECK.VerificationError):
            CHECK.parse_axioms(f"'{name}' depends on axioms: [Lean.ofReduceBool]", [name])

    def test_rejects_missing_result(self):
        with self.assertRaises(CHECK.VerificationError):
            CHECK.parse_axioms("", ["ArfunFormal.Geometry.valid"])

    def test_rejects_unexpected_result(self):
        with self.assertRaises(CHECK.VerificationError):
            CHECK.parse_axioms("'Unexpected.result' does not depend on any axioms", [])

    def test_rejects_duplicate_result(self):
        name = "ArfunFormal.Geometry.valid"
        line = f"'{name}' does not depend on any axioms\n"
        with self.assertRaises(CHECK.VerificationError):
            CHECK.parse_axioms(line + line, [name])


class RegistryTests(unittest.TestCase):
    def entry(self, name="ArfunFormal.Geometry.valid"):
        return {"name": name, "claim_es": "Una identidad.", "scope_es": "Sobre reales."}

    def test_accepts_documented_entry(self):
        entry = self.entry()
        self.assertEqual(CHECK.registry_names({"theorems": [entry]}), [entry["name"]])

    def test_rejects_empty_registry(self):
        with self.assertRaises(CHECK.VerificationError):
            CHECK.registry_names({"theorems": []})

    def test_rejects_duplicate_names(self):
        with self.assertRaises(CHECK.VerificationError):
            CHECK.registry_names({"theorems": [self.entry(), self.entry()]})

    def test_rejects_command_injection_in_name(self):
        with self.assertRaises(CHECK.VerificationError):
            CHECK.registry_names({"theorems": [self.entry("ArfunFormal.x\n#eval IO.println 1")]})

    def test_requires_scope(self):
        entry = self.entry()
        del entry["scope_es"]
        with self.assertRaises(CHECK.VerificationError):
            CHECK.registry_names({"theorems": [entry]})


class SourceCoverageTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "ArfunFormal").mkdir()
        (self.root / "ArfunFormal.lean").write_text("module\npublic import ArfunFormal.Geometry\n")
        self.source = self.root / "ArfunFormal" / "Geometry.lean"
        self.source.write_text(
            "module\npublic import Mathlib.Tactic\nnamespace ArfunFormal.Geometry\n"
            "theorem valid : True := by trivial\nend ArfunFormal.Geometry\n"
        )
        self.names = ["ArfunFormal.Geometry.valid"]

    def test_accepts_complete_source_registry(self):
        hashes, modules, imports = CHECK.inspect_sources(self.root, self.names)
        self.assertEqual(len(hashes), 2)
        self.assertEqual(modules, {"ArfunFormal", "ArfunFormal.Geometry"})
        self.assertEqual(imports, ["Mathlib.Tactic"])

    def test_rejects_unregistered_theorem(self):
        with self.assertRaises(CHECK.VerificationError):
            CHECK.inspect_sources(self.root, [])

    def test_rejects_registered_but_absent_theorem(self):
        with self.assertRaises(CHECK.VerificationError):
            CHECK.inspect_sources(self.root, self.names + ["ArfunFormal.Geometry.absent"])

    def test_rejects_module_not_imported(self):
        (self.root / "ArfunFormal.lean").write_text("module\n")
        with self.assertRaises(CHECK.VerificationError):
            CHECK.inspect_sources(self.root, self.names)

    def test_rejects_placeholder_in_source(self):
        self.source.write_text("namespace ArfunFormal.Geometry\ntheorem valid : False := by sorry\n")
        with self.assertRaises(CHECK.VerificationError):
            CHECK.inspect_sources(self.root, self.names)


class PinTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.rev = "a" * 40
        self.registry = {"lean_toolchain": "leanprover/lean4:v4.33.1", "mathlib_revision": self.rev}
        (self.root / "lean-toolchain").write_text(self.registry["lean_toolchain"] + "\n")
        (self.root / "lakefile.toml").write_text(
            'name = "test"\ndefaultTargets = ["ArfunFormal"]\n'
            '[leanOptions]\nautoImplicit = false\nwarningAsError = true\n'
            f'[[require]]\nname = "mathlib"\nrev = "{self.rev}"\n'
        )
        self.manifest = {"packages": [{"name": "mathlib", "type": "git", "rev": self.rev}]}
        self.write_manifest()

    def write_manifest(self):
        (self.root / "lake-manifest.json").write_text(json.dumps(self.manifest))

    def test_accepts_matching_pins(self):
        toolchain, deps = CHECK.validate_pins(self.root, self.registry)
        self.assertEqual(toolchain, self.registry["lean_toolchain"])
        self.assertEqual(deps, {"mathlib": self.rev})

    def test_rejects_wrong_toolchain(self):
        (self.root / "lean-toolchain").write_text("leanprover/lean4:stable\n")
        with self.assertRaises(CHECK.VerificationError):
            CHECK.validate_pins(self.root, self.registry)

    def test_rejects_mismatched_mathlib(self):
        self.manifest["packages"][0]["rev"] = "b" * 40
        self.write_manifest()
        with self.assertRaises(CHECK.VerificationError):
            CHECK.validate_pins(self.root, self.registry)

    def test_rejects_unpinned_transitive_dependency(self):
        self.manifest["packages"].append({"name": "other", "type": "git", "rev": "main"})
        self.write_manifest()
        with self.assertRaises(CHECK.VerificationError):
            CHECK.validate_pins(self.root, self.registry)


if __name__ == "__main__":
    unittest.main()
