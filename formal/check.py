import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tomllib


ALLOWED_AXIOMS = frozenset({"propext", "Classical.choice", "Quot.sound"})
NAME = re.compile(r"ArfunFormal(?:\.[A-Za-z_][A-Za-z0-9_']*)+")
NAMESPACE = re.compile(r"^namespace\s+(ArfunFormal(?:\.[A-Za-z_][A-Za-z0-9_']*)*)\s*$", re.M)
DECLARATION = re.compile(r"^(?:theorem|lemma)\s+([A-Za-z_][A-Za-z0-9_']*)\b", re.M)
IMPORT = re.compile(r"^(?:public\s+)?import\s+([A-Za-z_][A-Za-z0-9_.]*)\s*$", re.M)
FORBIDDEN = re.compile(r"\b(?:sorry|admit|axiom|unsafe|native_decide)\b|debug\.skipKernelTC")
AXIOM_OUTPUT = re.compile(
    r"'([^']+)' (?:depends on axioms:\s*\[([^\]]*)\]|does not depend on any axioms)",
    re.S,
)


class VerificationError(RuntimeError):
    pass


def registry_names(registry):
    entries = registry.get("theorems", [])
    names = [entry.get("name", "") for entry in entries]
    if not names or any(not NAME.fullmatch(name) for name in names):
        raise VerificationError("El registro debe contener nombres válidos de teoremas ArfunFormal.")
    if len(names) != len(set(names)):
        raise VerificationError("Hay teoremas duplicados en el registro.")
    for entry in entries:
        if not entry.get("claim_es") or not entry.get("scope_es"):
            raise VerificationError(f"Falta documentar el alcance de {entry['name']}.")
    return names


def parse_axioms(output, expected):
    found = {}
    for match in AXIOM_OUTPUT.finditer(output):
        name, values = match.groups()
        if name in found:
            raise VerificationError(f"Resultado duplicado de auditoría: {name}.")
        axioms = {value.strip() for value in (values or "").split(",") if value.strip()}
        extra = axioms - ALLOWED_AXIOMS
        if extra:
            raise VerificationError(f"{name} depende de axiomas no permitidos: {sorted(extra)}.")
        found[name] = sorted(axioms)
    if set(found) != set(expected):
        raise VerificationError(
            f"Auditoría incompleta: faltan {sorted(set(expected) - set(found))}; "
            f"sobran {sorted(set(found) - set(expected))}."
        )
    return found


def inspect_sources(root, expected):
    files = [root / "ArfunFormal.lean", *sorted((root / "ArfunFormal").rglob("*.lean"))]
    if len(files) < 2 or not all(path.is_file() for path in files):
        raise VerificationError("No se encontraron los módulos fuente de ArfunFormal.")
    declarations = set()
    imports = set()
    modules = set()
    hashes = {}
    for path in files:
        text = path.read_text(encoding="utf-8")
        if FORBIDDEN.search(text):
            raise VerificationError(f"Construcción no permitida en {path.relative_to(root)}.")
        hashes[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
        module = ".".join(path.relative_to(root).with_suffix("").parts)
        modules.add(module)
        imports.update(IMPORT.findall(text))
        local_names = DECLARATION.findall(text)
        if local_names:
            namespaces = NAMESPACE.findall(text)
            if namespaces != [module]:
                raise VerificationError(f"Use un único namespace {module} en {path.name}.")
            for name in local_names:
                full_name = f"{module}.{name}"
                if full_name in declarations:
                    raise VerificationError(f"Declaración duplicada: {full_name}.")
                declarations.add(full_name)
    if declarations != set(expected):
        raise VerificationError(
            f"Registro y fuentes difieren: sin registrar {sorted(declarations - set(expected))}; "
            f"sin declaración {sorted(set(expected) - declarations)}."
        )
    root_imports = set(IMPORT.findall(files[0].read_text(encoding="utf-8")))
    missing_imports = modules - {"ArfunFormal"} - root_imports
    if missing_imports:
        raise VerificationError(f"Faltan imports explícitos en ArfunFormal.lean: {sorted(missing_imports)}.")
    return hashes, modules, sorted(name for name in imports if name.startswith("Mathlib."))


def validate_pins(root, registry):
    toolchain = (root / "lean-toolchain").read_text(encoding="utf-8").strip()
    if toolchain != registry.get("lean_toolchain"):
        raise VerificationError("lean-toolchain no coincide con la versión registrada.")
    config = tomllib.loads((root / "lakefile.toml").read_text(encoding="utf-8"))
    required = [entry for entry in config.get("require", []) if entry["name"] == "mathlib"]
    if len(required) != 1 or required[0].get("rev") != registry.get("mathlib_revision"):
        raise VerificationError("La revisión de mathlib no coincide con el registro.")
    manifest = json.loads((root / "lake-manifest.json").read_text(encoding="utf-8"))
    packages = manifest.get("packages", [])
    mathlib = [package for package in packages if package["name"] == "mathlib"]
    if len(mathlib) != 1 or mathlib[0].get("rev") != registry["mathlib_revision"]:
        raise VerificationError("El lockfile no contiene la revisión fijada de mathlib.")
    if any(package.get("type") != "git" or not re.fullmatch(r"[0-9a-f]{40}", package.get("rev", ""))
           for package in packages):
        raise VerificationError("Todas las dependencias deben quedar fijadas a commits completos.")
    options = config.get("leanOptions", {})
    if options.get("warningAsError") is not True or options.get("autoImplicit") is not False:
        raise VerificationError("Se requieren warningAsError=true y autoImplicit=false.")
    if config.get("defaultTargets") != ["ArfunFormal"]:
        raise VerificationError("El objetivo predeterminado debe ser ArfunFormal.")
    return toolchain, {package["name"]: package["rev"] for package in packages}


def run(command, root, input_text=None, allow_failure=False, echo=True):
    print("+", shlex.join(command), flush=True)
    env = dict(os.environ)
    env.setdefault("LEAN_NUM_THREADS", "2")
    result = subprocess.run(
        command, cwd=root, env=env, input=input_text, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    if echo or (result.returncode and not allow_failure):
        print(result.stdout, end="", flush=True)
    if result.returncode and not allow_failure:
        raise VerificationError(f"Falló el comando anterior (salida {result.returncode}).")
    return result


def reject_controls(root):
    false_proof = "module\npublic import ArfunFormal\nexample : (0 : Nat) = 1 := by rfl\n"
    result = run(["lake", "env", "lean", "--stdin"], root, false_proof,
                 allow_failure=True, echo=False)
    if result.returncode == 0 or "error:" not in result.stdout:
        raise VerificationError("El control de prueba falsa no produjo el rechazo esperado.")
    controls = {"false_proof": "rejected"}
    probes = {
        "incomplete_proof": (
            "theorem pending : False := by sorry\n",
            "ArfunAuditProbe.pending",
            "sorryAx",
        ),
        "custom_axiom_dependency": (
            "axiom invented : False\ntheorem inherited : False := invented\n",
            "ArfunAuditProbe.inherited",
            "ArfunAuditProbe.invented",
        ),
    }
    for label, (body, name, expected_axiom) in probes.items():
        text = ("module\npublic import ArfunFormal\npublic section\nnamespace ArfunAuditProbe\n"
                + body + f"end ArfunAuditProbe\n#print axioms {name}\n")
        result = run(["lake", "env", "lean", "--stdin"], root, text, echo=False)
        if expected_axiom not in result.stdout:
            raise VerificationError(f"El control {label} no expuso el axioma esperado.")
        try:
            parse_axioms(result.stdout, [name])
        except VerificationError:
            controls[label] = "rejected"
        else:
            raise VerificationError(f"La auditoría aceptó el control no válido {label}.")
    return controls


def main(argv=None):
    parser = argparse.ArgumentParser(description="Verificar las pruebas formales de arfun.")
    parser.add_argument("--bootstrap", action="store_true",
                        help="Descargar la caché oficial de los imports de mathlib antes de compilar.")
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parent
    report_path = root / ".lake" / "verification-report.json"
    report = {"status": "failed", "scope": "Mathematical specification over real numbers; not Stan execution."}
    try:
        if shutil.which("lake") is None:
            raise VerificationError("No se encuentra Lake. Consulte MANUAL.es.md para instalar Elan.")
        registry = json.loads((root / "theorems.json").read_text(encoding="utf-8"))
        names = registry_names(registry)
        toolchain, dependencies = validate_pins(root, registry)
        hashes, modules, imports = inspect_sources(root, names)
        run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"], root)
        if args.bootstrap:
            run(["lake", "exe", "cache", "get", *imports], root)
        run(["lake", "--wfail", "build"], root)
        version = run(["lake", "env", "lean", "--short-version"], root).stdout.strip()
        if version != toolchain.rsplit(":", 1)[1].removeprefix("v"):
            raise VerificationError("La versión ejecutada de Lean no coincide con lean-toolchain.")
        commands = "module\npublic import ArfunFormal\n" + "\n".join(
            f"#print axioms {name}" for name in names
        ) + "\n"
        audit = run(["lake", "env", "lean", "-DwarningAsError=true", "--stdin"], root, commands)
        axioms = parse_axioms(audit.stdout, names)
        replay = run(["lake", "env", "leanchecker", "-v", "ArfunFormal"], root)
        checked_modules = set(re.findall(r"^replaying (\S+)\s*$", replay.stdout, re.M))
        if not modules.issubset(checked_modules):
            raise VerificationError(f"leanchecker no revisó todos los módulos: {sorted(modules - checked_modules)}.")
        controls = reject_controls(root)
        report.update(status="passed", lean_version=version, mathlib_revision=registry["mathlib_revision"],
                      dependencies=dependencies, theorem_count=len(names), axioms=axioms,
                      source_sha256=hashes, replayed_modules=sorted(checked_modules), negative_controls=controls)
        print(f"\nVERIFICADO: {len(names)} teoremas; {len(modules)} módulos; "
              f"{len(controls)} controles negativos rechazados.")
        print("Alcance: especificación matemática; no valida datos ecológicos ni todo el ejecutable Stan.")
    except Exception as error:
        report["error"] = str(error)
        raise
    finally:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print("Informe:", report_path)


if __name__ == "__main__":
    try:
        main()
    except (VerificationError, OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
