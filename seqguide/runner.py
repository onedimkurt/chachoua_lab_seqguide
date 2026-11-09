"""Runs workflow steps, logs everything and records provenance (tool versions, exact commands)."""

from __future__ import annotations

import json
import shlex
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

VERSION_CMDS = {
    "fastp": ["fastp", "--version"],
    "seqkit": ["seqkit", "version"],
    "unicycler": ["unicycler", "--version"],
    "spades.py": ["spades.py", "--version"],
    "flye": ["flye", "--version"],
    "medaka_consensus": ["medaka", "--version"],
    "minimap2": ["minimap2", "--version"],
    "samtools": ["samtools", "--version"],
    "bcftools": ["bcftools", "--version"],
}


@dataclass
class Step:
    name: str
    explanation: str                       # plain-language "what & why", shown in GUI and report
    cmds: list = field(default_factory=list)  # list of argv lists; >1 means a pipe cmd1 | cmd2 | ...
    stdout_to: Optional[str] = None        # write final stdout to this file
    func: Optional[Callable] = None        # Python step (called with no args)
    tools: list = field(default_factory=list)
    optional: bool = False                 # failure produces a warning instead of stopping


class StepFailed(RuntimeError):
    def __init__(self, step: Step, detail: str):
        super().__init__(f"Step '{step.name}' failed. {detail}")
        self.step = step
        self.detail = detail


class Runner:
    def __init__(self, workdir: Path, log: Optional[Callable[[str], None]] = None):
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.logfile = self.workdir / "run.log"
        self._log_cb = log
        self.records: list = []
        self.warnings: list = []

    def log(self, msg: str) -> None:
        stamp = time.strftime("%H:%M:%S")
        with open(self.logfile, "a") as fh:
            fh.write(f"[{stamp}] {msg}\n")
        if self._log_cb:
            self._log_cb(msg)

    @staticmethod
    def missing_tools(steps: list) -> list:
        need = sorted({t for s in steps for t in s.tools})
        return [t for t in need if shutil.which(t) is None]

    @staticmethod
    def tool_versions(tools: list) -> dict:
        out = {}
        for t in sorted(set(tools)):
            cmd = VERSION_CMDS.get(t, [t, "--version"])
            try:
                p = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                txt = (p.stdout or p.stderr).strip().splitlines()
                out[t] = txt[0] if txt else "unknown"
            except (OSError, subprocess.TimeoutExpired):
                out[t] = "not available"
        return out

    def run_one(self, step: Step, all_steps: list) -> None:
        """Run a single step (used by the GUI to show progress between steps)."""
        self.log(f"▶ {step.name}")
        t0 = time.time()
        rec = {"step": step.name, "explanation": step.explanation,
               "commands": [" | ".join(shlex.join(c) for c in step.cmds)] if step.cmds else ["(internal Python step)"]}
        try:
            if step.func:
                step.func()
            else:
                self._run_pipe(step)
            rec["status"] = "ok"
        except Exception as exc:  # noqa: BLE001
            rec["status"] = "failed"
            rec["error"] = str(exc)
            self.records.append(rec)
            if step.optional:
                self.warnings.append(f"{step.name}: {exc}")
                self.log(f"  ⚠ optional step failed: {exc}")
                return
            self.save_provenance(all_steps)
            if isinstance(exc, StepFailed):
                raise
            raise StepFailed(step, str(exc)) from exc
        rec["seconds"] = round(time.time() - t0, 1)
        self.records.append(rec)
        self.log(f"  ✓ done in {rec['seconds']} s")

    def run(self, steps: list, on_step: Optional[Callable[[int, int, Step], None]] = None) -> None:
        for i, step in enumerate(steps):
            if on_step:
                on_step(i, len(steps), step)
            self.run_one(step, steps)
        self.save_provenance(steps)

    def _run_pipe(self, step: Step) -> None:
        errfile = self.workdir / "logs" / (step.name.lower().replace(" ", "_") + ".stderr.txt")
        errfile.parent.mkdir(exist_ok=True)
        procs = []
        with open(errfile, "w") as err:
            prev = None
            for j, cmd in enumerate(step.cmds):
                last = j == len(step.cmds) - 1
                out = open(self.workdir / step.stdout_to, "wb") if (last and step.stdout_to) else (
                    subprocess.DEVNULL if last else subprocess.PIPE)
                p = subprocess.Popen(cmd, stdin=prev, stdout=out, stderr=err, cwd=self.workdir)
                if prev is not None:
                    prev.close()
                prev = p.stdout
                procs.append((p, out))
            codes = []
            for p, out in procs:
                codes.append(p.wait())
                if hasattr(out, "close") and out not in (subprocess.PIPE, subprocess.DEVNULL):
                    out.close()
        if any(codes):
            tail = errfile.read_text(errors="replace").strip().splitlines()[-15:]
            raise StepFailed(step, "Tool message:\n" + "\n".join(tail))

    def save_provenance(self, steps: list) -> None:
        tools = [t for s in steps for t in s.tools]
        versions = self.tool_versions(tools) if tools else {}
        try:
            from . import __version__
            versions["SeqGuide (built-in engine)" if not tools else "SeqGuide"] = __version__
        except ImportError:  # pragma: no cover
            pass
        data = {"steps": self.records, "tool_versions": versions, "warnings": self.warnings}
        (self.workdir / "provenance.json").write_text(json.dumps(data, indent=2))
