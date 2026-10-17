"""
MODULE: scripts/modal_cost.py
Lightweight cost tracking for Modal GPU usage across this repo's modal_*.py
scripts. Not billing-API-accurate (Modal's CLI doesn't expose a billing/usage
endpoint as of this writing) -- estimates cost as GPU-second rate x observed
container-seconds, which callers report themselves (see record_usage).

Rates below are Modal's published per-second GPU pricing (modal.com/pricing,
captured 2026-07-23) converted to $/hour. Re-check the pricing page if this
file gets old -- rates are not queried live.

INDEX:
  GPU_HOURLY_USD   -- $/GPU-hour by tier
  estimate_cost    -- (gpu, container_seconds, n_containers) -> USD
  record_usage     -- append one usage record to the JSONL log, return its cost
  print_summary    -- read the log, print per-script/per-gpu totals + grand total

CLI:
    python3 scripts/modal_cost.py                  # print summary of the default log
    python3 scripts/modal_cost.py --log <path>      # summarize a specific log
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_LOG = REPO_ROOT / "runs" / "modal_cost_log.jsonl"

# $/GPU-hour, derived from Modal's published $/sec rates x 3600.
GPU_HOURLY_USD = {
    "T4": 0.59,
    "L4": 0.80,
    "A10": 1.10,
    "L40S": 1.95,
    "A100-40GB": 2.10,
    "A100-80GB": 2.50,
    "RTX PRO 6000": 3.03,
    "H100": 3.95,
    "H200": 4.54,
    "B200": 6.25,
    "B300": 7.10,
}


def estimate_cost(gpu: str, container_seconds: float, n_containers: int = 1) -> float:
    """Cost = rate($/hr) / 3600 * container_seconds * n_containers.

    container_seconds is PER-CONTAINER wall time (e.g. one refold_score call's
    duration), not the local entrypoint's total wall-clock -- a fan-out of 10
    concurrent containers for 60s bills as 600 container-seconds, not 60."""
    rate = GPU_HOURLY_USD.get(gpu)
    if rate is None:
        raise ValueError(f"unknown GPU tier {gpu!r}; known: {sorted(GPU_HOURLY_USD)}")
    return round(rate / 3600.0 * container_seconds * n_containers, 4)


def record_usage(script: str, gpu: str, container_seconds: float, n_containers: int = 1,
                 entrypoint: str = "", app_id: str = "", note: str = "",
                 estimated: bool = True, log_path: Path = DEFAULT_LOG) -> dict:
    """Append one usage record to the JSONL cost log. Returns the record
    (including its computed cost) for the caller to print/report."""
    cost = estimate_cost(gpu, container_seconds, n_containers)
    record = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "script": script,
        "entrypoint": entrypoint,
        "gpu": gpu,
        "container_seconds": round(container_seconds, 1),
        "n_containers": n_containers,
        "gpu_hours": round(container_seconds * n_containers / 3600.0, 4),
        "rate_usd_per_hour": GPU_HOURLY_USD[gpu],
        "cost_usd": cost,
        "estimated": estimated,  # False only if sourced from Modal's own billed figure
        "app_id": app_id,
        "note": note,
    }
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as f:
        f.write(json.dumps(record) + "\n")
    return record


def _load(log_path: Path) -> list:
    if not log_path.exists():
        return []
    return [json.loads(line) for line in log_path.read_text().splitlines() if line.strip()]


def print_summary(log_path: Path = DEFAULT_LOG) -> None:
    records = _load(log_path)
    if not records:
        print(f"no usage records at {log_path}")
        return

    print(f"=== Modal GPU cost log: {log_path} ({len(records)} record(s)) ===\n")
    print(f"{'timestamp':<20} {'script':<28} {'gpu':<14} {'gpu_hrs':>8} {'cost_usd':>9}  note")
    for r in records:
        flag = "" if not r.get("estimated", True) else " (est.)"
        print(f"{r['timestamp']:<20} {r['script']:<28} {r['gpu']:<14} "
              f"{r['gpu_hours']:>8.4f} {r['cost_usd']:>9.4f}{flag}  {r.get('note', '')}")

    by_gpu: dict = {}
    by_script: dict = {}
    grand_total = 0.0
    for r in records:
        by_gpu.setdefault(r["gpu"], 0.0)
        by_gpu[r["gpu"]] += r["cost_usd"]
        by_script.setdefault(r["script"], 0.0)
        by_script[r["script"]] += r["cost_usd"]
        grand_total += r["cost_usd"]

    print("\n--- by GPU tier ---")
    for gpu, cost in sorted(by_gpu.items(), key=lambda x: -x[1]):
        rate = GPU_HOURLY_USD.get(gpu)
        print(f"  {gpu:<14} ${cost:.4f}  (rate: ${rate:.2f}/hr)" if rate else f"  {gpu:<14} ${cost:.4f}")

    print("\n--- by script ---")
    for script, cost in sorted(by_script.items(), key=lambda x: -x[1]):
        print(f"  {script:<28} ${cost:.4f}")

    print(f"\nGRAND TOTAL: ${grand_total:.4f}")
    print("\n--- $/GPU-hour reference (Modal, captured 2026-07-23) ---")
    for gpu, rate in sorted(GPU_HOURLY_USD.items(), key=lambda x: x[1]):
        print(f"  {gpu:<14} ${rate:.2f}/hr")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Summarize Modal GPU usage cost log.")
    p.add_argument("--log", default=str(DEFAULT_LOG), help="path to the JSONL cost log")
    args = p.parse_args(argv)
    print_summary(Path(args.log))
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
