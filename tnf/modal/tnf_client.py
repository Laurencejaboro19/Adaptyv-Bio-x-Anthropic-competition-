"""PDFA TNF launcher. Talks to the DEPLOYED Modal app 'pdfa-tnf-dual'. Run from the tnf_modal folder.
  python tnf_client.py start  --tag tnf --max-shards 1 --per-round 1 --goal 1 --budget 12
  python tnf_client.py start  --tag tnf --max-shards 14 --per-round 3 --goal 24 --budget 85
  python tnf_client.py status --tag tnf      |  logs --tag tnf  |  fetch --tag tnf"""
import argparse, json, pathlib, sys
import modal

APP = "pdfa-tnf-dual"
def F(name):
    return modal.Function.from_name(APP, name)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["check", "start", "status", "logs", "fetch"])
    ap.add_argument("--tag", default="tnf")
    ap.add_argument("--goal", type=int, default=24)
    ap.add_argument("--max-shards", type=int, default=14)
    ap.add_argument("--per-round", type=int, default=3)
    ap.add_argument("--budget", type=float, default=85.0)
    ap.add_argument("--deadline-h", type=float, default=14.0)
    ap.add_argument("--gate-json", default="")
    a = ap.parse_args()

    if a.cmd == "check":
        ok = True
        for f in ("tnf_human_AB.pdb", "tnf_mouse_AB.pdb", "tnf_human_trimer_12on.pdb",
                  "tnf_mouse_trimer_12on.pdb", "tnf_app.py", "tnf_client.py"):
            p = pathlib.Path(f)
            if not p.exists(): print("MISSING:", f); ok = False; continue
            if f.endswith(".pdb"):
                ch = {}
                for l in p.read_text().splitlines():
                    if l.startswith("ATOM") and l[12:16].strip() == "CA":
                        ch[l[21]] = ch.get(l[21], 0) + 1
                want = {"A": 146, "B": 146, "C": 146} if "human" in f else {"A": 145, "B": 145, "C": 145}
                if "_AB" in f: want.pop("C")
                good = ch == want
                ok &= good
                print(("OK   " if good else "BAD  ") + f, ch, "expected", want)
            else:
                print("OK   " + f, p.stat().st_size, "bytes")
        print("ALL FILES GOOD" if ok else "PROBLEM - fix before deploying")
    elif a.cmd == "start":
        call = F("controller").spawn(a.tag, a.goal, a.max_shards, a.per_round, a.budget,
                                     a.deadline_h, a.gate_json)
        print("CONTROLLER STARTED IN THE CLOUD  call id:", call.object_id, " tag:", a.tag)
        print("You can close this window and switch the PC off now.")
    elif a.cmd == "status":
        print(F("read_progress").remote(a.tag))
    elif a.cmd == "logs":
        print(F("tail_logs").remote(a.tag))
    elif a.cmd == "fetch":
        for n in ("designs.csv", "all_scored.csv", "top20.fasta", "progress.json"):
            b = F("read_file").remote(a.tag, n)
            if b:
                pathlib.Path(a.tag + "_" + n).write_bytes(b); print("wrote", a.tag + "_" + n, len(b), "bytes")
            else:
                print("not yet:", n)

if __name__ == "__main__":
    main()
