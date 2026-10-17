"""Emit a guarded minimal patch for PANGAEA's pinned original source."""
import argparse
import difflib
from pathlib import Path
from pinned_regression_check import correct, provenance, SOURCE

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--repo",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    before=provenance(args.repo)
    after=correct(before)
    diff="".join(difflib.unified_diff(
        before.splitlines(keepends=True),after.splitlines(keepends=True),
        fromfile="a/"+str(SOURCE),tofile="b/"+str(SOURCE),n=5))
    if not diff or diff.count("@@")<2:
        raise AssertionError("Patch was unexpectedly empty or incomplete")
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(diff)
    print("PATCH_READY",args.out,"lines",len(diff.splitlines()),flush=True)

if __name__=="__main__":
    main()
