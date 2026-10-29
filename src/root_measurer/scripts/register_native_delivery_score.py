"""Attach an immutable independent score to its exact local delivery."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from root_model_delivery_validation import register_score

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', required=True, type=Path)
    parser.add_argument('--score', required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(register_score(args.prefix, args.score), ensure_ascii=False))
