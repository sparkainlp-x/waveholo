# Reproducibility

## 5-minute reproduction
1. Clone: git clone https://github.com/sparkainlp-x/waveholo
2. Install: make install
3. Test: make test
4. Docs: make docs
5. Smoke: make smoke (expected exit code 2 with reduced seeds; checksums verify)

## Evidence status
- Historical slope check: REPORTED (6 stage positions)
- Full factorial benchmark: UNRUN
- Independent real-data rerun: UNRUN

## Claim gates
Smoke runner intentionally returns 2 when seed count does not meet the claim gate. This is fail-closed behavior, not a bug.
