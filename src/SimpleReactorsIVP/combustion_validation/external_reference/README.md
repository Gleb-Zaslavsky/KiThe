# Independent TW-04 Reference

This Python implementation independently assembles and integrates the 14-state
TW-04 system. It does not call KiThe, RustedSciThe, or read Rust-generated
expressions. The constants below are transcribed from the reviewed TW-01/02/03
fixtures and checked against `../tw_01.rs`, `../tw_02.rs`, and `../tw_03.rs`.

From this directory on Windows:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m tw_04.calibrate
.\.venv\Scripts\python.exe -m tw_04.generate_reference
```

Calibration and fixture regeneration are explicit operations. Ordinary Rust
tests only read the checked-in frozen profile and never rewrite it. Generated
metadata records the installed Python, NumPy, and SciPy versions and solver
settings. Keep `.venv` local; do not commit it.
