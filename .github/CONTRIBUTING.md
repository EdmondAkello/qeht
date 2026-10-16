# Contributing to QEHT

Thanks for taking an interest in QEHT. Bug reports, questions, test cases and code are all welcome.

## Reporting a problem or asking a question

Open an [issue](https://github.com/EdmondAkello/qeht/issues). For a bug, include:
- the QGIS version and operating system;
- the QEHT version (Plugins > Manage and Install Plugins);
- the tool and its parameters (the settings file that the pipeline writes is ideal);
- what you expected and what happened, with the Processing log;
- if you can, a small DEM or a synthetic case that reproduces it.

Please don't attach client or project data you are not free to share. A synthetic surface that shows the same behaviour is enough.

## Changing the code

1. Open an issue first for anything beyond a small fix, so the approach can be agreed before you write it.
2. Fork the repository and work on a branch (`fix/...`, `feat/...`, `docs/...`).
3. Open a pull request against `main` that refers to the issue. CI runs the pure test suites and Bandit on every pull request, and both must pass.

### Ground rules

- **No subprocesses and no pip packages.** QEHT runs in-process on NumPy and the GDAL bindings that ship with QGIS. This is the reason the plugin exists, and it is not negotiable.
- **`core/` never imports `qgis`.** GDAL imports inside `core/` stay inside the functions that need them. Only `processing_provider/` touches QGIS.
- **The design hydrology package schema (`qeht-heas-1`) only grows.** You may add fields, layers and metadata keys, but never rename or remove one. Every new field goes into `core/interop/field_dictionary.py` with its unit, meaning, method and downstream use. Before regenerating the golden fixture, run `python3 qeht/tools/dev/golden_diff.py`: it must report `REMOVED: none` and `CHANGED: none`.
- **No invented engineering values as defaults.** A proxy is labelled PROXY and an estimate is labelled ESTIMATE, and both are recorded in the run metadata.
- **Every change ships with tests.** Use a pure suite in `tests/test_<name>.py`, in the existing `check()` style, with analytic cases where possible. Add smoke checks in `tests/qgis_smoke.py` for new Processing tools, and update the README, DOCUMENTATION.md and CHANGELOG.md.

### Running the tests

The repository folder must be named `qeht`, because the tests import `qeht.*`. You need Python with NumPy and GDAL (the Python that ships with QGIS works).

```bash
git clone https://github.com/<you>/qeht.git qeht
bash qeht/tools/dev/ci_pure.sh          # the 25 pure suites, as CI runs them
bash qeht/tools/dev/run_tests.sh        # adds the QGIS smoke and end-to-end tests (needs QGIS)
bandit -q -r qeht -x qeht/tests
```

### AI-assisted contributions

If you are contributing from outside the project and an AI tool wrote a substantial part of your change, mention it in the pull request and name the tool. You remain responsible for the change: you have read it, tested it and can explain it. This keeps the project's record of how the software was written accurate.

## Code of conduct

Everyone taking part is expected to follow the [code of conduct](CODE_OF_CONDUCT.md).

## Maintainer

QEHT is maintained by Edmond Akello. Decisions on scope and design are made in issues and pull requests, in public.
