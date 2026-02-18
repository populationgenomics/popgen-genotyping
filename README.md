# popgen-genotyping
A pipeline to process genotyping microarray data using CPG Flow and analysis-runner.

## Abstract

Lorem ipsum.

```commandline
src
├── workflow_name
│   ├── __init__.py
│   ├── config_template.toml
│   ├── jobs
│   │   └── LogicForAStage.py
│   ├── main.py
│   ├── stages.py
│   └── utils.py
```

`workflow_name` occurs in a number of places ([pyproject.toml](pyproject.toml), [src](src), and the workflow name in the
template config file). It is intended that you remove this generic placeholder name, and replace it with the name of
your workflow.

`stages.py` contains Stages in the workflow, with the actual logic imported from files in `jobs`.

`stages.py` also links to the Pipeline Naming Conventions document, containing a number of recommendations for naming
Stages and other elements of the workflow.

`config_template.toml` is a template, indicating the settings which are mandatory for the pipeline to run. In
production-pipelines, many of these settings were satisfied by the cpg-workflows or per-workflow default TOML files. If
a pipeline is being migrated from production-pipelines, the previous default config TOML would be a better substitute.
