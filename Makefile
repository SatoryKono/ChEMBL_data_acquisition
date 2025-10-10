SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c
PYTHON ?= python3
VENV := .venv
PYTHON_BIN := $(VENV)/bin/python
PIP := $(VENV)/bin/pip
ACTIVITY_LIMIT ?= 25

ENV_FILES := $(wildcard .env .env.local)
ifneq ($(ENV_FILES),)
  include $(ENV_FILES)
  export $(shell sed -n 's/^\([A-Za-z_][A-Za-z0-9_]*\)=.*/\1/p' $(ENV_FILES))
endif

.PHONY: init lint test smoke test-report get-activities build release clean

init: $(PYTHON_BIN)

$(PYTHON_BIN): requirements.txt pyproject.toml $(wildcard .python-version)
	@if [ -f .python-version ]; then \
		$(PYTHON) -c "import pathlib, sys; expected = pathlib.Path('.python-version').read_text().strip(); actual = '.'.join(map(str, sys.version_info[:3])); assert actual == expected, f'Python {expected} required, but {actual} found'"; \
	fi
	$(PYTHON) -m venv $(VENV)
	$(PIP) install --upgrade pip
	$(PIP) install '.[dev]'

lint: $(PYTHON_BIN)
	$(VENV)/bin/ruff check .
	$(VENV)/bin/black --check .
	$(VENV)/bin/mypy

test: $(PYTHON_BIN)
	$(PYTHON_BIN) scripts/run_tests.py

smoke: $(PYTHON_BIN)
	CHEMBL_DA_BASE_PATH=$(PWD)/tests/resources/pipeline_inputs $(VENV)/bin/pytest tests/smoke -k "not testitem"

test-report: $(PYTHON_BIN)
	PYTHONHASHSEED=$${PYTHONHASHSEED:-0} \
	CHEMBL_DA_BASE_PATH=$(PWD)/tests/resources/pipeline_inputs \
	$(PYTHON_BIN) scripts/run_tests.py

get-activities: $(PYTHON_BIN)
	$(PYTHON_BIN) scripts/get_activities.py --limit $(ACTIVITY_LIMIT) --dry-run

build: $(PYTHON_BIN)
	rm -rf dist
	$(PYTHON_BIN) -m build

release: build
	$(PYTHON_BIN) -m twine check dist/*

clean:
	rm -rf $(VENV) build dist .mypy_cache .pytest_cache
