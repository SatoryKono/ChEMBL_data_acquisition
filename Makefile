SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c
PYTHON ?= python3
VENV := .venv
PYTHON_BIN := $(VENV)/bin/python
PIP := $(VENV)/bin/pip
DEV_SENTINEL := $(VENV)/.dev-deps-installed
ACTIVITY_LIMIT ?= 25

DOC_PROTOCOL_SOURCES := docs/en/PROTOCOL_EN.md docs/ru/PROTOCOL_RU.md
DOC_PROTOCOL_OUTPUT := docs/ChEMBL_Data_Acquisition_Protocol_v2.1.docx

ENV_FILES := $(wildcard .env .env.local)
ifneq ($(ENV_FILES),)
  include $(ENV_FILES)
  export $(shell sed -n 's/^\([A-Za-z_][A-Za-z0-9_]*\)=.*/\1/p' $(ENV_FILES))
endif

.PHONY: init lint test smoke test-report get-activities build release clean

init: $(PYTHON_BIN)

$(PYTHON_BIN): requirements.txt pyproject.toml
	@if [ -f .python-version ]; then \
		$(PYTHON) -c "import pathlib, sys; expected = pathlib.Path('.python-version').read_text().strip(); actual = '.'.join(map(str, sys.version_info[:3])); assert actual == expected, f'Python {expected} required, but {actual} found'"; \
	fi
	$(PYTHON) -m venv $(VENV)
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements-lock.txt
	$(PIP) install --no-deps -e .

$(DEV_SENTINEL): $(PYTHON_BIN) requirements-dev.txt
	$(PIP) install -r requirements-dev.txt
	touch $(DEV_SENTINEL)

lint: $(PYTHON_BIN)
	$(VENV)/bin/ruff check .
	$(VENV)/bin/black --check .
	$(VENV)/bin/mypy

test: $(DEV_SENTINEL)
	$(PYTHON_BIN) scripts/run_tests.py

smoke: $(PYTHON_BIN)
        PYTHONHASHSEED=$${PYTHONHASHSEED:-0} \
        CHEMBL_DA_BASE_PATH=$(PWD)/tests/resources/pipeline_inputs \
        $(VENV)/bin/pytest -m "smoke"

test-report: $(DEV_SENTINEL)
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
	rm -rf $(VENV) build dist .pytest_cache
	find . -name '.mypy_cache' -prune -exec rm -rf {} +
	find . -type d -name '__pycache__' -prune -exec rm -rf {} +

.PHONY: protocol-docx
protocol-docx: $(PYTHON_BIN)
	$(PYTHON_BIN) scripts/convert_md_to_docx.py $(DOC_PROTOCOL_SOURCES) \
	  -o $(DOC_PROTOCOL_OUTPUT) --number-sections --toc
	@echo "DOCX artefact ready at $(DOC_PROTOCOL_OUTPUT) (git-ignored deliverable)."
