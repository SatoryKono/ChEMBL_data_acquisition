SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c
PYTHON_VERSION := $(shell cat .python-version)
PYTHON ?= python3
VENV := .venv
PYTHON_BIN := $(VENV)/bin/python
PIP := $(VENV)/bin/pip

ENV_FILES := $(wildcard .env .env.local)
ifneq ($(ENV_FILES),)
  include $(ENV_FILES)
  export $(shell sed -n 's/^\([A-Za-z_][A-Za-z0-9_]*\)=.*/\1/p' $(ENV_FILES))
endif

.PHONY: init lint test smoke build release clean

init: $(PYTHON_BIN)

$(PYTHON_BIN): .python-version requirements.txt pyproject.toml
	@test -f .python-version
	@$(PYTHON) -c 'import pathlib, sys; expected = pathlib.Path(".python-version").read_text().strip(); actual = ".".join(map(str, sys.version_info[:3])); assert actual == expected, f"Python {expected} required, but {actual} found"'
	$(PYTHON) -m venv $(VENV)
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt
	$(PIP) install -e .[dev]

lint: $(PYTHON_BIN)
	$(VENV)/bin/ruff check .
	$(VENV)/bin/black --check .
	$(VENV)/bin/mypy

test: $(PYTHON_BIN)
	$(VENV)/bin/pytest

smoke: $(PYTHON_BIN)
	CHEMBL_DA_BASE_PATH=$(PWD)/tests/data $(VENV)/bin/pytest tests/smoke -k "not testitem"

build: $(PYTHON_BIN)
	rm -rf dist
	$(PYTHON_BIN) -m build

release: build
	$(PYTHON_BIN) -m twine check dist/*

clean:
	rm -rf $(VENV) build dist .mypy_cache .pytest_cache
