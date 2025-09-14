# Configuration Reference

Все скрипты читают настройки из `config.yaml` в корне проекта.
Структура файла валидируется по `config.schema.json`.  
Параметры могут быть переопределены переменными окружения
и флагами командной строки.

## Основные разделы

### `api`

| Ключ | Значение по умолчанию | Описание |
|------|-----------------------|----------|
| `chembl_base_url` | `https://www.ebi.ac.uk/chembl/api/data` | Базовый URL ChEMBL API |
| `pubchem_base_url` | `https://pubchem.ncbi.nlm.nih.gov/rest/pug` | REST‑сервисы PubChem |
| `uniprot_base_url` | `https://rest.uniprot.org` | API UniProt |
| `eutils_base_url` | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils` | NCBI E‑utilities |
| `semanticscholar_base_url` | `https://api.semanticscholar.org/graph/v1` | Semantic Scholar |
| `openalex_base_url` | `https://api.openalex.org` | OpenAlex |
| `crossref_base_url` | `https://api.crossref.org` | CrossRef |
| `gtp_base_url` | `https://www.guidetopharmacology.org` | Guide to Pharmacology |

### `timeouts`

| Ключ | Значение по умолчанию | Описание |
|------|-----------------------|----------|
| `connect` | `10` | таймаут установления соединения (сек) |
| `read`    | `30` | таймаут ожидания ответа (сек) |

### `rate_limits`

| Ключ | Значение по умолчанию | Описание |
|------|-----------------------|----------|
| `max_requests_per_second` | `5`   | лимит запросов в секунду |
| `max_retries`             | `3`   | число повторов при временных ошибках |
| `backoff_factor`          | `0.3` | коэффициент экспоненциального ожидания |

### `chembl`

| Ключ | Значение по умолчанию | Описание |
|------|-----------------------|----------|
| `cache_ttl`     | `3600` | время жизни кэша ответов API (сек) |
| `cache_maxsize` | `1024` | максимальное число записей в кэше |

### `output`

| Ключ | Значение по умолчанию | Описание |
|------|-----------------------|----------|
| `data_dir` | `data` | каталог выходных данных |
| `logs_dir` | `logs` | каталог логов |
| `tmp_dir`  | `tmp`  | временные файлы |

## Переопределения через окружение

Шаблон переменной:

```
CHEMBL_DA__<SECTION>__<KEY>=value
```

Примеры:

```bash
export CHEMBL_DA__API__CHEMBL_BASE_URL=https://mirror/api
export CHEMBL_DA__OUTPUT__DATA_DIR=/mnt/datasets
```

Короткие алиасы:

| Переменная            | Соответствие              |
|-----------------------|---------------------------|
| `CHEMBL_DA_BASE`      | `CHEMBL_DA__API__CHEMBL_BASE_URL` |
| `CHEMBL_DA_TIMEOUT_CONNECT` | `CHEMBL_DA__API__CONNECT` |
| `CHEMBL_DA_TIMEOUT_READ`    | `CHEMBL_DA__API__READ`    |
| `CHEMBL_DA_OUTDIR`    | `CHEMBL_DA__OUTPUT__DATA_DIR` |
| `CHEMBL_DA_LOG_LEVEL` | `CHEMBL_DA__LOG__LEVEL` |

Неизвестные переменные игнорируются с предупреждением.

## Переопределения через CLI

Каждый ключ может быть указан в виде `--section.key=value`, например:

```bash
python scripts/get_activity_data.py \
    --api.chembl_base_url https://mirror/api \
    --output.data_dir results
```

CLI‑параметры имеют самый высокий приоритет, затем переменные окружения,
последним — значения из `config.yaml`.

## Советы по использованию

* При работе на многопользовательских системах переназначайте `output.data_dir`.
* Для краткосрочных экспериментов можно опустить `config.yaml` —
  скрипты используют значения по умолчанию.
* Все пути в конфигурации интерпретируются относительно корня проекта,
  если не указано иное.

