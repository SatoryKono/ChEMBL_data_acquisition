# Configuration Reference

Все скрипты читают настройки из `config.yaml` в корне проекта.
Файл валидируется по `config.schema.json`, поэтому любые неизвестные ключи
считаются ошибкой. Параметры можно переопределять через переменные окружения и
флаги командной строки — приоритет следующий: `config.yaml` < окружение < CLI.

## Структура `config.yaml`

| Раздел        | Назначение                                                                      |
|---------------|---------------------------------------------------------------------------------| 
| `sources`     | Подключение к удалённым источникам данных (ChEMBL, UniProt, CrossRef и др.).     |
| `local`       | Пути к локальным ресурсам и каталогам, используемым при обработке.              |
| `system`      | Общие настройки журнала, повторов и агрегации документации.                     |

Ниже перечислены ключевые параметры. Если значение чувствительно (ключ, токен,
контактный e-mail), задавайте его через переменную окружения.

### `sources.chembl`

| Ключ                               | Значение по умолчанию | Описание |
|------------------------------------|-----------------------|----------|
| `api.chembl_base`                  | `https://www.ebi.ac.uk/chembl/api/data` | Базовый URL API. |
| `api.user_agent`                   | `chembl-da/0.1 (mailto:contact@example.org)` | User-Agent с контактами. Настройте через `CHEMBL_DA_BASE`/`CHEMBL_DA_USER_AGENT`. |
| `api.rps` / `api.burst`            | `20`                  | Ограничения скорости для лимитера. |
| `cache.cache_ttl` / `cache_maxsize`| `3600` / `1024`       | TTL и размер кэша ответов API. |
| `pipelines.activity.*`             | см. значения          | Параметры загрузки активностей. |
| `pipelines.assay.*`                | см. значения          | Параметры загрузки ассайев. |
| `pipelines.testitem.*`             | см. значения          | Параметры загрузки тест-объектов. |
| `pipelines.document.*`             | см. значения          | Настройки выгрузки публикаций ChEMBL/PubMed. |
| `pipelines.target.*`               | см. значения          | Пути и ограничения для модулей целей. |

### `sources.openalex`, `sources.crossref`

Таблицы имеют одинаковую структуру: `base`, `timeout_connect`, `timeout_read`,
`retries`, `rps`, `burst`, `mailto`. Персональный e-mail обязателен.

### `sources.uniprot`

| Ключ                     | Значение по умолчанию | Описание |
|--------------------------|-----------------------|----------|
| `api.base`               | `https://rest.uniprot.org` | REST API UniProt. |
| `api.delay`              | `0.25`                | Пауза между запросами. |
| `mapping.base`           | `https://rest.uniprot.org/idmapping` | Сервис ID‑маппинга. |
| `mapping.poll_interval`  | `0.5`                 | Интервал повторного опроса задачи. |
| `mapping.timeout`        | `300`                 | Максимальное время ожидания результата. |
| `mapping.cache_maxsize`  | `128`                 | Размер кэша ответов. |

### `sources.iuphar`, `sources.pubchem`, `sources.pubmed`, `sources.semantic_scholar`

Каждый блок содержит базовый URL и ограничения по времени/скорости. Для
IUPHAR используются глобальные правила ретраев (`system.retry`), поэтому в
секцию источника вынесены только лимиты RPS/`burst`. Для PubChem
дополнительно задаётся `delay` и `cache_ttl` (секунды). Для PubMed fallback
кодировки CSV (`utf-8-sig`, `cp1251`, `latin1`) зашиты в библиотеку и не
редактируются через конфигурацию.

### `local.resources`

| Ключ              | Значение по умолчанию         | Описание |
|-------------------|-------------------------------|----------|
| `dictionary_dir`  | `dictionary`                  | Корень словарей. |
| `iuphar_target_csv` / `iuphar_family_csv` | см. значения | CSV c сопоставлениями IUPHAR. |
| `uniprot_data_dir`| `dictionary/uniprot`          | Кэш UniProt JSON. |
| `organism_csv`    | `dictionary/_Target/targets_type.csv` | Справочник организмов. |

### `local.io`

| Ключ          | Значение по умолчанию | Описание |
|---------------|-----------------------|----------|
| `output_dir`  | `data/output`         | Каталог для итоговых файлов. |
| `cache_dir`   | `.cache`              | Каталог HTTP-кэша. |
| `csv_sep`     | `,`                   | Разделитель в CSV. |
| `csv_encoding`| `utf-8-sig`           | Кодировка CSV. |
| `exist_ok`    | `true`                | Создавать каталоги при отсутствии. |

### `local.init`

Содержит пути к исходным Excel-файлам и директорию вывода подготовленных данных.

### `system`

| Ключ                     | Значение по умолчанию | Описание |
|--------------------------|-----------------------|----------|
| `log.level`              | `INFO`                | Уровень логирования. |
| `log.format` / `datefmt` | см. значения          | Формат сообщений. |
| `rate.global_rps/burst`  | `100`                 | Глобальные лимиты для токен-бакета. |
| `rate.limiter_cache_*`   | `128` / `600`         | Настройки кэша лимитеров. |
| `retry.max_attempts`     | `3`                   | Число повторов при ошибках. |
| `retry.backoff_factor`   | `0.5`                 | Экспоненциальная задержка. |
| `retry.status_forcelist` | `[429,500,502,503,504]` | Коды для повторов. |
| `doc_type.*`             | см. значения          | Весовые коэффициенты классификации документов. |

## Переменные окружения

Переменные формируются по шаблону `CHEMBL_DA__<SECTION>__...__<KEY>`.
Например:

```bash
export CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT="project/1.0 (mailto:me@example.com)"
export CHEMBL_DA__LOCAL__IO__OUTPUT_DIR=/mnt/datasets
```

Для часто используемых настроек доступны короткие алиасы:

| Переменная            | Соответствие                                      |
|-----------------------|---------------------------------------------------|
| `CHEMBL_DA_BASE`      | `sources.chembl.api.chembl_base`                  |
| `CHEMBL_DA_RPS`       | `sources.chembl.api.rps`                          |
| `CHEMBL_DA_BURST`     | `sources.chembl.api.burst`                        |
| `CHEMBL_DA_TIMEOUT_CONNECT` | `sources.chembl.api.timeout_connect`       |
| `CHEMBL_DA_TIMEOUT_READ`    | `sources.chembl.api.timeout_read`          |
| `CHEMBL_DA_CACHE_DIR` | `local.io.cache_dir`                              |
| `CHEMBL_DA_OUTDIR`    | `local.io.output_dir`                             |
| `CHEMBL_DA_CACHE_TTL` | `sources.chembl.cache.cache_ttl`                  |
| `CHEMBL_DA_CACHE_MAXSIZE` | `sources.chembl.cache.cache_maxsize`         |
| `CHEMBL_DA_GLOBAL_RPS`    | `system.rate.global_rps`                     |
| `CHEMBL_DA_GLOBAL_BURST`  | `system.rate.global_burst`                   |

## CLI-переопределения

Любой параметр можно передать в виде `--path.to.value`. Примеры:

```bash
python scripts/get_activity_data.py \
  --sources.chembl.api.user_agent "project/1.0 (mailto:me@example.com)" \
  --system.retry.max_attempts 5
```

## Рекомендации по безопасности

* Не коммитьте персональные e-mail или ключи в `config.yaml` — используйте `.env`
  или переменные окружения.
* Каталоги из `local.*` лучше переназначать на изолированное хранилище при
  работе в многопользовательской среде.
* Значения `rps`/`burst` соблюдайте в соответствии с политиками провайдеров.
* В репозитории есть `.env.example` с базовым набором переменных для контактов
  API. Скопируйте его в `.env` и укажите свои значения.
