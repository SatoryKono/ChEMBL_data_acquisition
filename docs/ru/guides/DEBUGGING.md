# Памятка по отладке

Частые проблемы и способы их устранения при работе с конвейерами.

## 1. Ошибки CLI

- `get_document_data` требует явного режима — передавайте `--mode <chembl|pubmed|all>`
  или позиционную команду совместимости. `--mode all` (или позиционная команда
  `all`) запускает объединённый прогон. `get_target_data` по-прежнему требует
  явную подкоманду (`chembl`, `uniprot`, `iuphar`, `all`).
- Флаг `--limit 0` полностью отключает конвейер — удобен для пропуска этапа, но
  может стать неожиданностью.
- Используйте абсолютные пути или комбинацию `--base-path` + `--input-dir` / `--output-dir`,
  чтобы не зависеть от текущего каталога.

## 2. Проблемы с кодировкой и разделителями

Симптомы: Pandera сообщает об отсутствующих колонках, `csv.Error: iterator should
return strings, not bytes`.

- Проверьте кодировку: `file -bi data/input/document.csv`; при необходимости задайте
  `--encoding`.
- Для TSV явно задайте `--sep '\t'`. Автоопределение поддерживает табуляцию и
  точку с запятой, но явный флаг работает надёжнее.

## 3. Сбои API

- Изучите логи с `event=api_request` и `status_code`. Постоянные 404 означают
  некорректный идентификатор или устаревший словарь.
- При 429/500 отслеживайте счётчик `retry_attempt` — при необходимости уменьшите
  `batch_size` или настройте `sources.*.rps` и `backoff_factor`.
- Для PubChem убедитесь, что `sources.pubchem.user_agent` содержит валидный
  e-mail (требование API).

## 4. Отсутствующие словари

- Таргеты используют CSV из `config/dictionary/_target`. Запустите
  `python scripts/get_target_data.py all --print-config`, чтобы проверить пути.
- Для тестовых объектов убедитесь, что есть файл `molecule_hierarchy.csv`. При
  необходимости пересоберите кэш: `python -m library.integration.molecule_catalog --help`.

## 5. Action type и границы

- Если `action_type=unknown`, проверьте, есть ли метрика в `activity_enrichment.action_type.metrics`
  и разрешена ли она в allowlist.
- Отрицательные `standard_value` сигнализируют о некорректном `relation`; модуль
  `activity_bounds` по умолчанию обрезает значения до нуля.

## 6. Недетерминированный результат

- Используйте одну и ту же версию Python и `requirements-lock.txt`.
- Очистите кэши (`rm -rf $CHEMBL_DA_BASE_PATH/cache`) и повторите запуск.
- Примените `check-determinism --explain`, чтобы увидеть различающиеся файлы.

## 7. Советы для CI

- В контейнерах задайте `CHEMBL_DA_BASE_PATH=/tmp/chembl-da`, чтобы избежать проблем
  с правами.
- Публикуйте `reports/test_report.json` и `reports/test_summary.md` как артефакты.

При повторяющихся ошибках соберите логи (`logs/*.log` по умолчанию или
`<base>/logs` при установленной `CHEMBL_DA_BASE_PATH`), файл метаданных и
проблемные строки CSV, затем создайте issue с описанием уже проделанных шагов.
