# Activity flags provenance — 20251005

Ниже представлена детализация происхождения и вычисления флагов финальной таблицы `extended.output.activity_20251005.csv`. Для каждого признака приведены источники данных, порядок правил, соответствие историческому Power Query и фактическая Python-реализация.

## unknown_chirality

| Поле | Где создаётся/меняется | Источники (таблицы/колонки) | Алгоритм/правило (пошагово) | Приоритет/порядок правил | Тип и допуски | Значения по умолчанию / NULL | Детерминизм | M-эквивалент | Python-реализация |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| unknown_chirality | `library/postprocessing/activity_extended.py:L497-L504`; объяснение `library/postprocessing/_explain_activity_flags.py:L95-L152` | Колонка `nstereo` из входного `output.activity_*.csv` | 1. Привести `nstereo` к `Int64` (`_safe_to_int`).<br>2. Флаг = `nstereo != 1`.<br>3. Пустые значения → `True`.<br>4. Колонка удаляется из набора, остаётся только флаг. | Сначала используется фактическое значение `nstereo`; при его отсутствии подставляется `True`. | `boolean` (`pd.BooleanDtype`) | Если `nstereo` отсутствует или `NA`, выставляется `True`. | Вход сортируется `helpers.sort_power_query`, затем фиксировано преобразуется; отсутствие случайности. | Шаг Power Query `Changed Type` и `Filled Down` в историческом листе: после объединения молекул столбец `unknown_chirality` вычислялся из `NSTEREO` и заполнялся вниз. | `_prepare_unknown_chirality` + `explain_unknown_chirality` фиксируют логику и комментарии. |

**Пример данных**

| activity_chembl_id | raw_nstereo | unknown_chirality |
| --- | --- | --- |
| ACT-101 | 1 | False |
| ACT-102 | 1 | False |
| ACT-103 | 2 | True |
| ACT-201 | — | True |
| ACT-301 | 0 | True |

**Краевые случаи**
- `nstereo` = 1 ⇒ флаг строго `False`.
- Любое другое числовое значение (`0`, `2`, `-1`) ⇒ `True`.
- Полностью отсутствующая колонка `nstereo` ⇒ создаётся `True` для всех записей.
- Строковые представления (`"1"`, `"0"`) допустимы благодаря `_safe_to_int`.

**Unit-тест кейсы**
- `tests/postprocessing/test_activity_flags_provenance.py::test_explain_unknown_chirality__cases`
- Косвенно: `tests/postprocessing/test_activity_extended.py::test_transform_activity_frame__parses_activity_properties_flags`

## multmol_assay

| Поле | Где создаётся/меняется | Источники | Алгоритм/правило | Приоритет/порядок правил | Тип и допуски | Значения по умолчанию / NULL | Детерминизм | M-эквивалент | Python-реализация |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| multmol_assay | Первичная загрузка из CSV, затем `_apply_multimol_logic` (`activity_extended.py:L507-L535`); объяснение `*_explain_activity_flags.py:L155-L220` | Колонка `multmol_assay` входного файла; комбинация ключей (`salt_chembl_id`, `molecule_chembl_id`, `target_chembl_id`, `assay_chembl_id`, `standard_type`); `unknown_chirality` | 1. Считываем исходный флаг (`_safe_to_bool`).<br>2. Считаем размер групп по ключам.<br>3. Если `unknown_chirality=False`, исходный флаг пуст и размер группы > 1, все записи этой группы помечаются как True (мульти-молекулярный тест).<br>4. Итог = логическое «ИЛИ» исходного и вычисленного флага. | Приоритет у ручного значения: если пользователь указал True, оно сохраняется. Затем применяется групповой детектор. | `boolean` | Пустые значения → `False`; при срабатывании группового правила становятся `True`. | Группы сортируются `helpers.sort_power_query`; `drop_duplicates` на финальном этапе фиксирует порядок. | Power Query шаги `Group Rows` и `Custom` в legacy workbook объединяли активности по ключам и ставили флаг «MULTMOL_ASSAY_SAME». | `_apply_multimol_logic` и `explain_multmol_assay` документируют групповую логику и причину. |

**Пример данных**

| activity_chembl_id | raw_multmol_assay | unknown_chirality | финальный multmol_assay |
| --- | --- | --- | --- |
| ACT-101 | — | False | True (дубликаты ASSAY-X) |
| ACT-102 | 1 | False | True (исходное значение) |
| ACT-103 | 0 | True | False |
| ACT-201 | 0 | True/NA | False |

**Краевые случаи**
- Запись с `unknown_chirality=True` никогда не поднимает флаг через групповое правило.
- Если группа >1, но исходный флаг явно `False`, итог станет `True` (привязка ко всей группе).
- Несуществующие ключевые колонки → `ActivityExtendedError` (гарантия целостности).

**Unit-тест кейсы**
- `tests/postprocessing/test_activity_flags_provenance.py::test_explain_multmol_assay__duplicate_trigger`
- `tests/postprocessing/test_activity_extended.py::test_apply_multimol_logic__marks_duplicate_multimol_assay`

## exact_data_citation

| Поле | Где создаётся/меняется | Источники | Алгоритм/правило | Приоритет/порядок правил | Тип | Значения по умолчанию / NULL | Детерминизм | M-эквивалент | Python-реализация |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| exact_data_citation | Переименование (`activity_extended.py:L577-L594`), нормализация `_compute_citation_flags` (`L740-L755`); объяснение `*_explain_activity_flags.py:L223-L234` | Колонки `exact_cited_activity` и (если присутствует) `exact_data_citation` | 1. Переименовать `exact_cited_activity` в целевой столбец.<br>2. Привести к `boolean` (`_safe_to_bool`) и заполнить `False` для `NA`.<br>3. В `is_citation` участвует при расчёте частоты цитирования. | Если присутствуют оба столбца (`exact_cited_activity` и уже подготовленный `exact_data_citation`), используется первый ненулевой. | `boolean` | `NA` → `False`. Невалидные значения вызывают `ValueError` и логируются. | Те же операции выполняются для каждой строки; без случайности. | Power Query шаг `Changed Type` + расчёт признаков цитирования (флаг `exact_cited_activity`). | `_compute_citation_flags` + `explain_exact_data_citation`. |

**Пример данных**

| activity_chembl_id | raw_exact_cited_activity | exact_data_citation |
| --- | --- | --- |
| ACT-101 | 1 | True |
| ACT-102 | 0 | False |
| ACT-103 | 0 | False |
| ACT-201 | 1 | True |

**Краевые случаи**
- Строковые значения (`"TRUE"`, `"false"`) корректно приводятся.
- Если столбец отсутствует, создаётся `False` (и `is_citation` учитывает только другие признаки).

**Unit-тест кейсы**
- `tests/postprocessing/test_activity_flags_provenance.py::test_boolean_flag_explainers__source_priority`
- Косвенно: `tests/postprocessing/test_activity_extended.py::test_transform_activity_frame__parses_activity_properties_flags`

## higly_correlated_assay / highly_correlated_assay

| Поле | Где создаётся/меняется | Источники | Алгоритм/правило | Приоритет | Тип | Значения по умолчанию / NULL | Детерминизм | M-эквивалент | Python-реализация |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| higly_correlated_assay | Переименование `higly_correlated_cit` → `higly_correlated_assay` (`activity_extended.py:L577-L594`), нормализация `_compute_citation_flags`; alias создаётся там же | `higly_correlated_cit`, `higly_correlated_assay`, `highly_correlated_assay` (если уже присутствует) | 1. Переименовать входную колонку.<br>2. Создать alias `highly_correlated_assay` с теми же значениями (`activity_extended.py:L589-L593`).<br>3. Привести к `boolean`, `NA` → `False`. | Источник приоритизируется в порядке: `higly_correlated_cit` → `higly_correlated_assay` → `highly_correlated_assay`. Alias всегда копия целевого флага. | `boolean` | `NA` → `False`. | Инвариантное поведение через `helpers.sort_power_query`; alias синхронизирован детерминированно. | Power Query столбец `higly_correlated_cit` (орфография наследована), который далее использовался в QA-классификаторе. | `_compute_citation_flags`, `_rename_columns`; объяснение `*_explain_activity_flags.py:L237-L255`. |

**Пример данных**

| activity_chembl_id | raw_higly_correlated_cit | higly_correlated_assay | highly_correlated_assay |
| --- | --- | --- | --- |
| ACT-102 | 1 | True | True |
| ACT-103 | 0 | False | False |
| ACT-201 | 0 | False | False |

**Краевые случаи**
- Alias гарантированно совпадает, тесты контролируют это.
- При отсутствии исходной колонки оба признака заполняются `False`.
- Орфография `higly_` сохранена ради обратной совместимости, но отчёт фиксирует наличие alias.

**Unit-тест кейсы**
- `tests/postprocessing/test_activity_flags_provenance.py::test_boolean_flag_explainers__source_priority`
- `tests/postprocessing/test_activity_extended.py::test_transform_activity_frame__parses_activity_properties_flags`

## shuffled_assay

| Поле | Где создаётся/меняется | Источники | Алгоритм/правило | Приоритет | Тип | Значения по умолчанию / NULL | Детерминизм | M-эквивалент | Python-реализация |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| shuffled_assay | Переименование `shuffled_cit` (`activity_extended.py:L577-L594`), нормализация `_compute_citation_flags`; объяснение `*_explain_activity_flags.py:L259-L270` | `shuffled_cit`, `shuffled_assay` | 1. Переименовать колонку.<br>2. Привести к булевому типу, `NA` → `False`. | В порядке появления: `shuffled_cit` → `shuffled_assay`. | `boolean` | `NA` → `False`. | Детеминированный благодаря фиксированным операциям преобразования и сортировке входа. | Power Query шаг `Changed Type` для столбца `shuffled_cit`. | `_compute_citation_flags`, `explain_shuffled_assay`. |

**Пример данных**

| activity_chembl_id | raw_shuffled_cit | shuffled_assay |
| --- | --- | --- |
| ACT-102 | 1 | True |
| ACT-101 | 0 | False |
| ACT-103 | 0 | False |

**Краевые случаи**
- Пустой столбец => все `False`.
- Невалидные значения фиксируются через исключение, что делает проблему заметной в QA.

**Unit-тест кейсы**
- `tests/postprocessing/test_activity_flags_provenance.py::test_boolean_flag_explainers__source_priority`

## review

| Поле | Где создаётся/меняется | Источники | Алгоритм/правило | Приоритет | Тип | Значения по умолчанию / NULL | Детерминизм | M-эквивалент | Python-реализация |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| review | Переименование `review_doc` (`activity_extended.py:L577-L594`); при мерже документов добавляется `review` из словаря (`L790-L812`); объяснение `*_explain_activity_flags.py:L273-L284` | `review_doc` (из основного экспорта), `review` (из `document.csv`) | 1. При объединении с `document.csv` значения приводятся к `boolean` (`_merge_document_metadata`).<br>2. После переименования `review_doc` → `review` дублирующий столбец из словаря перезаписывается (сохраняется приоритет «сырых» данных).<br>3. `_compute_citation_flags` нормализует тип и заполняет `False`. | Используется `review_doc`, если оно не `NA`; иначе — словарное значение. | `boolean` | При отсутствии обоих источников → `False`. | Детеминизм обеспечивается фиксированным порядком merge и сортировкой `helpers.sort_power_query`. | Power Query шаги `Merged Queries` + `Changed Type` в историческом workbook, где `review_doc` замещал `review`. | `_merge_document_metadata`, `_compute_citation_flags`, `explain_review`. |

**Пример данных**

| activity_chembl_id | raw_review_doc | словарный review | итоговый review |
| --- | --- | --- | --- |
| ACT-301 | 1 | FALSE | True |
| ACT-101 | 0 | TRUE | False |
| ACT-401 | — | FALSE | False |

**Краевые случаи**
- Если документ отсутствует в словаре и `review_doc` пустой, флаг = `False`.
- При конфликте (словари говорят `TRUE`, `review_doc`=0) побеждает `review_doc` (отражено в объяснении).

**Unit-тест кейсы**
- `tests/postprocessing/test_activity_flags_provenance.py::test_boolean_flag_explainers__source_priority`
- `tests/postprocessing/test_activity_extended.py::test_transform_activity_frame__parses_activity_properties_flags`

## rounded_data_citation

| Поле | Где создаётся/меняется | Источники | Алгоритм/правило | Приоритет | Тип | Значения по умолчанию / NULL | Детерминизм | M-эквивалент | Python-реализация |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rounded_data_citation | Переименование `approx_cited_activity` (`activity_extended.py:L577-L594`); `_compute_citation_flags`; объяснение `*_explain_activity_flags.py:L287-L298` | `approx_cited_activity`, `rounded_data_citation` (если уже есть) | 1. Переименовать исходное поле (и удалить дубликаты).<br>2. Привести к булевому типу, заполнить `False`. | Используется первое доступное значение: `approx_cited_activity` → `rounded_data_citation`. | `boolean` | `NA` → `False`. | Детеминизм за счёт фиксированного переименования и сортировки. | Power Query столбец `approx_cited_activity`. | `_compute_citation_flags`, `explain_rounded_data_citation`. |

**Пример данных**

| activity_chembl_id | raw_approx_cited_activity | rounded_data_citation |
| --- | --- | --- |
| ACT-101 | 1 | True |
| ACT-102 | 0 | False |
| ACT-103 | 0 | False |

**Краевые случаи**
- При наличии обоих столбцов (как в CSV) приоритет у `approx_cited_activity`; исторический `rounded_data_citation` используется только если первый пустой.

**Unit-тест кейсы**
- `tests/postprocessing/test_activity_flags_provenance.py::test_boolean_flag_explainers__source_priority`

## high_citation_rate

| Поле | Где создаётся/меняется | Источники | Алгоритм/правило | Приоритет | Тип | Значения по умолчанию / NULL | Детерминизм | M-эквивалент | Python-реализация |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| high_citation_rate | `_annotate_high_citation` (`activity_extended.py:L758-L787`); объяснение `*_explain_activity_flags.py:L301-L363` | `document_chembl_id`, агрегированные `is_citation` (из пяти флагов), `config/dictionary/_activity/citation_fraction.csv` | 1. Для каждого документа подсчитать `n_citation` и `n_non_citation` по признаку `is_citation`.<br>2. Исключить документы без обеих категорий.<br>3. Подтянуть пороги `K_min_significant` по `N = n_citation + n_non_citation`.<br>4. `True`, если порог существует и `n_citation >= K_min_significant`; иначе `False`.<br>5. Логировать результат через `_log_flag_summary`. | Порог из словаря имеет высший приоритет; если строка отсутствует, флаг автоматически `False`. | `boolean` | Документы без смешанных наблюдений или без порога → `False`. | Группировка и merge выполняются в фиксированном порядке; `helpers.sort_power_query` обеспечивает стабильность. | Power Query таблица `citation_fraction` использовалась аналогично для QA (критерий `K_min`). | `_annotate_high_citation`, `explain_high_citation_rate`. |

**Пример данных**

| document_chembl_id | n_citation | n_non_citation | N | threshold | high_citation_rate |
| --- | --- | --- | --- | --- | --- |
| DOC-1 | 2 | 1 | 3 | 2 | True |
| DOC-2 | 2 | 0 | 2 | 1 | False (нет non-citation) |
| DOC-3 | 0 | 1 | 1 | — | False |
| DOC-4 | 0 | 1 | 1 | — | False |

**Краевые случаи**
- Если `N` отсутствует в словаре `citation_fraction`, флаг = `False` и объяснение отражает отсутствие порога.
- Документ с только цитатами либо только нецитатами не считается «high». 
- Словарь порогов контролируется конфигурацией; изменения требуют обновления золотых файлов.

**Unit-тест кейсы**
- `tests/postprocessing/test_activity_flags_provenance.py::test_explain_high_citation_rate__thresholds`
- `tests/postprocessing/test_activity_flags_provenance.py::test_process_activity_extended__golden_snapshot`

## original_activity_approx & original_activity_exact

| Поле | Где создаётся/меняется | Источники | Алгоритм/правило | Приоритет | Тип | Значения по умолчанию / NULL | Детерминизм | M-эквивалент | Python-реализация |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| original_activity_approx / original_activity_exact | Переносятся напрямую из входного CSV (`_ensure_required_input_columns`, `activity_extended.py:L538-L575`), тип не принудительно меняется; объяснение `*_explain_activity_flags.py:L366-L403` | Колонки `original_activity_approx`, `original_activity_exact` | 1. При отсутствии колонок создаются `Float64` серии и сразу же конвертируются в пустые строки (`_ensure_required_input_columns`).<br>2. Значения не нормализуются, только приводятся к `string` при объяснении.<br>3. В финальном экспорте сохраняются как текст (`object`). | Единственный источник — сам CSV; приоритетов нет. | `object` (строковое представление) | `NA` сохраняется (`<NA>`). | Значения не модифицируются, порядок строк фиксируется сортировкой и `dedupe_final`. | Power Query просто копировал текстовые метки (`original_activity_exact` / `approx`). | `_ensure_required_input_columns` (заполнение), `_select_and_cast` (сохранение), `explain_original_activity_flags`. |

**Пример данных**

| activity_chembl_id | raw_original_activity_approx | raw_original_activity_exact | итоговый approx | итоговый exact |
| --- | --- | --- | --- | --- |
| ACT-101 | approx | exact | approx | exact |
| ACT-102 | approx | — | approx | — |
| ACT-103 | — | — | — | — |
| ACT-202 | approx | exact | approx | exact |
| ACT-301 | — | exact | — | exact |

**Краевые случаи**
- Пустые строки/`NA` остаются пустыми; downstream-логика должна учитывать отсутствие данных.
- Не проводится нормализация регистра; значения передаются как есть.

**Unit-тест кейсы**
- `tests/postprocessing/test_activity_flags_provenance.py::test_explain_original_activity_flags__strings`
- `tests/postprocessing/test_activity_flags_provenance.py::test_activity_flags_output_schema__dtypes`

## Дополнительные сведения

- **Логирование метрик.** `_log_flag_summary` (`activity_extended.py:L440-L459`, вызов `L1288-L1293`) записывает в журнал количество `True/False/NA` и доли по каждому полю, что фиксирует детерминизм и облегчает QA.
- **Golden snapshot.** Файл `tests/golden/extended.output.activity_20251005.csv` фиксирует эталонный результат; тест `test_process_activity_extended__golden_snapshot` сравнивает и валидирует SHA-256 (`7ed915aec9d6cbf75585722750474d5b7444f4070ad1721f972d3490e237671b`).
- **Идемпотентность.** Все группировки выполняются через `helpers.sort_power_query`, а финальный порядок строк закреплён `dedupe_final` (`activity_extended.py:L900-L938`).
