# Схема данных

## Входные таблицы

### activity.csv
- **Назначение:** исходный список идентификаторов активностей ChEMBL, который читает CLI `get_activity_data` перед обращением к API.【F:config.yaml†L25-L31】【356235†L1-L5】

| Колонка | Тип данных | Источник | Описание |
| --- | --- | --- | --- |
| activity_chembl_id | строка | Внешний CSV с идентификаторами | Идентификатор активности в базе ChEMBL, используется для поиска записей `/activity`.| 

### assay.csv
- **Назначение:** перечень идентификаторов биологических тестов ChEMBL, которые обрабатываются скриптом `get_assay_data`.【F:config.yaml†L32-L36】【29bb29†L1-L5】

| Колонка | Тип данных | Источник | Описание |
| --- | --- | --- | --- |
| assay_chembl_id | строка | Внешний CSV с идентификаторами | ChEMBL ID теста, ключ для запроса `/assay`.|

### documents.csv
- **Назначение:** список идентификаторов документов ChEMBL, являющийся входом для `get_document_data` (режимы `chembl`/`all`).【F:config.yaml†L49-L56】【fca225†L1-L5】

| Колонка | Тип данных | Источник | Описание |
| --- | --- | --- | --- |
| document_chembl_id | строка | Внешний CSV с идентификаторами | Основной идентификатор публикации в ChEMBL.| 

### targets.csv
- **Назначение:** входной набор идентификаторов целей ChEMBL для пайплайна `get_target_data`.【F:config.yaml†L62-L80】【086d4d†L1-L5】

| Колонка | Тип данных | Источник | Описание |
| --- | --- | --- | --- |
| target_chembl_id | строка | Внешний CSV с идентификаторами | ChEMBL ID белковой мишени, используется при извлечении `/target`.|

### testitem.csv
- **Назначение:** перечень идентификаторов молекул ChEMBL, которые обогащаются структурной и PubChem информацией в `get_testitem_data`.【F:config.yaml†L37-L41】【85074e†L1-L5】

| Колонка | Тип данных | Источник | Описание |
| --- | --- | --- | --- |
| molecule_chembl_id | строка | Внешний CSV с идентификаторами | ChEMBL ID молекулы, служит ключом для `/molecule`.|

## Выходные таблицы

### activity.csv (обработанный экспорт)
- **Назначение:** нормализованный набор экспериментальных значений активности, полученный через ChEMBL API и дополненный метаданными пайплайна.【F:scripts/get_activity_data.py†L63-L220】【F:library/chembl_assay.py†L62-L111】【F:schemas/activities.py†L16-L56】【F:library/pipeline_metadata.py†L60-L84】

| Колонка | Тип данных | Источник | Описание |
| --- | --- | --- | --- |
| activity_id | строка/число | ChEMBL `/activity` | Внутренний идентификатор измерения активности, дублирует `ACTIVITY_ID` из API.|
| molecule_chembl_id | строка | ChEMBL `/activity` | Идентификатор молекулы, для которой снята активность.|
| assay_chembl_id | строка | ChEMBL `/activity` | Ссылка на тест (`assay_chembl_id`), в котором измерена активность.|
| activity_comment | строка | ChEMBL `/activity` | Текстовый комментарий к эксперименту или результату.|
| assay_description | строка | ChEMBL `/activity` | Описание связанного теста.|
| assay_variant_accession | строка | ChEMBL `/activity` | UniProt-акцессия варианта белка, если указана в записи активности.|
| assay_variant_mutation | строка | ChEMBL `/activity` | Информация о мутации белка в тесте.|
| bao_format | строка | ChEMBL `/activity` | Идентификатор формата BAO (BioAssay Ontology).|
| bao_label | строка | ChEMBL `/activity` | Метка BAO, описывающая тип результата.|
| data_validity_comment | строка | ChEMBL `/activity` | Комментарий о валидности данных (например, предупреждения).|
| data_validity_description | строка | ChEMBL `/activity` | Подробное описание статуса валидности.|
| document_chembl_id | строка | ChEMBL `/activity` | Ссылка на исходный документ-публикацию.|
| pchembl_value | строка/число | ChEMBL `/activity` | Логарифмическое значение pChEMBL, если рассчитано источником.|
| potential_duplicate | булево/строка | ChEMBL `/activity` | Флаг возможного дублирования записи.|
| qudt_units | строка | ChEMBL `/activity` | Единицы измерения в терминах QUDT (если доступны).|
| record_id | строка/число | ChEMBL `/activity` | Идентификатор записи в исходной базе.|
| relation | строка | ChEMBL `/activity` | Оператор отношения (`<`, `<=`, `=`, `>=` и т. п.) для исходного значения.|
| src_assay_id | строка/число | ChEMBL `/activity` | Внутренний идентификатор теста в источнике.|
| src_id | строка/число | ChEMBL `/activity` | Идентификатор источника данных (база/партнёр).|
| standard_relation | строка | ChEMBL `/activity` | Нормализованный оператор отношения после обработки.|
| standard_units | строка | ChEMBL `/activity` | Нормализованные единицы измерения.|
| type | строка | ChEMBL `/activity` | Тип исходного измерения (как в API).|
| units | строка | ChEMBL `/activity` | Единицы исходного значения `value`.|
| value | строка/число | ChEMBL `/activity` | Первичное численное значение, как возвращает API.|
| standard_type | строка | ChEMBL `/activity` | Нормализованный тип (ограничен значениями `IC50` или `Ki`).|
| standard_value | число (float) | ChEMBL `/activity` | Нормализованное числовое значение (молярные единицы), гарантированно неотрицательное.|
| pipeline_version | строка | Пайплайн | Версия пакета `chembl-data-acquisition`, добавляется при экспорте.|
| timestamp_utc | строка (ISO 8601) | Пайплайн | Метка времени выполнения экспорта в UTC.|

### assay.csv (обработанный экспорт)
- **Назначение:** агрегированный набор описаний биологических тестов с подсчётом числа тестов на ту же пару документ/мишень и метаданными пайплайна.【F:scripts/get_assay_data.py†L47-L167】【F:library/chembl_assay.py†L22-L59】【F:library/assay_postprocessing.py†L1-L41】【F:schemas/assays.py†L1-L85】【F:library/pipeline_metadata.py†L60-L84】

| Колонка | Тип данных | Источник | Описание |
| --- | --- | --- | --- |
| assay_chembl_id | строка | ChEMBL `/assay` | Основной идентификатор теста.|
| ASSAY_ID | строка | ChEMBL `/assay` | Дополнительный идентификатор теста из первичного источника.|
| Target TYPE | строка | ChEMBL `/assay` | Тип биологической мишени (клеточный, белковый и т. д.).|
| accession | строка | ChEMBL `/assay` | UniProt-акцессия мишени, если известна.|
| acts_per_assay_step5 | строка/число | ChEMBL `/assay` | Количество активностей, зафиксированных на шаге 5 (тип поля сохраняется как в API).|
| assay_cell_type | строка | ChEMBL `/assay` | Тип клетки, использованной в эксперименте.|
| assay_subcellular_fraction | строка | ChEMBL `/assay` | Фракция клетки (например, мембрана), на которой проводился тест.|
| assay_tissue | строка | ChEMBL `/assay` | Ткань или орган, откуда взят материал.|
| bao_format | строка | ChEMBL `/assay` | Идентификатор формата BAO.|
| cited_assay_corr | строка/булево | ChEMBL `/assay` | Флаг, указывает цитируется ли тест как коррелированный.|
| description | строка | ChEMBL `/assay` | Текстовое описание теста.|
| document_chembl_id | строка | ChEMBL `/assay` | Ссылка на документ-источник.|
| error_assay_corr | строка/булево | ChEMBL `/assay` | Признак ошибки в пометке корреляции.|
| higly_correlated_cit | строка/булево | ChEMBL `/assay` | Флаг высокой коррелированности по цитированию.|
| isoform | строка | ChEMBL `/assay` | Номер/описание белкового изоформа.|
| month | строка/число | ChEMBL `/assay` | Месяц публикации, как указан в исходных данных.|
| mutation | строка | ChEMBL `/assay` | Описание мутации мишени.|
| shuffled_cit | строка/булево | ChEMBL `/assay` | Индикатор «перемешанной» цитаты.|
| shuffled_target_assay | строка/булево | ChEMBL `/assay` | Индикатор перемешанной пары мишень/тест.|
| substrate_name | строка | ChEMBL `/assay` | Название субстрата.|
| target_chembl_id | строка | ChEMBL `/assay` | Ссылка на мишень ChEMBL.|
| target_name | строка | ChEMBL `/assay` | Человеко-читаемое имя мишени.|
| version | строка/число | ChEMBL `/assay` | Версия записи теста.|
| year | строка/число | ChEMBL `/assay` | Год публикации.|
| assay_with_same_target | число (int) | Произведено `postprocess_assays` | Количество тестов с той же комбинацией `document_chembl_id` и `target_chembl_id`.|
| pipeline_version | строка | Пайплайн | Версия пакета `chembl-data-acquisition`.|
| timestamp_utc | строка (ISO 8601) | Пайплайн | Метка времени выгрузки.|

### documents.csv (обработанный экспорт)
- **Назначение:** консолидированное библиографическое описание публикаций с объединением данных ChEMBL, PubMed, Semantic Scholar, OpenAlex и Crossref, а также вычисленными классами публикаций.【F:scripts/get_document_data.py†L18-L884】【F:library/document_pipeline.py†L20-L119】【F:schemas/documents.py†L14-L119】【F:library/document_postprocessing.py†L18-L154】【F:library/pipeline_metadata.py†L60-L84】

| Колонка | Тип данных | Источник | Описание |
| --- | --- | --- | --- |
| document_chembl_id | строка | ChEMBL `/document` | Основной идентификатор публикации.|
| title | строка | ChEMBL `/document` | Название публикации.|
| abstract | строка | ChEMBL `/document` | Аннотация.|
| doi | строка | ChEMBL `/document` и агрегаторы | DOI в исходном формате (позже нормализуется).|
| year | число/строка | ChEMBL `/document` | Год публикации.|
| journal | строка | ChEMBL `/document` | Полное название журнала.|
| journal_abbrev | строка | ChEMBL `/document` | Сокращение журнала.|
| volume | число/строка | ChEMBL `/document` | Номер тома.|
| issue | число/строка | ChEMBL `/document` | Номер выпуска.|
| first_page | число/строка | ChEMBL `/document` | Первая страница.|
| last_page | число/строка | ChEMBL `/document` | Последняя страница.|
| pubmed_id | число/строка | ChEMBL `/document` | PubMed ID, если ChEMBL его содержит.|
| authors | строка | ChEMBL `/document` | Список авторов.|
| source | строка | ChEMBL `/document` | Источник записи (обычно `ChEMBL`).|
| doi_normalised | строка | Производное поле | DOI приведённый к каноническому виду.|
| publication_types_normalised | строка | Производное поле | Объединённые типы публикаций из всех источников.|
| publication_type_score_review | число (int) | Производное поле | Балл принадлежности к обзорам (scoring).|
| publication_type_score_experimental | число (int) | Производное поле | Балл экспериментальной публикации.|
| publication_type_score_unknown | число (int) | Производное поле | Балл неопределённого типа.|
| publication_class | строка | Производное поле | Итоговая классификация (например, `review`/`experimental`).|
| PubMed.PMID | число/строка | PubMed API | Основной идентификатор PubMed.|
| PubMed.DOI | строка | PubMed API | DOI из PubMed после нормализации.|
| PubMed.ArticleTitle | строка | PubMed API | Название статьи из PubMed.|
| PubMed.Abstract | строка | PubMed API | Аннотация из PubMed.|
| PubMed.JournalTitle | строка | PubMed API | Название журнала в PubMed.|
| PubMed.Volume | число/строка | PubMed API | Номер тома из PubMed.|
| PubMed.Issue | число/строка | PubMed API | Номер выпуска из PubMed.|
| PubMed.StartPage | число/строка | PubMed API | Стартовая страница.|
| PubMed.EndPage | число/строка | PubMed API | Конечная страница.|
| PubMed.PublicationType | строка | PubMed API | Сырые термины типов публикаций.|
| PubMed.MeSH_Descriptors | строка | PubMed API | Перечень MeSH-дескрипторов.|
| PubMed.MeSH_Qualifiers | строка | PubMed API | Перечень MeSH-квалификаторов.|
| PubMed.ChemicalList | строка | PubMed API | Перечень химических сущностей.|
| PubMed.DayRevised | число/строка | PubMed API | День последней правки записи.|
| PubMed.MonthRevised | число/строка | PubMed API | Месяц последней правки.|
| PubMed.YearRevised | число/строка | PubMed API | Год последней правки.|
| PubMed.YearCompleted | число/строка | PubMed API | Год завершения индексации.|
| PubMed.MonthCompleted | число/строка | PubMed API | Месяц завершения индексации.|
| PubMed.DayCompleted | число/строка | PubMed API | День завершения индексации.|
| PubMed.Error | строка | PubMed API | Сообщение об ошибке запроса (если было).|
| PubMed.ISSN | строка | PubMed API | ISSN журнала.|
| scholar.PMID | число/строка | Semantic Scholar | Привязанный PubMed ID.|
| scholar.Venue | строка | Semantic Scholar | Площадка/журнал в Semantic Scholar.|
| scholar.PublicationTypes | строка | Semantic Scholar | Типы публикаций.|
| scholar.SemanticScholarId | строка | Semantic Scholar | Внутренний ID Semantic Scholar.|
| scholar.ExternalIds | строка | Semantic Scholar | Внешние идентификаторы (JSON).|
| scholar.DOI | строка | Semantic Scholar | DOI из Semantic Scholar.|
| scholar.Error | строка | Semantic Scholar | Диагностика ошибок вызова.|
| OpenAlex.PublicationTypes | строка | OpenAlex | Типы публикаций (список).|
| OpenAlex.TypeCrossref | строка | OpenAlex/Crossref | Классификация Crossref, возвращаемая OpenAlex.|
| OpenAlex.Genre | строка | OpenAlex | Жанр публикации.|
| OpenAlex.Id | строка | OpenAlex | Идентификатор записи OpenAlex.|
| OpenAlex.Venue | строка | OpenAlex | Площадка публикации по данным OpenAlex.|
| OpenAlex.MeshDescriptors | строка | OpenAlex | MeSH-дескрипторы от OpenAlex.|
| OpenAlex.MeshQualifiers | строка | OpenAlex | MeSH-квалификаторы от OpenAlex.|
| OpenAlex.Error | строка | OpenAlex | Сообщение об ошибке.|
| crossref.Type | строка | Crossref | Тип публикации из Crossref.|
| crossref.Subtype | строка | Crossref | Подтип публикации.|
| crossref.Title | строка | Crossref | Заголовок по Crossref.|
| crossref.Subtitle | строка | Crossref | Подзаголовок.|
| crossref.Subject | строка | Crossref | Тематические рубрики.|
| crossref.Error | строка | Crossref | Сообщение об ошибке запроса.|
| date_code | строка | Постобработка | Служебный код даты (формируется при нормализации).|
| Index | число/строка | Постобработка | Порядковый индекс из исходной таблицы.|
| PubMed.is_review | булево/строка | Постобработка | Флаг «обзор» на основе PubMed терминов.|
| scholar.is_review | булево/строка | Постобработка | Флаг «обзор» из Semantic Scholar.|
| OpenAlex.is_review | булево/строка | Постобработка | Флаг «обзор» по данным OpenAlex.|
| pipeline_version | строка | Пайплайн | Версия пакета `chembl-data-acquisition`.|
| timestamp_utc | строка (ISO 8601) | Пайплайн | Время формирования выгрузки.|

### target.csv (финализированный экспорт)
- **Назначение:** унифицированная таблица белковых мишеней с объединением атрибутов ChEMBL, UniProt и IUPHAR, приведённая к детерминированному порядку колонок и форматам, совместимым с исходной BI-процессингом.【F:scripts/get_target_data.py†L31-L1148】【F:library/chembl_target.py†L185-L315】【F:library/target_postprocessing.py†L181-L442】【F:schemas/targets.py†L17-L214】【F:library/pipeline_metadata.py†L60-L84】

| Колонка | Тип данных | Источник | Описание |
| --- | --- | --- | --- |
| target_chembl_id | строка | ChEMBL `/target` | Основной идентификатор мишени.|
| uniprot_id_primary | строка | Постобработка UniProt | Основная UniProt-акцессия (нормализованный `uniProtkbId`).|
| uniprot_ids_all | строка | Постобработка UniProt | Слитый список всех доступных UniProt идентификаторов (основной, вторичные, карта ID).|
| isoform_ids | строка | UniProt (через merge) | Список идентификаторов изоформ, если доступны.|
| isoform_names | строка | UniProt (через merge) | Названия изоформ.|
| isoform_synonyms | строка | UniProt (через merge) | Синонимы изоформ (нижний регистр).|
| hgnc_id | строка | HGNC (через cross-references) | Идентификатор HGNC без префикса.|
| gene_symbol | строка | Постобработка | Итоговый символ гена (верхний регистр).|
| protein_name_canonical | строка | UniProt/ChEMBL | Каноническое название белка.|
| protein_name_alt | строка | UniProt/ChEMBL | Альтернативные названия (pipe-разделённый список).|
| organism | строка | ChEMBL/UniProt | Род (genus) организма, для которого описана мишень.|
| taxon_id | строка | ChEMBL | Идентификатор таксона NCBI.|
| lineage_superkingdom | строка | UniProt taxon merge | Надцарство организма.|
| lineage_phylum | строка | UniProt taxon merge | Тип организма.|
| lineage_class | строка | UniProt taxon merge | Класс организма.|
| sequence_length | строка | UniProt | Длина аминокислотной последовательности.|
| features_signal_peptide | строка | UniProt features | Наличие сигнального пептида (нижний регистр).|
| features_transmembrane | строка | UniProt features | Признак трансмембранных сегментов (нижний регистр).|
| features_topology | строка | UniProt features | Топология (например, внутриклеточная/внеклеточная).|
| ptm_glycosylation | строка | UniProt PTM | Аннотации гликозилирования.|
| ptm_lipidation | строка | UniProt PTM | Аннотации липидирования.|
| ptm_disulfide_bond | строка | UniProt PTM | Аннотации дисульфидных мостиков.|
| ptm_modified_residue | строка | UniProt PTM | Прочие модифицированные остатки.|
| xref_chembl | строка | Постобработка | Повторяет `target_chembl_id` для совместимости.|
| xref_uniprot | строка | Постобработка | Повторяет `uniprot_id_primary`.|
| xref_ensembl | строка | ChEMBL/UniProt | Список Ensembl ID из cross-references.|
| xref_iuphar | строка | IUPHAR | Идентификатор мишени в IUPHAR.|
| gtop_target_id | строка | Guide to PHARMACOLOGY | Идентификатор объекта в базе GToP.|
| gtop_synonyms | строка | Guide to PHARMACOLOGY | Синонимы из GToP (нижний регистр).|
| gtop_natural_ligands_n | строка | Guide to PHARMACOLOGY | Количество естественных лигандов (как в источнике).|
| gtop_interactions_n | строка | Guide to PHARMACOLOGY | Количество взаимодействий (как в источнике).|
| gtop_function_text_short | строка | Guide to PHARMACOLOGY | Краткое описание функции.|
| uniprot_last_update | строка | UniProt | Дата последнего обновления записи.|
| uniprot_version | строка | UniProt | Версия UniProt записи.|
| pipeline_version | строка | Пайплайн | Версия пакета `chembl-data-acquisition`.|
| timestamp_utc | строка (ISO 8601) | Пайплайн | Время формирования выгрузки.|
| pfam | строка | UniProt | Список доменов Pfam из UniProt.|
| interpro | строка | UniProt | Список доменов InterPro из UniProt.|
| xref_pdb | строка | UniProt | Ссылки на структуры PDB.|
| xref_alphafold | строка | UniProt/AlphaFold | Ссылки на модели AlphaFold.|
| hgnc_name | строка | HGNC | Полное имя гена в HGNC.|
| uniProtkbId | строка | UniProt | Исходное поле UniProt ID до нормализации.|
| secondaryAccessions | строка | UniProt | Список вторичных акцессий UniProt.|
| recommendedName | строка | UniProt | Рекомендованное название белка (как в UniProt).|
| geneName | строка | UniProt | Основной символ гена из UniProt.|
| secondaryAccessionNames | строка | UniProt | Названия, связанные со вторичными акцессиями.|
| molecular_function | строка | UniProt GO | GO-аннотации молекулярной функции.|
| cellular_component | строка | UniProt GO | GO-аннотации клеточной локализации.|
| subcellular_location | строка | UniProt | Локализация белка.|
| topology | строка | UniProt | Топология трансмембранного белка.|
| transmembrane | строка | UniProt | Индикатор трансмембранных сегментов (после нормализации).|
| intramembrane | строка | UniProt | Индикатор внутримембранных сегментов.|
| glycosylation | строка | UniProt | Аннотация гликозилирования (нижний регистр).|
| lipidation | строка | UniProt | Аннотация липидирования.|
| disulfide_bond | строка | UniProt | Аннотация дисульфидных связей.|
| modified_residue | строка | UniProt | Прочие модификации остатков.|
| phosphorylation | строка | UniProt | Аннотации фосфорилирования.|
| acetylation | строка | UniProt | Аннотации ацетилирования.|
| ubiquitination | строка | UniProt | Аннотации убиквитинирования.|
| signal_peptide | строка | UniProt | Аннотации сигнального пептида (нижний регистр).|
| propeptide | строка | UniProt | Аннотации пропептида.|
| GuidetoPHARMACOLOGY | строка | IUPHAR/GToP | Основной ID в Guide to PHARMACOLOGY.|
| family | строка | IUPHAR | Название семейства по IUPHAR.|
| SUPFAM | строка | UniProt/IUPHAR | Суперсемейство (как в источнике).|
| PROSITE | строка | UniProt | Доменные мотивы PROSITE.|
| InterPro | строка | UniProt | Домены InterPro (альтернативный источник).|
| Pfam | строка | UniProt | Домены Pfam (альтернативное поле).|
| PRINTS | строка | UniProt | Мотивы PRINTS.|
| TCDB | строка | UniProt | Классификация транспортеров TCDB.|
| pref_name | строка | ChEMBL | Преферируемое имя мишени в ChEMBL.|
| target_type | строка | ChEMBL | Тип мишени по ChEMBL (`SINGLE PROTEIN`, `PROTEIN FAMILY` и т. п.).|
| tax_id | строка | ChEMBL | Значение `tax_id` из ChEMBL (строковое представление).|
| species_group_flag | строка | ChEMBL | Флаг группировки видов.|
| target_components | строка (JSON) | ChEMBL | Сериализованное описание компонент мишени.|
| protein_classifications | строка (JSON) | ChEMBL | Иерархия белковых классов.|
| cross_references | строка (JSON) | ChEMBL | Внешние ссылки из записи ChEMBL.|
| gene_symbol_list | строка | Постобработка | Объединённый список символов гена (pipe-разделитель, верхний регистр).|
| protein_synonym_list | строка | Постобработка | Сводный список синонимов белка (нижний регистр).|
| reactions | строка | ChEMBL | Связанные реакции (если присутствуют).|
| reaction_ec_numbers | строка | Постобработка | Список EC-номеров реакций, собранный из синонимов и xref.|
| protein_class_pred_L1 | строка | IUPHAR доп.данные | Предсказанный класс белка (уровень 1).|
| protein_class_pred_L2 | строка | IUPHAR доп.данные | Предсказанный класс (уровень 2).|
| protein_class_pred_L3 | строка | IUPHAR доп.данные | Предсказанный класс (уровень 3).|
| protein_class_pred_rule_id | строка | IUPHAR доп.данные | Идентификатор правила классификации.|
| protein_class_pred_evidence | строка | IUPHAR доп.данные | Тип доказательств классификации.|
| protein_class_pred_confidence | строка | IUPHAR доп.данные | Оценка уверенности классификации.|
| iuphar_target_id | строка | IUPHAR | Идентификатор мишени в IUPHAR.|
| iuphar_family_id | строка | IUPHAR | Идентификатор семейства в IUPHAR.|
| iuphar_type | строка | IUPHAR | Тип объекта по IUPHAR.|
| iuphar_class | строка | IUPHAR | Класс объекта.|
| iuphar_subclass | строка | IUPHAR | Подкласс объекта.|
| iuphar_chain | строка | IUPHAR | Идентификатор цепи/подъединицы.|
| iuphar_name | строка | IUPHAR | Официальное название в IUPHAR.|
| iuphar_full_id_path | строка | IUPHAR | Полный путь идентификаторов в иерархии.|
| iuphar_full_name_path | строка | IUPHAR | Полный путь названий в иерархии.|

### testitem.csv (обработанный экспорт)
- **Назначение:** описание лекарственных веществ/соединений из ChEMBL, дополненное структурными полями и обогащением из PubChem, нормализованное и снабжённое метаданными пайплайна.【F:scripts/get_testitem_data.py†L36-L114】【F:library/chembl_assay.py†L91-L111】【F:schemas/testitems.py†L12-L37】【F:library/pipeline_metadata.py†L60-L84】

| Колонка | Тип данных | Источник | Описание |
| --- | --- | --- | --- |
| molecule_chembl_id | строка | ChEMBL `/molecule` | Основной идентификатор молекулы.|
| black_box_warning | строка/булево | ChEMBL `/molecule` | Флаг наличия «black box» предупреждения.|
| first_approval | строка/дата | ChEMBL `/molecule` | Год/дата первого одобрения препарата (как в источнике).|
| max_phase | строка | ChEMBL `/molecule` | Максимальная стадия клинических испытаний.|
| canonical_smiles | строка | ChEMBL `/molecule` | Канонический SMILES из ChEMBL.|
| standard_inchi | строка | ChEMBL `/molecule` | Стандартный InChI.|
| standard_inchi_key | строка | ChEMBL `/molecule` | Ключ InChIKey.|
| molecule_type | строка | ChEMBL `/molecule` | Тип молекулы (small molecule, antibody и т. п.).|
| oral | строка/булево | ChEMBL `/molecule` | Признак оральной формы.|
| parenteral | строка/булево | ChEMBL `/molecule` | Признак парентерального применения.|
| pref_name | строка | ChEMBL `/molecule` | Преферируемое название.|
| pubchem_canonical_smiles | строка | PubChem | Канонический SMILES, полученный по CID.|
| pubchem_cid | строка/число | PubChem | Идентификатор PubChem CID для SMILES.|
| pubchem_inchi | строка | PubChem | Стандартный InChI из PubChem.|
| pubchem_inchikey | строка | PubChem | InChIKey из PubChem.|
| pubchem_isomeric_smiles | строка | PubChem | Изомерический SMILES из PubChem.|
| pubchem_iupac_name | строка | PubChem | IUPAC-название из PubChem.|
| pubchem_molecular_formula | строка | PubChem | Молекулярная формула из PubChem.|
| structure_type | строка | ChEMBL `/molecule` | Тип структуры (например, `MOLFILE`).|
| topical | строка/булево | ChEMBL `/molecule` | Признак топического применения.|
| pipeline_version | строка | Пайплайн | Версия пакета `chembl-data-acquisition`.|
| timestamp_utc | строка (ISO 8601) | Пайплайн | Время формирования выгрузки.|

