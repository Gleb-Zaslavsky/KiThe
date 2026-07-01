//! I/O слой для `TGASeries`.
//!
//! Этот модуль держит сериализацию и десериализацию серии экспериментов
//! в CSV-формате, включая метаданные, привязки колонок и восстановление
//! `TGASeries` из сохранённого файла.

use super::*;
use crate::Kinetics::experimental_kinetics::column_provenance::ColumnProvenance;

#[derive(Clone, Debug, Default)]
struct SeriesExperimentRecord {
    id: String,
    heating_rate: Option<f64>,
    isothermal_temperature: Option<f64>,
    comment: Option<String>,
    len: usize,
    columns: Vec<SeriesColumnRecord>,
    binds: HashMap<String, Option<String>>,
}

#[derive(Clone, Debug)]
struct SeriesColumnRecord {
    column_name: String,
    header_name: String,
    unit: Unit,
    origin: ColumnOrigin,
}

pub(super) fn to_csv_series_impl(series: &TGASeries, path: &Path) -> Result<(), TGADomainError> {
    let mut collected: Vec<DataFrame> = Vec::with_capacity(series.experiments.len());
    let mut max_height = 0usize;
    for exp in &series.experiments {
        let df = exp.dataset.frame.clone().collect()?;
        max_height = max_height.max(df.height());
        collected.push(df);
    }

    let mut meta_lines: Vec<String> = Vec::new();
    let mut export_columns: Vec<Column> = Vec::new();
    let mut header_counts: HashMap<String, usize> = HashMap::new();

    for (idx, exp) in series.experiments.iter().enumerate() {
        let df = &collected[idx];
        let exp_id = exp.meta.id.clone();
        let exp_id_token = escape_meta_field(&exp_id);

        meta_lines.push(format!(
            "#META\tEXP\t{}\t{}\t{}\t{}",
            exp_id_token,
            encode_opt_f64(exp.meta.heating_rate),
            encode_opt_f64(exp.meta.isothermal_temperature),
            encode_opt_string(exp.meta.comment.as_deref()),
        ));
        meta_lines.push(format!("#META\tLEN\t{}\t{}", exp_id_token, df.height()));

        for (bind_key, bind_value) in [
            ("time", exp.dataset.schema.time.as_ref()),
            ("temperature", exp.dataset.schema.temperature.as_ref()),
            ("mass", exp.dataset.schema.mass.as_ref()),
            ("alpha", exp.dataset.schema.alpha.as_ref()),
            ("dm_dt", exp.dataset.schema.dm_dt.as_ref()),
            ("eta", exp.dataset.schema.eta.as_ref()),
            ("deta_dt", exp.dataset.schema.deta_dt.as_ref()),
            ("dalpha_dt", exp.dataset.schema.dalpha_dt.as_ref()),
            ("dT_dt", exp.dataset.schema.dT_dt.as_ref()),
            ("E", exp.dataset.schema.E.as_ref()),
        ] {
            meta_lines.push(format!(
                "#META\tBIND\t{}\t{}\t{}",
                exp_id_token,
                bind_key,
                encode_opt_string(bind_value.map(String::as_str)),
            ));
        }

        for col in df.columns() {
            let col_name = col.name().to_string();
            let header_name =
                build_unique_series_header(&mut header_counts, &exp_id, &col_name, idx);
            let meta = exp.dataset.schema.columns.get(&col_name);
            let unit = meta.map(|m| m.unit).unwrap_or(Unit::Unknown);
            let origin = meta.map(|m| m.origin).unwrap_or(ColumnOrigin::Imported);

            meta_lines.push(format!(
                "#META\tCOL\t{}\t{}\t{}\t{}\t{}",
                exp_id_token,
                escape_meta_field(&col_name),
                escape_meta_field(&header_name),
                unit_tag(unit),
                origin_tag(origin),
            ));

            let casted = col.cast(&DataType::Float64)?;
            let values = casted.f64()?;
            let mut padded: Vec<Option<f64>> = Vec::with_capacity(max_height);
            for row in 0..max_height {
                if row < values.len() {
                    padded.push(values.get(row));
                } else {
                    padded.push(None);
                }
            }
            export_columns.push(Series::new(header_name.clone().into(), padded).into());
        }
    }

    let mut file = std::fs::File::create(path).map_err(|e| {
        TGADomainError::InvalidOperation(format!(
            "Failed to create series CSV '{}': {}",
            path.display(),
            e
        ))
    })?;
    writeln!(file, "# KiThe TGA Series").map_err(|e| {
        TGADomainError::InvalidOperation(format!(
            "Failed to write series CSV '{}': {}",
            path.display(),
            e
        ))
    })?;
    writeln!(file, "# format_version=1").map_err(|e| {
        TGADomainError::InvalidOperation(format!(
            "Failed to write series CSV '{}': {}",
            path.display(),
            e
        ))
    })?;
    for line in &meta_lines {
        writeln!(file, "{line}").map_err(|e| {
            TGADomainError::InvalidOperation(format!(
                "Failed to write series CSV '{}': {}",
                path.display(),
                e
            ))
        })?;
    }

    if !export_columns.is_empty() {
        let mut df = DataFrame::new(max_height, export_columns)?;
        CsvWriter::new(&mut file)
            .include_header(true)
            .finish(&mut df)?;
    }

    Ok(())
}

pub(super) fn from_csv_series_impl(path: &Path) -> Result<TGASeries, TGADomainError> {
    let text = std::fs::read_to_string(path).map_err(|e| {
        TGADomainError::InvalidOperation(format!(
            "Failed to read series CSV '{}': {}",
            path.display(),
            e
        ))
    })?;

    let mut records: HashMap<String, SeriesExperimentRecord> = HashMap::new();
    let mut order: Vec<String> = Vec::new();
    let mut csv_data_lines: Vec<String> = Vec::new();

    for line in text.lines() {
        if let Some(meta_payload) = line.strip_prefix("#META\t") {
            parse_series_meta_line(meta_payload, &mut records, &mut order)?;
            continue;
        }
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        csv_data_lines.push(line.to_string());
    }

    let data_df = if csv_data_lines.is_empty() {
        None
    } else {
        let tmp = tempfile::NamedTempFile::new().map_err(|e| {
            TGADomainError::InvalidOperation(format!(
                "Failed to create temp file while reading '{}': {}",
                path.display(),
                e
            ))
        })?;
        std::fs::write(tmp.path(), csv_data_lines.join("\n")).map_err(|e| {
            TGADomainError::InvalidOperation(format!(
                "Failed to write temp CSV while reading '{}': {}",
                path.display(),
                e
            ))
        })?;
        let plpath = PlRefPath::try_from_path(tmp.path())?;
        Some(
            LazyCsvReader::new(plpath)
                .with_has_header(true)
                .finish()?
                .collect()?,
        )
    };

    let mut series = TGASeries::new();
    for id in order {
        let rec = records.get(&id).cloned().ok_or_else(|| {
            TGADomainError::InvalidOperation(format!(
                "Series metadata is inconsistent for experiment id '{}'",
                id
            ))
        })?;

        let mut exp_columns: Vec<Column> = Vec::new();
        if let Some(df) = data_df.as_ref() {
            for c in &rec.columns {
                let column = df.column(&c.header_name).map_err(|_| {
                    TGADomainError::InvalidOperation(format!(
                        "Series CSV is missing expected header '{}'",
                        c.header_name
                    ))
                })?;
                let mut owned = column.clone();
                owned.rename(c.column_name.clone().into());
                exp_columns.push(owned);
            }
        }

        let mut exp_df = if exp_columns.is_empty() {
            DataFrame::default()
        } else {
            let height = exp_columns.iter().map(|c| c.len()).max().unwrap_or(0);
            DataFrame::new(height, exp_columns)?
        };
        if rec.len < exp_df.height() {
            exp_df = exp_df.slice(0, rec.len);
        }
        if rec.len > exp_df.height() {
            return Err(TGADomainError::InvalidOperation(format!(
                "Series CSV for experiment '{}' has less data rows than declared length {}",
                id, rec.len
            )));
        }

        let mut schema_columns: HashMap<String, ColumnMeta> = HashMap::new();
        for c in &rec.columns {
            schema_columns.insert(
                c.column_name.clone(),
                ColumnMeta::new(
                    c.column_name.clone(),
                    c.unit,
                    c.origin,
                    ColumnNature::Unknown,
                    ColumnProvenance::manual(
                        c.column_name.clone(),
                        "from_csv_series",
                        Some("restored from series csv".to_string()),
                    ),
                ),
            );
        }

        let mut dataset = TGADataset {
            frame: exp_df.lazy(),
            schema: TGASchema {
                columns: schema_columns,
                time: rec.binds.get("time").cloned().unwrap_or(None),
                temperature: rec.binds.get("temperature").cloned().unwrap_or(None),
                mass: rec.binds.get("mass").cloned().unwrap_or(None),
                alpha: rec.binds.get("alpha").cloned().unwrap_or(None),
                dm_dt: rec.binds.get("dm_dt").cloned().unwrap_or(None),
                eta: rec.binds.get("eta").cloned().unwrap_or(None),
                deta_dt: rec.binds.get("deta_dt").cloned().unwrap_or(None),
                dalpha_dt: rec.binds.get("dalpha_dt").cloned().unwrap_or(None),
                dT_dt: rec.binds.get("dT_dt").cloned().unwrap_or(None),
                E: rec.binds.get("E").cloned().unwrap_or(None),
                R2: rec.binds.get("R2").cloned().unwrap_or(None),
            },
            oneframeplot: None,
            history_of_operations: History::new(),
            undo_stack: Vec::new(),
            undo_snapshot_latch: false,
        };
        dataset.initialize_column_provenance();
        dataset.log_operation(
            "from_csv_series",
            crate::Kinetics::experimental_kinetics::one_experiment_dataset::AffectedColumns::All,
            None,
            format!(
                "Restored experiment '{}' from series CSV '{}'",
                rec.id,
                path.display()
            ),
            false,
        );
        let experiment = TGAExperiment {
            dataset,
            meta: ExperimentMeta {
                id: rec.id,
                heating_rate: rec.heating_rate,
                isothermal_temperature: rec.isothermal_temperature,
                comment: rec.comment,
            },
        };
        series.push(experiment)?;
    }

    Ok(series)
}

fn parse_series_meta_line(
    payload: &str,
    records: &mut HashMap<String, SeriesExperimentRecord>,
    order: &mut Vec<String>,
) -> Result<(), TGADomainError> {
    let parts: Vec<&str> = payload.split('\t').collect();
    if parts.is_empty() {
        return Ok(());
    }

    match parts[0] {
        "EXP" => {
            if parts.len() != 5 {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Invalid EXP metadata line: '{}'",
                    payload
                )));
            }
            let id = unescape_meta_field(parts[1])?;
            if !records.contains_key(&id) {
                order.push(id.clone());
            }
            let rec = records
                .entry(id.clone())
                .or_insert_with(|| SeriesExperimentRecord {
                    id: id.clone(),
                    ..Default::default()
                });
            rec.id = id;
            rec.heating_rate = decode_opt_f64(parts[2])?;
            rec.isothermal_temperature = decode_opt_f64(parts[3])?;
            rec.comment = decode_opt_string(parts[4])?;
        }
        "LEN" => {
            if parts.len() != 3 {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Invalid LEN metadata line: '{}'",
                    payload
                )));
            }
            let id = unescape_meta_field(parts[1])?;
            if !records.contains_key(&id) {
                order.push(id.clone());
            }
            let rec = records
                .entry(id.clone())
                .or_insert_with(|| SeriesExperimentRecord {
                    id,
                    ..Default::default()
                });
            rec.len = parts[2].parse::<usize>().map_err(|e| {
                TGADomainError::InvalidOperation(format!(
                    "Invalid LEN value '{}' in metadata: {}",
                    parts[2], e
                ))
            })?;
        }
        "BIND" => {
            if parts.len() != 4 {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Invalid BIND metadata line: '{}'",
                    payload
                )));
            }
            let id = unescape_meta_field(parts[1])?;
            if !records.contains_key(&id) {
                order.push(id.clone());
            }
            let rec = records
                .entry(id.clone())
                .or_insert_with(|| SeriesExperimentRecord {
                    id,
                    ..Default::default()
                });
            rec.binds
                .insert(parts[2].to_string(), decode_opt_string(parts[3])?);
        }
        "COL" => {
            if parts.len() != 6 {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Invalid COL metadata line: '{}'",
                    payload
                )));
            }
            let id = unescape_meta_field(parts[1])?;
            if !records.contains_key(&id) {
                order.push(id.clone());
            }
            let rec = records
                .entry(id.clone())
                .or_insert_with(|| SeriesExperimentRecord {
                    id,
                    ..Default::default()
                });
            rec.columns.push(SeriesColumnRecord {
                column_name: unescape_meta_field(parts[2])?,
                header_name: unescape_meta_field(parts[3])?,
                unit: parse_unit_tag(parts[4])?,
                origin: parse_origin_tag(parts[5])?,
            });
        }
        _ => {}
    }

    Ok(())
}

fn build_unique_series_header(
    header_counts: &mut HashMap<String, usize>,
    exp_id: &str,
    col_name: &str,
    exp_idx: usize,
) -> String {
    let left = sanitize_for_header(exp_id);
    let right = sanitize_for_header(col_name);
    let mut base = format!("{left}_{right}");
    if left == "unnamed" {
        base = format!("exp{exp_idx}_{right}");
    }

    let count = header_counts.entry(base.clone()).or_insert(0);
    *count += 1;
    if *count == 1 {
        base
    } else {
        format!("{}_{}", base, count)
    }
}

fn sanitize_for_header(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len());
    for ch in raw.chars() {
        if ch.is_ascii_alphanumeric() || ch == '_' {
            out.push(ch);
        } else {
            out.push('_');
        }
    }
    let out = out.trim_matches('_').to_string();
    if out.is_empty() {
        "unnamed".to_string()
    } else {
        out
    }
}

fn unit_tag(unit: Unit) -> &'static str {
    match unit {
        Unit::Second => "Second",
        Unit::Minute => "Minute",
        Unit::Hour => "Hour",
        Unit::Kelvin => "Kelvin",
        Unit::Celsius => "Celsius",
        Unit::MilliVolt => "MilliVolt",
        Unit::Milligram => "Milligram",
        Unit::MilligramPerSecond => "MilligramPerSecond",
        Unit::MilligramPerMinute => "MilligramPerMinute",
        Unit::MilligramPerHour => "MilligramPerHour",
        Unit::KelvinPerSecond => "KelvinPerSecond",
        Unit::KelvinPerMinute => "KelvinPerMinute",
        Unit::KelvinPerHour => "KelvinPerHour",
        Unit::CelsiusPerSecond => "CelsiusPerSecond",
        Unit::CelsiusPerMinute => "CelsiusPerMinute",
        Unit::CelsiusPerHour => "CelsiusPerHour",
        Unit::PerSecond => "PerSecond",
        Unit::PerMinute => "PerMinute",
        Unit::PerHour => "PerHour",
        Unit::Gram => "Gram",
        Unit::Dimensionless => "Dimensionless",
        Unit::Unknown => "Unknown",
    }
}

fn parse_unit_tag(tag: &str) -> Result<Unit, TGADomainError> {
    match tag {
        "Second" => Ok(Unit::Second),
        "Minute" => Ok(Unit::Minute),
        "Hour" => Ok(Unit::Hour),
        "Kelvin" => Ok(Unit::Kelvin),
        "Celsius" => Ok(Unit::Celsius),
        "MilliVolt" => Ok(Unit::MilliVolt),
        "Milligram" => Ok(Unit::Milligram),
        "MilligramPerSecond" => Ok(Unit::MilligramPerSecond),
        "MilligramPerMinute" => Ok(Unit::MilligramPerMinute),
        "MilligramPerHour" => Ok(Unit::MilligramPerHour),
        "KelvinPerSecond" => Ok(Unit::KelvinPerSecond),
        "KelvinPerMinute" => Ok(Unit::KelvinPerMinute),
        "KelvinPerHour" => Ok(Unit::KelvinPerHour),
        "CelsiusPerSecond" => Ok(Unit::CelsiusPerSecond),
        "CelsiusPerMinute" => Ok(Unit::CelsiusPerMinute),
        "CelsiusPerHour" => Ok(Unit::CelsiusPerHour),
        "PerSecond" => Ok(Unit::PerSecond),
        "PerMinute" => Ok(Unit::PerMinute),
        "PerHour" => Ok(Unit::PerHour),
        "Gram" => Ok(Unit::Gram),
        "Dimensionless" => Ok(Unit::Dimensionless),
        "Unknown" => Ok(Unit::Unknown),
        "-" => Ok(Unit::Dimensionless),
        _ => Err(TGADomainError::InvalidOperation(format!(
            "Unknown unit tag '{}'",
            tag
        ))),
    }
}

fn origin_tag(origin: ColumnOrigin) -> &'static str {
    match origin {
        ColumnOrigin::Raw => "Raw",
        ColumnOrigin::PolarsDerived => "PolarsDerived",
        ColumnOrigin::NumericDerived => "NumericDerived",
        ColumnOrigin::Imported => "Imported",
    }
}

fn parse_origin_tag(tag: &str) -> Result<ColumnOrigin, TGADomainError> {
    match tag {
        "Raw" => Ok(ColumnOrigin::Raw),
        "PolarsDerived" => Ok(ColumnOrigin::PolarsDerived),
        "NumericDerived" => Ok(ColumnOrigin::NumericDerived),
        "Imported" => Ok(ColumnOrigin::Imported),
        _ => Err(TGADomainError::InvalidOperation(format!(
            "Unknown column origin tag '{}'",
            tag
        ))),
    }
}

fn encode_opt_f64(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("1:{x}"),
        None => "0".to_string(),
    }
}

fn decode_opt_f64(token: &str) -> Result<Option<f64>, TGADomainError> {
    if token == "0" {
        return Ok(None);
    }
    let raw = token.strip_prefix("1:").ok_or_else(|| {
        TGADomainError::InvalidOperation(format!("Invalid optional-f64 token '{}'", token))
    })?;
    raw.parse::<f64>().map(Some).map_err(|e| {
        TGADomainError::InvalidOperation(format!("Invalid optional-f64 value '{}': {}", raw, e))
    })
}

fn encode_opt_string(v: Option<&str>) -> String {
    match v {
        Some(x) => format!("1:{}", escape_meta_field(x)),
        None => "0".to_string(),
    }
}

fn decode_opt_string(token: &str) -> Result<Option<String>, TGADomainError> {
    if token == "0" {
        return Ok(None);
    }
    let raw = token.strip_prefix("1:").ok_or_else(|| {
        TGADomainError::InvalidOperation(format!("Invalid optional-string token '{}'", token))
    })?;
    Ok(Some(unescape_meta_field(raw)?))
}

fn escape_meta_field(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len());
    for ch in raw.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            _ => out.push(ch),
        }
    }
    out
}

fn unescape_meta_field(raw: &str) -> Result<String, TGADomainError> {
    let mut out = String::with_capacity(raw.len());
    let mut chars = raw.chars();
    while let Some(ch) = chars.next() {
        if ch != '\\' {
            out.push(ch);
            continue;
        }

        let next = chars.next().ok_or_else(|| {
            TGADomainError::InvalidOperation(format!("Invalid escaped metadata token '{}'", raw))
        })?;
        match next {
            '\\' => out.push('\\'),
            't' => out.push('\t'),
            'n' => out.push('\n'),
            'r' => out.push('\r'),
            _ => {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Unsupported escape sequence '\\{}' in metadata token '{}'",
                    next, raw
                )));
            }
        }
    }
    Ok(out)
}
