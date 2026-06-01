use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};

/// A feature-by-sample count table.
///
/// Layout matches the rsomics diversity family (scikit-bio / QIIME / phyloseq):
/// first column is the feature (OTU/taxon) ID, the header row names the samples,
/// cell `[feature][sample]` is the count.
pub struct CountTable {
    pub feature_ids: Vec<String>,
    pub sample_names: Vec<String>,
    /// One count vector per sample (column-major), each of length `feature_ids.len()`.
    pub columns: Vec<Vec<f64>>,
}

impl CountTable {
    /// # Errors
    /// Errors on a missing header, a ragged row, a non-numeric cell, or a
    /// negative count.
    pub fn parse<R: BufRead>(reader: R, delim: char) -> Result<CountTable> {
        let mut lines = reader.lines();
        let header = loop {
            match lines.next() {
                Some(line) => {
                    let line = line.map_err(RsomicsError::Io)?;
                    if line.trim().is_empty() || line.starts_with('#') {
                        continue;
                    }
                    break line;
                }
                None => return Err(RsomicsError::InvalidInput("empty count table".into())),
            }
        };
        let sample_names: Vec<String> = header
            .split(delim)
            .skip(1)
            .map(|s| s.trim().to_string())
            .collect();
        if sample_names.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "header has no sample columns (need feature-ID column + ≥1 sample)".into(),
            ));
        }
        let n = sample_names.len();
        let mut feature_ids = Vec::new();
        let mut columns: Vec<Vec<f64>> = vec![Vec::new(); n];
        for (row_idx, line) in lines.enumerate() {
            let line = line.map_err(RsomicsError::Io)?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let mut fields = line.split(delim);
            let feature = fields.next().unwrap_or("").trim().to_string();
            let mut seen = 0usize;
            for (col, field) in fields.enumerate() {
                if col >= n {
                    return Err(RsomicsError::InvalidInput(format!(
                        "row {} (feature '{feature}') has more columns than the header",
                        row_idx + 2
                    )));
                }
                let count: f64 = field.trim().parse().map_err(|_| {
                    RsomicsError::InvalidInput(format!(
                        "row {} (feature '{feature}'), sample '{}': '{}' is not a numeric count",
                        row_idx + 2,
                        sample_names[col],
                        field.trim()
                    ))
                })?;
                if count < 0.0 {
                    return Err(RsomicsError::InvalidInput(format!(
                        "row {} (feature '{feature}'), sample '{}': counts cannot be negative",
                        row_idx + 2,
                        sample_names[col]
                    )));
                }
                columns[col].push(count);
                seen += 1;
            }
            if seen != n {
                return Err(RsomicsError::InvalidInput(format!(
                    "row {} (feature '{feature}') has {seen} count columns, header has {n}",
                    row_idx + 2
                )));
            }
            feature_ids.push(feature);
        }
        Ok(CountTable {
            feature_ids,
            sample_names,
            columns,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_columns() {
        let t = CountTable::parse(
            "feature\tS1\tS2\nOTU1\t4\t0\nOTU2\t0\t10\n".as_bytes(),
            '\t',
        )
        .unwrap();
        assert_eq!(t.sample_names, ["S1", "S2"]);
        assert_eq!(t.feature_ids, ["OTU1", "OTU2"]);
        assert_eq!(t.columns[0], [4.0, 0.0]);
        assert_eq!(t.columns[1], [0.0, 10.0]);
    }

    #[test]
    fn negative_count_errors() {
        assert!(CountTable::parse("feature\tA\nOTU1\t-2\n".as_bytes(), '\t').is_err());
    }

    #[test]
    fn ragged_row_errors() {
        assert!(CountTable::parse("feature\tA\tB\nOTU1\t4\n".as_bytes(), '\t').is_err());
    }
}
