use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use rsomics_phylo_tree::Tree;

use rsomics_unifrac::{Config, Mode, run};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-unifrac", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Count table (feature-by-sample TSV); reads stdin when "-" or omitted.
    #[arg(default_value = "-")]
    input: PathBuf,

    /// Rooted Newick tree whose tips are the OTU/taxon IDs.
    #[arg(long)]
    tree: PathBuf,

    /// Weighted (quantitative) UniFrac instead of unweighted (qualitative).
    #[arg(long, default_value_t = false)]
    weighted: bool,

    /// Branch-length-normalize the weighted distance into [0, 1] (needs --weighted).
    #[arg(long, default_value_t = false)]
    normalized: bool,

    /// Treat the input table as comma-separated instead of tab-separated.
    #[arg(long, default_value_t = false)]
    csv: bool,

    /// Output path; writes stdout when "-".
    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        if self.normalized && !self.weighted {
            return Err(RsomicsError::InvalidInput(
                "--normalized applies only to weighted UniFrac (pass --weighted)".into(),
            ));
        }
        let mode = match (self.weighted, self.normalized) {
            (false, _) => Mode::Unweighted,
            (true, false) => Mode::Weighted,
            (true, true) => Mode::WeightedNormalized,
        };

        let newick = fs::read_to_string(&self.tree)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", self.tree.display())))?;
        let tree = Tree::from_newick(&newick)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", self.tree.display())))?;

        let cfg = Config {
            mode,
            delim: if self.csv { ',' } else { '\t' },
        };

        let reader: Box<dyn std::io::BufRead> = if self.input.as_os_str() == "-" {
            Box::new(BufReader::new(std::io::stdin().lock()))
        } else {
            Box::new(BufReader::new(File::open(&self.input).map_err(|e| {
                RsomicsError::InvalidInput(format!("{}: {e}", self.input.display()))
            })?))
        };
        let mut out: Box<dyn Write> = if self.output == "-" {
            Box::new(BufWriter::new(std::io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                File::create(&self.output).map_err(RsomicsError::Io)?,
            ))
        };
        run(reader, &mut out, &tree, &cfg)?;
        out.flush().map_err(RsomicsError::Io)
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Pairwise UniFrac phylogenetic beta-diversity from a count table and a rooted tree.",
    origin: Some(Origin {
        upstream: "scikit-bio skbio.diversity.beta_diversity (un/weighted_unifrac)",
        upstream_license: "BSD-3-Clause",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1128/AEM.71.12.8228-8235.2005"),
    }),
    usage_lines: &["[table.tsv] --tree tree.nwk [--weighted [--normalized]] [-o dm.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "tree",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("path"),
                required: true,
                default: None,
                description: "Rooted Newick tree whose tips are the OTU/taxon IDs.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "weighted",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: Some("false"),
                description: "Weighted (quantitative) UniFrac instead of unweighted.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "normalized",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: Some("false"),
                description: "Branch-length-normalize the weighted distance into [0, 1].",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "csv",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: Some("false"),
                description: "Parse the table as comma-separated.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("path"),
                required: false,
                default: Some("-"),
                description: "Output path (- for stdout).",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Unweighted UniFrac distance matrix",
            command: "rsomics-unifrac counts.tsv --tree tree.nwk",
        },
        Example {
            description: "Weighted-normalized, to a file",
            command: "rsomics-unifrac counts.tsv --tree tree.nwk --weighted --normalized -o dm.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
