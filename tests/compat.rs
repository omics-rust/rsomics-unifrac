use std::path::PathBuf;
use std::process::Command;

fn golden(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn run_binary(table: &PathBuf, tree: &PathBuf, flags: &[&str]) -> String {
    let bin = env!("CARGO_BIN_EXE_rsomics-unifrac");
    let out = Command::new(bin)
        .arg(table)
        .arg("--tree")
        .arg(tree)
        .args(flags)
        .output()
        .expect("run rsomics-unifrac");
    assert!(out.status.success(), "binary failed: {out:?}");
    String::from_utf8(out.stdout).unwrap()
}

/// Square matrix → flat (id, id, value) triples, header-aware.
fn parse(text: &str) -> Vec<(String, String, f64)> {
    let mut lines = text.lines().filter(|l| !l.is_empty());
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .skip(1)
        .map(str::to_string)
        .collect();
    let mut out = Vec::new();
    for line in lines {
        let mut f = line.split('\t');
        let row = f.next().unwrap().to_string();
        for (col, cell) in f.enumerate() {
            out.push((row.clone(), header[col].clone(), cell.parse().unwrap()));
        }
    }
    out
}

fn assert_close(got: &str, want: &str, tol: f64) {
    let g = parse(got);
    let w = parse(want);
    assert_eq!(g.len(), w.len());
    for ((gr, gc, gv), (wr, wc, wv)) in g.iter().zip(&w) {
        assert_eq!((gr, gc), (wr, wc));
        assert!((gv - wv).abs() < tol, "{gr}/{gc}: ours {gv} vs skbio {wv}");
    }
}

const MODES: [(&str, &[&str]); 3] = [
    ("unweighted.tsv", &[]),
    ("weighted.tsv", &["--weighted"]),
    ("weighted_normalized.tsv", &["--weighted", "--normalized"]),
];

#[test]
fn matches_committed_skbio_golden() {
    let table = golden("table.tsv");
    let tree = golden("tree.nwk");
    for (expected, flags) in MODES {
        let got = run_binary(&table, &tree, flags);
        let want = std::fs::read_to_string(golden(expected)).unwrap();
        assert_close(&got, &want, 1e-6);
    }
}

#[test]
fn matches_live_skbio() {
    let python = std::env::var("SKBIO_PYTHON").unwrap_or_else(|_| "python3".to_string());
    let probe = Command::new(&python).args(["-c", "import skbio"]).status();
    match probe {
        Ok(s) if s.success() => {}
        _ => {
            eprintln!("SKIP matches_live_skbio: scikit-bio not importable via '{python}'");
            return;
        }
    }

    let table = golden("table.tsv");
    let tree = golden("tree.nwk");
    let scratch = std::env::temp_dir().join("rsomics_unifrac_live");
    let script = r#"
import sys
import numpy as np
from skbio import TreeNode
from skbio.diversity import beta_diversity
table, treepath, metric, norm = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
tree = TreeNode.read([open(treepath).read().strip()])
lines = [l for l in open(table).read().splitlines() if l and not l.startswith('#')]
samples = lines[0].split('\t')[1:]
taxa, rows = [], []
for l in lines[1:]:
    f = l.split('\t')
    taxa.append(f[0]); rows.append([int(x) for x in f[1:]])
counts = np.array(rows).T
kw = {'normalized': True} if norm == '1' else {}
dm = beta_diversity(metric, counts, ids=samples, taxa=taxa, tree=tree, **kw)
out = ['\t' + '\t'.join(dm.ids)]
for i, rid in enumerate(dm.ids):
    out.append(rid + '\t' + '\t'.join(repr(float(v)) for v in dm.data[i]))
print('\n'.join(out))
"#;
    let cases: [(&str, &str, &[&str]); 3] = [
        ("unweighted_unifrac", "0", &[]),
        ("weighted_unifrac", "0", &["--weighted"]),
        ("weighted_unifrac", "1", &["--weighted", "--normalized"]),
    ];
    for (metric, norm, flags) in cases {
        let oracle = Command::new(&python)
            .args(["-c", script])
            .arg(&table)
            .arg(&tree)
            .arg(metric)
            .arg(norm)
            .env("MPLCONFIGDIR", &scratch)
            .output()
            .expect("skbio output");
        assert!(oracle.status.success(), "skbio failed: {oracle:?}");
        let want = String::from_utf8(oracle.stdout).unwrap();
        let got = run_binary(&table, &tree, flags);
        assert_close(&got, &want, 1e-9);
    }
}
