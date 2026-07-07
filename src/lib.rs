use std::collections::HashMap;
use std::io::{BufRead, Write};

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};
use rsomics_phylo_tree::Tree;

mod table;
pub use table::CountTable;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Unweighted,
    Weighted,
    WeightedNormalized,
}

pub struct Config {
    pub mode: Mode,
    pub delim: char,
}

/// The tree flattened for UniFrac: a postorder node list with branch lengths,
/// each node's child span, and the root-distance of every tip.
///
/// Counts are pushed to tips then summed up the tree (postorder reduce, per
/// Hamady, Lozupone & Knight 2010, Fast UniFrac), so every node carries the
/// total count of the leaves below it.
struct UniFracTree {
    branch_length: Vec<f64>,
    /// Parent slot of each postorder node; the root maps to itself.
    parent: Vec<usize>,
    /// Root-distance of each node, masked to 0 at internal nodes.
    tip_root_distance: Vec<f64>,
    /// Postorder slot of each tip, keyed by tip name.
    tip_slot: HashMap<String, usize>,
    n_nodes: usize,
}

impl UniFracTree {
    fn build(tree: &Tree) -> Result<UniFracTree> {
        if tree.nodes[tree.root].children.len() > 2 {
            return Err(RsomicsError::InvalidInput(
                "tree is not rooted (root has more than two children)".into(),
            ));
        }

        let n_nodes = tree.nodes.len();
        let mut postorder = Vec::with_capacity(n_nodes);
        let mut slot_of = vec![0usize; n_nodes];
        let mut stack = vec![(tree.root, false)];
        while let Some((id, expanded)) = stack.pop() {
            if expanded {
                slot_of[id] = postorder.len();
                postorder.push(id);
            } else {
                stack.push((id, true));
                for &c in &tree.nodes[id].children {
                    stack.push((c, false));
                }
            }
        }

        let mut branch_length = vec![0.0f64; n_nodes];
        let mut parent = vec![0usize; n_nodes];
        let mut tip_root_distance = vec![0.0f64; n_nodes];
        let mut tip_slot = HashMap::new();

        for (slot, &id) in postorder.iter().enumerate() {
            let node = &tree.nodes[id];
            parent[slot] = node.parent.map_or(slot, |p| slot_of[p]);
            match node.branch_length {
                Some(bl) => branch_length[slot] = bl,
                None if id == tree.root => {}
                None => {
                    return Err(RsomicsError::InvalidInput(
                        "every non-root node must have a branch length".into(),
                    ));
                }
            }
            if node.children.is_empty() {
                let name = node
                    .name
                    .as_deref()
                    .ok_or_else(|| RsomicsError::InvalidInput("a tip has no name".into()))?;
                if tip_slot.insert(name.to_string(), slot).is_some() {
                    return Err(RsomicsError::InvalidInput(format!(
                        "duplicate tip name '{name}' in the tree"
                    )));
                }
            }
        }

        // The root's own pendant length is part of every tip-to-root distance
        // (skbio folds root_bl into `_tip_distances`), so seed it before the
        // accumulation carries it down; the mask below zeroes the root itself.
        let root_slot = n_nodes - 1;
        tip_root_distance[root_slot] = branch_length[root_slot];

        // Preorder accumulation of root-distance (parent precedes child), masked
        // to leave only tips nonzero, per skbio's `_tip_distances`.
        for slot in (0..n_nodes).rev() {
            if parent[slot] != slot {
                tip_root_distance[slot] = branch_length[slot] + tip_root_distance[parent[slot]];
            }
        }
        for (slot, &id) in postorder.iter().enumerate() {
            if !tree.nodes[id].children.is_empty() {
                tip_root_distance[slot] = 0.0;
            }
        }

        Ok(UniFracTree {
            branch_length,
            parent,
            tip_root_distance,
            tip_slot,
            n_nodes,
        })
    }

    /// Push each sample's tip counts into postorder slots, then sum up the tree.
    /// Returns one descendant-summed node-count vector per sample plus the
    /// per-sample tip total.
    fn node_counts(&self, table: &CountTable, tip_of_row: &[usize]) -> (Vec<Vec<f64>>, Vec<f64>) {
        let n_samples = table.sample_names.len();
        let mut per_sample = vec![vec![0.0f64; self.n_nodes]; n_samples];
        let mut totals = vec![0.0f64; n_samples];

        for (s, counts) in per_sample.iter_mut().enumerate() {
            let col = &table.columns[s];
            for (row, &c) in col.iter().enumerate() {
                if c != 0.0 {
                    counts[tip_of_row[row]] += c;
                    totals[s] += c;
                }
            }
            // Postorder push: each node adds its total to its parent (children
            // always precede their parent, so one forward pass suffices).
            for slot in 0..self.n_nodes {
                let p = self.parent[slot];
                if p != slot {
                    counts[p] += counts[slot];
                }
            }
        }
        (per_sample, totals)
    }

    fn unweighted(&self, u: &[f64], v: &[f64]) -> f64 {
        let mut unique = 0.0;
        let mut observed = 0.0;
        for slot in 0..self.n_nodes {
            let up = u[slot] > 0.0;
            let vp = v[slot] > 0.0;
            if up || vp {
                let bl = self.branch_length[slot];
                observed += bl;
                if up != vp {
                    unique += bl;
                }
            }
        }
        if observed == 0.0 {
            0.0
        } else {
            unique / observed
        }
    }

    fn weighted(&self, u: &[f64], ut: f64, v: &[f64], vt: f64, normalized: bool) -> f64 {
        if normalized && ut == 0.0 && vt == 0.0 {
            return 0.0;
        }
        let mut wu = 0.0;
        let mut norm = 0.0;
        for slot in 0..self.n_nodes {
            let up = if ut > 0.0 { u[slot] / ut } else { 0.0 };
            let vp = if vt > 0.0 { v[slot] / vt } else { 0.0 };
            wu += self.branch_length[slot] * (up - vp).abs();
            if normalized {
                norm += self.tip_root_distance[slot] * (up + vp);
            }
        }
        if normalized { wu / norm } else { wu }
    }
}

/// A symmetric pairwise distance matrix over the samples, row-major dense.
pub struct DistanceMatrix {
    ids: Vec<String>,
    data: Vec<f64>,
}

impl DistanceMatrix {
    fn compute(table: &CountTable, tree: &UniFracTree, mode: Mode) -> DistanceMatrix {
        let n = table.sample_names.len();
        let tip_of_row: Vec<usize> = table.feature_ids.iter().map(|t| tree.tip_slot[t]).collect();
        let (counts, totals) = tree.node_counts(table, &tip_of_row);

        let pairs: Vec<(usize, usize)> = (0..n)
            .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
            .collect();
        let upper: Vec<f64> = pairs
            .par_iter()
            .map(|&(i, j)| match mode {
                Mode::Unweighted => tree.unweighted(&counts[i], &counts[j]),
                Mode::Weighted => {
                    tree.weighted(&counts[i], totals[i], &counts[j], totals[j], false)
                }
                Mode::WeightedNormalized => {
                    tree.weighted(&counts[i], totals[i], &counts[j], totals[j], true)
                }
            })
            .collect();

        let mut data = vec![0.0f64; n * n];
        for (&(i, j), &d) in pairs.iter().zip(&upper) {
            data[i * n + j] = d;
            data[j * n + i] = d;
        }
        DistanceMatrix {
            ids: table.sample_names.clone(),
            data,
        }
    }

    /// scikit-bio `DistanceMatrix` TSV (LSMat): empty top-left cell, sample IDs
    /// as the header, then one labelled row per sample with `repr(float)` cells.
    ///
    /// # Errors
    /// Propagates write errors.
    pub fn write_tsv<W: Write>(&self, mut out: W) -> Result<()> {
        let n = self.ids.len();
        for id in &self.ids {
            write!(out, "\t{id}").map_err(RsomicsError::Io)?;
        }
        writeln!(out).map_err(RsomicsError::Io)?;
        let mut row = String::new();
        for i in 0..n {
            row.clear();
            row.push_str(&self.ids[i]);
            for j in 0..n {
                row.push('\t');
                push_pyrepr(&mut row, self.data[i * n + j]);
            }
            writeln!(out, "{row}").map_err(RsomicsError::Io)?;
        }
        Ok(())
    }
}

/// Append `x` as Python's `repr(float)`: shortest round-trip decimal, integer
/// floats keep a trailing `.0`, NaN renders lowercase `nan`.
fn push_pyrepr(buf: &mut String, x: f64) {
    use std::fmt::Write;
    if x.is_nan() {
        buf.push_str("nan");
        return;
    }
    if x.is_infinite() {
        buf.push_str(if x < 0.0 { "-inf" } else { "inf" });
        return;
    }
    let start = buf.len();
    let _ = write!(buf, "{x}");
    if !buf[start..].contains(['.', 'e', 'E']) {
        buf.push_str(".0");
    }
}

/// # Errors
/// Propagates parse, tree-build, and write errors; errors if a table taxon is
/// not a tip in the tree.
pub fn run<R: BufRead, W: Write>(reader: R, out: W, tree: &Tree, cfg: &Config) -> Result<()> {
    let table = CountTable::parse(reader, cfg.delim)?;
    let uf = UniFracTree::build(tree)?;
    for taxon in &table.feature_ids {
        if !uf.tip_slot.contains_key(taxon) {
            return Err(RsomicsError::InvalidInput(format!(
                "taxon '{taxon}' from the count table is not a tip in the tree"
            )));
        }
    }
    let dm = DistanceMatrix::compute(&table, &uf, cfg.mode);
    dm.write_tsv(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const NWK: &str = "(((((U1:0.5,U2:0.5):0.5,U3:1.0):1.0):0.0,(U4:0.75,(U5:0.5,((U6:0.33,U7:0.62):0.5,U8:0.5):0.5):0.5):1.25):0.0)root;";

    fn doc_table() -> &'static str {
        "feature\tu\tv\nU1\t1\t0\nU2\t0\t1\nU3\t0\t1\nU4\t4\t6\nU5\t1\t0\nU6\t2\t1\nU7\t3\t0\nU8\t0\t0\n"
    }

    fn uv(mode: Mode) -> f64 {
        let tree = Tree::from_newick(NWK).unwrap();
        let cfg = Config { mode, delim: '\t' };
        let mut out = Vec::new();
        run(std::io::Cursor::new(doc_table()), &mut out, &tree, &cfg).unwrap();
        let s = String::from_utf8(out).unwrap();
        let line = s.lines().nth(1).unwrap();
        line.split('\t').nth(2).unwrap().parse().unwrap()
    }

    #[test]
    fn unweighted_matches_skbio() {
        assert!((uv(Mode::Unweighted) - 0.369_230_769_230_769_25).abs() < 1e-12);
    }

    #[test]
    fn weighted_matches_skbio() {
        assert!((uv(Mode::Weighted) - 1.543_434_343_434_343_2).abs() < 1e-12);
    }

    #[test]
    fn weighted_normalized_matches_skbio() {
        assert!((uv(Mode::WeightedNormalized) - 0.327_503_429_355_281_15).abs() < 1e-12);
    }

    #[test]
    fn diagonal_is_zero_and_symmetric() {
        let tree = Tree::from_newick(NWK).unwrap();
        let cfg = Config {
            mode: Mode::Weighted,
            delim: '\t',
        };
        let table = CountTable::parse(doc_table().as_bytes(), '\t').unwrap();
        let uf = UniFracTree::build(&tree).unwrap();
        let dm = DistanceMatrix::compute(&table, &uf, cfg.mode);
        let n = dm.ids.len();
        for i in 0..n {
            assert_eq!(dm.data[i * n + i], 0.0);
            for j in 0..n {
                assert_eq!(dm.data[i * n + j], dm.data[j * n + i]);
            }
        }
    }

    fn root_bl_tree_uv(mode: Mode) -> f64 {
        let tree = Tree::from_newick("((A:1.0,B:2.0):0.5,C:3.0)root:9.0;").unwrap();
        let cfg = Config { mode, delim: '\t' };
        let mut out = Vec::new();
        run(
            std::io::Cursor::new("feature\tx\ty\nA\t1\t0\nB\t0\t2\nC\t2\t1\n"),
            &mut out,
            &tree,
            &cfg,
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        s.lines()
            .nth(1)
            .unwrap()
            .split('\t')
            .nth(2)
            .unwrap()
            .parse()
            .unwrap()
    }

    // The root's stored branch length must enter the weighted-normalized
    // denominator (skbio folds root_bl into every tip-to-root distance).
    #[test]
    fn weighted_normalized_folds_root_branch_length() {
        assert!(
            (root_bl_tree_uv(Mode::WeightedNormalized) - 0.122_302_158_273_381_3).abs() < 1e-12
        );
    }

    #[test]
    fn root_branch_length_leaves_other_modes_untouched() {
        assert!((root_bl_tree_uv(Mode::Weighted) - 2.833_333_333_333_333).abs() < 1e-12);
        assert!((root_bl_tree_uv(Mode::Unweighted) - 0.193_548_387_096_774_2).abs() < 1e-12);
    }

    #[test]
    fn duplicate_feature_rejected() {
        let tree = Tree::from_newick("((A:1.0,B:2.0):0.5,C:3.0)root;").unwrap();
        let cfg = Config {
            mode: Mode::WeightedNormalized,
            delim: '\t',
        };
        let mut out = Vec::new();
        let table = "feature\tx\ty\nA\t1\t0\nB\t0\t2\nA\t2\t1\n";
        assert!(run(std::io::Cursor::new(table), &mut out, &tree, &cfg).is_err());
    }

    #[test]
    fn duplicate_sample_rejected() {
        let tree = Tree::from_newick("((A:1.0,B:2.0):0.5,C:3.0)root;").unwrap();
        let cfg = Config {
            mode: Mode::Weighted,
            delim: '\t',
        };
        let mut out = Vec::new();
        let table = "feature\tx\tx\nA\t1\t0\nB\t0\t2\nC\t2\t1\n";
        assert!(run(std::io::Cursor::new(table), &mut out, &tree, &cfg).is_err());
    }

    #[test]
    fn empty_samples_give_zero() {
        let tree = Tree::from_newick("((A:1.0,B:2.0):0.5,C:3.0)root;").unwrap();
        let table = CountTable::parse(
            "feature\tx\ty\nA\t0\t0\nB\t0\t0\nC\t0\t0\n".as_bytes(),
            '\t',
        )
        .unwrap();
        let uf = UniFracTree::build(&tree).unwrap();
        for mode in [Mode::Unweighted, Mode::Weighted, Mode::WeightedNormalized] {
            let dm = DistanceMatrix::compute(&table, &uf, mode);
            assert_eq!(dm.data[1], 0.0);
        }
    }

    #[test]
    fn unknown_taxon_rejected() {
        let tree = Tree::from_newick("((A:1.0,B:2.0):0.5,C:3.0)root;").unwrap();
        let cfg = Config {
            mode: Mode::Unweighted,
            delim: '\t',
        };
        let mut out = Vec::new();
        let table = "feature\tx\nA\t1\nZ\t1\n";
        assert!(run(std::io::Cursor::new(table), &mut out, &tree, &cfg).is_err());
    }

    #[test]
    fn trifurcating_root_rejected() {
        let tree = Tree::from_newick("(A:1.0,B:2.0,C:3.0)root;").unwrap();
        let cfg = Config {
            mode: Mode::Unweighted,
            delim: '\t',
        };
        let mut out = Vec::new();
        let table = "feature\tx\ty\nA\t1\t0\nB\t1\t1\nC\t0\t1\n";
        assert!(run(std::io::Cursor::new(table), &mut out, &tree, &cfg).is_err());
    }
}
