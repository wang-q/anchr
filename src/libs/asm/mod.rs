//! De novo assembly and read-refinement algorithms.

pub mod assemble;
pub mod dfa;
pub mod extend;
pub mod multik;
pub mod refine;
pub mod table;

/// `Tadpole.isJunction(max, second)`: depth-ratio branch detection.
///
/// Shared by the refine (ecc) and assemble (contig/unitig) modes; the
/// formula lives here so the two option types cannot drift apart.
pub(crate) fn is_junction(
    max: u32,
    second: u32,
    branch_mult1: f32,
    branch_mult2: f32,
    branch_lower_const: usize,
    min_count_extend: usize,
) -> bool {
    if second < 1
        || (second as f32) * branch_mult1 < max as f32
        || (second <= branch_lower_const as u32
            && (max as f32) >= (min_count_extend as f32).max(second as f32 * branch_mult2))
    {
        return false;
    }
    true
}
