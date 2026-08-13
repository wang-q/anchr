//! Per-position base grouping, replicating fastqc 0.12.1 `BaseGroup`.
//!
//! The installed fastqc 0.12.1 (`cbp`) differs from the `FastQC-master`
//! sources: groups start with 9 ungrouped positions, then use a fixed
//! interval from `getLinearInterval` (no progressive widening), with a
//! special `10..interval-1` group when the interval exceeds 10. This was
//! recovered from the installed `.class` bytecode and verified against
//! golden output for lengths 75/76/108/150/300/2000.

/// 1-based inclusive base ranges, e.g. `(1, 1)`, `(10, 11)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseGroup {
    pub start: u32,
    pub end: u32,
}

impl std::fmt::Display for BaseGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.start == self.end {
            write!(f, "{}", self.start)
        } else {
            write!(f, "{}-{}", self.start, self.end)
        }
    }
}

/// Replicates fastqc `BaseGroup.makeBaseGroups(max_length)`.
pub fn make_base_groups(max_length: u32) -> Vec<BaseGroup> {
    if max_length <= 75 {
        return (1..=max_length)
            .map(|i| BaseGroup { start: i, end: i })
            .collect();
    }

    let interval = linear_interval(max_length);
    let mut groups = Vec::new();
    let mut start = 1;
    while start <= max_length {
        let mut end = start + interval - 1;
        if start < 10 {
            end = start; // first nine positions ungrouped
        } else if start == 10 && interval > 10 {
            end = interval - 1; // e.g. 10-49 for interval 50
        }
        if end > max_length {
            end = max_length;
        }
        groups.push(BaseGroup { start, end });
        if start < 10 {
            start += 1; // first nine positions ungrouped (step by one)
        } else if start == 10 && interval > 10 {
            start = interval; // skip to the interval boundary
        } else {
            start += interval;
        }
    }
    groups
}

/// Replicates fastqc `BaseGroup.getLinearInterval`: smallest interval from
/// {2, 5, 10} × 10^k whose `9 + ceil((len-9)/interval)` group count < 75.
fn linear_interval(length: u32) -> u32 {
    let bases = [2u32, 5, 10];
    let mut multiplier = 1u32;
    loop {
        for &b in &bases {
            let interval = b * multiplier;
            let group_count = 9 + (length - 9).div_ceil(interval);
            if group_count < 75 {
                return interval;
            }
        }
        multiplier *= 10;
        assert!(
            multiplier < 10_000_000,
            "no sensible interval for length {length}"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ungrouped_up_to_75() {
        assert_eq!(make_base_groups(75).len(), 75);
        assert_eq!(make_base_groups(75)[74], BaseGroup { start: 75, end: 75 });
    }

    #[test]
    fn length_108_matches_golden() {
        let groups = make_base_groups(108);
        // 1..9 ungrouped, then pairs, final single
        assert_eq!(groups.len(), 59);
        assert_eq!(groups[0], BaseGroup { start: 1, end: 1 });
        assert_eq!(groups[8], BaseGroup { start: 9, end: 9 });
        assert_eq!(groups[9], BaseGroup { start: 10, end: 11 });
        assert_eq!(
            groups[57],
            BaseGroup {
                start: 106,
                end: 107
            }
        );
        assert_eq!(
            groups[58],
            BaseGroup {
                start: 108,
                end: 108
            }
        );
    }

    #[test]
    fn length_150_interval_five() {
        let groups = make_base_groups(150);
        assert_eq!(groups.len(), 38);
        assert_eq!(groups[9], BaseGroup { start: 10, end: 14 });
        assert_eq!(
            groups[37],
            BaseGroup {
                start: 150,
                end: 150
            }
        );
    }

    #[test]
    fn length_2000_special_ten_group() {
        let groups = make_base_groups(2000);
        assert_eq!(groups.len(), 50);
        assert_eq!(groups[9], BaseGroup { start: 10, end: 49 });
        assert_eq!(groups[10], BaseGroup { start: 50, end: 99 });
        assert_eq!(
            groups[49],
            BaseGroup {
                start: 2000,
                end: 2000
            }
        );
    }
}
