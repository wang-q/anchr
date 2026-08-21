//! BubblePopper over the contig graph (assemble.BubblePopper): divergent
//! paths that reconverge are merged so variant sites do not fragment the
//! main path.

use super::contig::{
    contig_cmp, remove_all_edges, set_associate, set_used, Contig, EdgeRef, DEAD_END, LOOP,
};
use super::AssembleOptions;
use std::collections::HashMap;

/// BubblePopper over the contig graph (assemble.BubblePopper).
struct BubblePopper {
    contigs: Vec<Contig>,
    dest_map: HashMap<usize, Vec<EdgeRef>>,
    k: usize,
    min_len: usize,
    center: usize,
    dest: usize,
    last_mutual_dest: i64,
    last_mutual_dest_orientation: i64,
    expansions: usize,
    contigs_absorbed: usize,
}

impl BubblePopper {
    fn dest_to_edge_map(&self) -> HashMap<usize, Vec<EdgeRef>> {
        let mut map: HashMap<usize, Vec<EdgeRef>> = HashMap::new();
        for c in &self.contigs {
            if c.used || c.associate {
                continue;
            }
            for e in &c.left_edges {
                map.entry(e.borrow().destination)
                    .or_default()
                    .push(e.clone());
            }
            for e in &c.right_edges {
                map.entry(e.borrow().destination)
                    .or_default()
                    .push(e.clone());
            }
        }
        map
    }

    fn expand(&mut self, center_id: usize) -> usize {
        self.center = center_id;
        let mut count = 0;
        while self.expand_right_simple() {
            count += 1;
        }
        while self.contigs[center_id].right_forward_branch() && self.expand_right() {
            count += 1;
            while self.expand_right_simple() {
                count += 1;
            }
        }
        let left_ok = {
            let c = &self.contigs[center_id];
            (c.left_code != LOOP && c.left_code != DEAD_END && !c.left_edges.is_empty())
                || c.left_forward_branch()
        };
        if left_ok {
            let inbound = self.dest_map.get(&center_id).cloned();
            self.contigs[center_id].flip(inbound.as_deref());
            while self.expand_right_simple() {
                count += 1;
            }
            while self.contigs[center_id].right_forward_branch() && self.expand_right() {
                count += 1;
                while self.expand_right_simple() {
                    count += 1;
                }
            }
        }
        count
    }

    fn expand_right_simple(&mut self) -> bool {
        let center_id = self.center;
        let outbound = self.contigs[center_id].right_edges.clone();
        if outbound.is_empty() || self.contigs[center_id].right_code == LOOP || outbound.len() > 1 {
            return false;
        }
        let left_edge = outbound[0].clone();
        let dest_id = left_edge.borrow().destination;
        let dest_right = left_edge.borrow().dest_right();
        if self.contigs[dest_id].used || dest_id == center_id {
            return false;
        }
        let (outbound_right, right_code) = {
            let d = &self.contigs[dest_id];
            if dest_right {
                (d.right_edges.clone(), d.right_code)
            } else {
                (d.left_edges.clone(), d.left_code)
            }
        };
        if right_code == LOOP {
            return false;
        }
        if !outbound_right.is_empty() {
            if outbound_right.len() > 1 {
                return false;
            }
            if outbound_right[0].borrow().destination != center_id {
                return false;
            }
        }
        if self.count_inbound(center_id, true) > 1 {
            return false;
        }
        if self.count_inbound(dest_id, dest_right) > 1 {
            return false;
        }
        if dest_right {
            let inbound = self.dest_map.get(&dest_id).cloned();
            self.contigs[dest_id].flip(inbound.as_deref());
        }
        self.merge(center_id, dest_id, left_edge)
    }

    fn count_inbound(&self, id: usize, dest_right: bool) -> usize {
        self.dest_map
            .get(&id)
            .map(|v| {
                v.iter()
                    .filter(|e| e.borrow().dest_right() == dest_right)
                    .count()
            })
            .unwrap_or(0)
    }

    fn merge(&mut self, left_id: usize, right_id: usize, left_edge: EdgeRef) -> bool {
        let k = self.k;
        let original_left_len = self.contigs[left_id].bases.len();
        let mut bb: Vec<u8> = self.contigs[left_id].bases.clone();
        {
            let eb = left_edge.borrow();
            if eb.bases.len() > 1 {
                bb.extend_from_slice(&eb.bases[..eb.bases.len() - 1]);
            }
        }
        bb.extend_from_slice(&self.contigs[right_id].bases);
        self.contigs[left_id].bases = bb;
        self.contigs[left_id].right_edges.clear();
        let right_right = self.contigs[right_id].right_edges.clone();
        if right_right.is_empty() {
            self.contigs[left_id].right_edges = Vec::new();
        } else {
            for e in &right_right {
                e.borrow_mut().origin = left_id;
            }
            self.contigs[left_id].right_edges = right_right;
        }
        self.redirect_edges(right_id, left_id, true);
        let inbound_right = self.dest_map.get(&right_id).cloned();
        set_used(right_id, inbound_right.as_deref(), &mut self.contigs);
        let right_len = self.contigs[right_id].bases.len();
        let (right_max_cov, right_min_cov, right_code, right_ratio, right_coverage) = {
            let r = &self.contigs[right_id];
            (
                r.max_cov,
                r.min_cov,
                r.right_code,
                r.right_ratio,
                r.coverage,
            )
        };
        {
            let left = &mut self.contigs[left_id];
            left.max_cov = left.max_cov.max(right_max_cov);
            left.min_cov = left.min_cov.min(right_min_cov);
            left.right_code = right_code;
            left.right_ratio = right_ratio;
            let coverage_sum = left.coverage as f64 * (original_left_len - k + 1) as f64
                + right_coverage as f64 * (right_len - k + 1) as f64;
            left.coverage = (coverage_sum / (left.bases.len() - k + 1) as f64) as f32;
        }
        if self.is_loop(left_id) {
            self.contigs[left_id].left_code = LOOP;
            self.contigs[left_id].right_code = LOOP;
            let inbound = self.dest_map.get(&left_id).cloned();
            remove_all_edges(left_id, inbound.as_deref(), &mut self.contigs);
        }
        self.expansions += 1;
        self.contigs_absorbed += 1;
        true
    }

    fn redirect_edges(&mut self, from: usize, to: usize, dest_right: bool) {
        if from == to {
            return;
        }
        let Some(inbound_from) = self.dest_map.remove(&from) else {
            return;
        };
        let mut inbound_to = self.dest_map.get(&to).cloned().unwrap_or_default();
        for e in &inbound_from {
            if e.borrow().dest_right() == dest_right {
                e.borrow_mut().destination = to;
                inbound_to.push(e.clone());
            }
        }
        if inbound_to.is_empty() {
            self.dest_map.remove(&to);
        } else {
            self.dest_map.insert(to, inbound_to);
        }
    }

    fn is_loop(&self, id: usize) -> bool {
        let c = &self.contigs[id];
        if c.left_code == LOOP && c.right_code == LOOP {
            return true;
        }
        if c.left_edges.len() != 1 || c.right_edges.len() != 1 {
            return false;
        }
        for e in &c.left_edges {
            let e = e.borrow();
            if e.destination != id || !e.dest_right() {
                return false;
            }
        }
        for e in &c.right_edges {
            let e = e.borrow();
            if e.destination != id || e.dest_right() {
                return false;
            }
        }
        if let Some(inbound) = self.dest_map.get(&id) {
            for e in inbound {
                if e.borrow().origin != id {
                    return false;
                }
            }
        }
        true
    }

    fn expand_right(&mut self) -> bool {
        let center_id = self.center;
        self.dest = usize::MAX;
        self.last_mutual_dest = -1;
        self.last_mutual_dest_orientation = -1;
        if !self.contigs[center_id].right_forward_branch()
            || self.contigs[center_id].right_edges.is_empty()
        {
            return false;
        }
        let outbound = self.contigs[center_id].right_edges.clone();
        let Some(left_mid_edge) = self.find_representative_mid_edge(&outbound) else {
            return false;
        };
        let mid_id = left_mid_edge.borrow().destination;
        if self.contigs[mid_id].bases.len() < self.min_len {
            return false;
        }
        let mutual_dest = self.find_mutual_dest(&outbound);
        let mutual_dest_orientation = self.last_mutual_dest_orientation;
        let mutual_dest_right = (mutual_dest_orientation & 2) == 2;
        if mutual_dest < 0 || mutual_dest_orientation < 0 {
            return false;
        }
        let dest_id = mutual_dest as usize;
        if self.contigs[dest_id].used || dest_id == center_id {
            return false;
        }
        if mutual_dest_right && !self.contigs[dest_id].right_forward_branch() {
            return false;
        }
        if !mutual_dest_right && !self.contigs[dest_id].left_forward_branch() {
            return false;
        }
        let dest_outbound = {
            let d = &self.contigs[dest_id];
            if mutual_dest_right {
                d.right_edges.clone()
            } else {
                d.left_edges.clone()
            }
        };
        if dest_outbound.is_empty() {
            return false;
        }
        let mutual_dest2 = self.find_mutual_dest(&dest_outbound);
        if mutual_dest2 < 0 || mutual_dest2 as usize != center_id {
            return false;
        }
        let Some(mid_nodes) = self.fetch_mid_nodes(&outbound, true) else {
            return false;
        };
        // `mid_nodes_concur` compares each mid's right-dest against
        // `self.dest`; the reference assigns `dest` before the concurrency
        // check, so set it here or the check sees the `usize::MAX` reset and
        // always fails (silently disabling indirect bubble pops).
        self.dest = dest_id;
        if !self.mid_nodes_concur(&mid_nodes) {
            return false;
        }
        if mutual_dest_right {
            let inbound = self.dest_map.get(&dest_id).cloned();
            self.contigs[dest_id].flip(inbound.as_deref());
        }
        let right_mid_edge = self.contigs[mid_id].get_right_edge(dest_id, Some(1));
        let Some(right_mid_edge) = right_mid_edge else {
            return false;
        };
        self.pop(
            center_id,
            dest_id,
            mid_id,
            left_mid_edge,
            right_mid_edge,
            &mid_nodes,
        )
    }

    fn find_representative_mid_edge(&self, edges: &[EdgeRef]) -> Option<EdgeRef> {
        let mut mid_edge: Option<EdgeRef> = None;
        let mut mid_len = 0usize;
        for e in edges {
            let c = &self.contigs[e.borrow().destination];
            let clen = c.bases.len();
            match &mid_edge {
                None => {
                    mid_edge = Some(e.clone());
                    mid_len = clen;
                }
                Some(me) => {
                    let me_depth = me.borrow().depth;
                    let e_depth = e.borrow().depth;
                    if clen >= self.min_len
                        && (mid_len < self.min_len
                            || e_depth > me_depth
                            || (e_depth == me_depth && clen > mid_len))
                    {
                        mid_edge = Some(e.clone());
                        mid_len = clen;
                    }
                }
            }
        }
        mid_edge
    }

    fn find_mutual_dest(&mut self, edges: &[EdgeRef]) -> i64 {
        self.last_mutual_dest = -2;
        self.last_mutual_dest_orientation = -1;
        for e in edges {
            let mid_id = e.borrow().destination;
            if mid_id == self.center {
                return -1;
            }
            let outbound = {
                let mid = &self.contigs[mid_id];
                if e.borrow().dest_right() {
                    mid.left_edges.clone()
                } else {
                    mid.right_edges.clone()
                }
            };
            for o in &outbound {
                let ob = o.borrow();
                if self.last_mutual_dest < 0 {
                    self.last_mutual_dest = ob.destination as i64;
                    self.last_mutual_dest_orientation = (ob.orientation & 2) as i64;
                } else if self.last_mutual_dest != ob.destination as i64
                    || self.last_mutual_dest_orientation != (ob.orientation & 2) as i64
                {
                    return -1;
                }
            }
        }
        self.last_mutual_dest
    }

    fn fetch_mid_nodes(
        &mut self,
        outbound: &[EdgeRef],
        flip_as_needed: bool,
    ) -> Option<Vec<usize>> {
        let mut mid_nodes: Vec<usize> = Vec::new();
        for e in outbound {
            let mid_id = e.borrow().destination;
            if mid_nodes.contains(&mid_id) {
                return None;
            }
            if self.contigs[mid_id].used {
                return None;
            }
            mid_nodes.push(mid_id);
            if flip_as_needed && e.borrow().dest_right() {
                let inbound = self.dest_map.get(&mid_id).cloned();
                self.contigs[mid_id].flip(inbound.as_deref());
            }
        }
        Some(mid_nodes)
    }

    fn mid_nodes_concur(&self, mid_nodes: &[usize]) -> bool {
        let center_id = self.center;
        let dest_id = self.dest;
        let mut left_dest: i64 = -1;
        let mut right_dest: i64 = -1;
        for &mid_id in mid_nodes {
            let c = &self.contigs[mid_id];
            if c.left_edges.is_empty() || c.right_edges.is_empty() {
                return false;
            }
            for e in &c.left_edges {
                let eb = e.borrow();
                if left_dest < 0 {
                    left_dest = eb.destination as i64;
                } else if left_dest != eb.destination as i64 {
                    return false;
                }
                if eb.origin == eb.destination {
                    return false;
                }
            }
            for e in &c.right_edges {
                let eb = e.borrow();
                if right_dest < 0 {
                    right_dest = eb.destination as i64;
                } else if right_dest != eb.destination as i64 {
                    return false;
                }
                if eb.origin == eb.destination {
                    return false;
                }
            }
            let incoming = self.dest_map.get(&mid_id);
            let Some(incoming) = incoming else {
                return false;
            };
            for e in incoming {
                let origin = e.borrow().origin;
                if origin != center_id && origin != dest_id {
                    return false;
                }
            }
        }
        if left_dest >= 0 && left_dest as usize != center_id {
            return false;
        }
        if right_dest >= 0 && right_dest as usize != dest_id {
            return false;
        }
        left_dest >= 0 && right_dest >= 0
    }

    fn pop(
        &mut self,
        left_id: usize,
        right_id: usize,
        mid_id: usize,
        left_mid_edge: EdgeRef,
        right_mid_edge: EdgeRef,
        mid_nodes: &[usize],
    ) -> bool {
        let k = self.k;
        let original_left_len = self.contigs[left_id].bases.len();
        let mut bb: Vec<u8> = self.contigs[left_id].bases.clone();
        {
            let eb = left_mid_edge.borrow();
            if eb.bases.len() > 1 {
                bb.extend_from_slice(&eb.bases[..eb.bases.len() - 1]);
            }
        }
        {
            let mid = &self.contigs[mid_id];
            let lim = mid.bases.len() - k + 1;
            if k - 1 < lim {
                bb.extend_from_slice(&mid.bases[k - 1..lim]);
            }
        }
        {
            let eb = right_mid_edge.borrow();
            if eb.bases.len() > 1 {
                bb.extend_from_slice(&eb.bases[..eb.bases.len() - 1]);
            }
        }
        bb.extend_from_slice(&self.contigs[right_id].bases);
        self.contigs[left_id].bases = bb;
        self.contigs[left_id].right_edges.clear();
        let right_right = self.contigs[right_id].right_edges.clone();
        if right_right.is_empty() {
            self.contigs[left_id].right_edges = Vec::new();
        } else {
            for e in &right_right {
                e.borrow_mut().origin = left_id;
            }
            self.contigs[left_id].right_edges = right_right;
        }
        self.redirect_edges(right_id, left_id, true);
        let inbound_right = self.dest_map.get(&right_id).cloned();
        set_used(right_id, inbound_right.as_deref(), &mut self.contigs);
        for &c in mid_nodes {
            let inbound = self.dest_map.get(&c).cloned();
            if c == mid_id {
                set_used(c, inbound.as_deref(), &mut self.contigs);
            } else {
                set_associate(c, inbound.as_deref(), &mut self.contigs);
            }
        }
        let right_len = self.contigs[right_id].bases.len();
        let (right_max_cov, right_min_cov, right_code, right_ratio, right_coverage) = {
            let r = &self.contigs[right_id];
            (
                r.max_cov,
                r.min_cov,
                r.right_code,
                r.right_ratio,
                r.coverage,
            )
        };
        let (mid_max_cov, mid_min_cov) = {
            let m = &self.contigs[mid_id];
            (m.max_cov, m.min_cov)
        };
        {
            let left = &mut self.contigs[left_id];
            left.max_cov = left.max_cov.max(right_max_cov).max(mid_max_cov);
            left.min_cov = left.min_cov.min(right_min_cov).min(mid_min_cov);
            left.right_code = right_code;
            left.right_ratio = right_ratio;
            let coverage_sum = left.coverage as f64 * (original_left_len - k + 1) as f64
                + right_coverage as f64 * (right_len - k + 1) as f64;
            left.coverage = (coverage_sum / (left.bases.len() - k + 1) as f64) as f32;
        }
        if self.is_loop(left_id) {
            self.contigs[left_id].left_code = LOOP;
            self.contigs[left_id].right_code = LOOP;
            let inbound = self.dest_map.get(&left_id).cloned();
            remove_all_edges(left_id, inbound.as_deref(), &mut self.contigs);
        }
        self.expansions += 1;
        self.contigs_absorbed += 1 + mid_nodes.len();
        true
    }

    fn remove_dead_edges(&self, c: &mut Contig) {
        c.left_edges.retain(|e| {
            let d = e.borrow().destination;
            let dc = &self.contigs[d];
            !(dc.used || dc.associate)
        });
        c.right_edges.retain(|e| {
            let d = e.borrow().destination;
            let dc = &self.contigs[d];
            !(dc.used || dc.associate)
        });
    }
}

/// `Tadpole.popBubbles`: one bubble-popping pass, then deterministic sort and
/// renumbering.
pub(crate) fn pop_bubbles(contigs: &mut Vec<Contig>, opts: &AssembleOptions) {
    let dest_map = {
        let mut map: HashMap<usize, Vec<EdgeRef>> = HashMap::new();
        for c in contigs.iter() {
            if c.used || c.associate {
                continue;
            }
            for e in &c.left_edges {
                map.entry(e.borrow().destination)
                    .or_default()
                    .push(e.clone());
            }
            for e in &c.right_edges {
                map.entry(e.borrow().destination)
                    .or_default()
                    .push(e.clone());
            }
        }
        map
    };
    let mut bp = BubblePopper {
        contigs: std::mem::take(contigs),
        dest_map,
        k: opts.k,
        min_len: 2 * opts.k - 1,
        center: 0,
        dest: usize::MAX,
        last_mutual_dest: -1,
        last_mutual_dest_orientation: -1,
        expansions: 0,
        contigs_absorbed: 0,
    };
    for i in 0..bp.contigs.len() {
        let c = &bp.contigs[i];
        if !c.used && (c.left_forward_branch() || c.right_forward_branch()) {
            bp.expand(i);
        }
    }
    let dest_map2 = bp.dest_to_edge_map();
    let mut temp: Vec<Contig> = Vec::new();
    for i in 0..bp.contigs.len() {
        if bp.contigs[i].used {
            continue;
        }
        let mut c = bp.contigs[i].clone();
        bp.remove_dead_edges(&mut c);
        temp.push(c);
    }
    temp.sort_by(contig_cmp);
    for (new_id, c) in temp.iter_mut().enumerate() {
        let inbound = dest_map2.get(&c.id).cloned();
        c.renumber(new_id, inbound.as_deref());
    }
    *contigs = temp;
}
