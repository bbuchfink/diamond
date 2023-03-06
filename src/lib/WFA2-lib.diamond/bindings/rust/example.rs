/// Include the generated bindings into a separate module.
#[allow(non_upper_case_globals)]
#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
#[allow(unused)]
mod wfa {
    include!(concat!(env!("OUT_DIR"), "/bindings_wfa.rs"));
}

/// Compute the affine alignment score between `a` and `b` with the given substitution,
/// gap-open, and gap-extend penalties.
fn linear_score(a: &[i8], b: &[i8], sub: i32, indel: i32) -> i32 {
    unsafe {
        let mut attributes = wfa::wavefront_aligner_attr_default;
        // Do not use a heuristic (enabled by default).
        attributes.heuristic.strategy = wfa::wf_heuristic_strategy_wf_heuristic_none;
        // Only compute the score (not a path).
        attributes.alignment_scope = wfa::alignment_scope_t_compute_score;

        // Set the cost model and parameters.
        attributes.distance_metric = wfa::distance_metric_t_gap_affine;
        attributes.affine_penalties.mismatch = mismatch as i32;
        attributes.affine_penalties.gap_opening = gap_open as i32;
        attributes.affine_penalties.gap_extension = gap_extend as i32;

        // Initialize the aligner object.
        // This should be reused for multiple queries.
        let wf_aligner = wfa::wavefront_aligner_new(&mut attributes);

        // Do the alignment.
        let status = wfa::wavefront_align(
            wf_aligner,
            a.as_ptr(),
            a.len() as i32,
            b.as_ptr(),
            b.len() as i32,
        );
        assert_eq!(status, 0);

        let score = (*wf_aligner).cigar.score;

        // Clean up memory.
        wfa::wavefront_aligner_delete(wf_aligner);
        score
    }
}
