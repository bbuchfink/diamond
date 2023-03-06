extern "C" {
#include <wavefront/wavefront_align.h>
}

#include <string>
#include <fstream>

#include <streambuf>

void do_align(std::string_view const pattern, std::string_view const ref)
{
  // Configure alignment attributes
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.alignment_scope = compute_alignment;
  attributes.distance_metric = gap_affine;
  attributes.affine_penalties.mismatch = 3;
  attributes.affine_penalties.gap_opening = 5;
  attributes.affine_penalties.gap_extension = 1;
  attributes.alignment_form.span = alignment_end2end;
  attributes.alignment_form.pattern_begin_free = 0;
  attributes.alignment_form.pattern_end_free = 0;
  attributes.alignment_form.text_begin_free = 0;
  attributes.alignment_form.text_end_free = 0;
  attributes.memory_mode = wavefront_memory_ultralow;

  attributes.heuristic.strategy = wf_heuristic_wfadaptive;
  attributes.heuristic.min_wavefront_length = 10;
  attributes.heuristic.max_distance_threshold = 50;
  attributes.heuristic.steps_between_cutoffs = 1;
  attributes.memory_mode = wavefront_memory_ultralow;

  wavefront_aligner_t * const wf_aligner = wavefront_aligner_new(&attributes);

  int res = wavefront_align(wf_aligner, pattern.data(), pattern.size(), ref.data(), ref.size());

  cigar_print_pretty(stderr,wf_aligner->cigar,pattern.data(),pattern.size(),ref.data(),ref.size());
  fprintf(stderr,"Alignment Score %d\nResult:%d\n", wf_aligner->cigar->score, res);

  assert(res != -1);
  wavefront_aligner_delete(wf_aligner);
}

int main()
{
    std::ifstream qrystr{"qry.txt"};
    std::string qry{std::istreambuf_iterator<char>(qrystr),
                    std::istreambuf_iterator<char>()};
    std::ifstream refstr{"ref.txt"};
    std::string ref{std::istreambuf_iterator<char>(refstr),
                    std::istreambuf_iterator<char>()};

    do_align(qry, ref);
}
